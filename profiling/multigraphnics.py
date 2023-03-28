import networkx as nx
from multiphenics import *
from fenics import *
from xii import *

'''
Copyright (C) 2022-2023 by Ingeborg Gjerde
This file is a part of the graphnics project (https://arxiv.org/abs/2212.02916)
You can freely redistribute it and/or modify it under the terms of the GNU General Public License, version 3.0, provided that the above copyright notice is kept intact and that the source code is made available under an open-source license.
'''


"""
The Graphnics class constructs fenics meshes from networkx directed graphs.
"""


# Marker tags for inward/outward pointing bifurcation nodes and boundary nodes
BIF_IN = 1
BIF_OUT = 2
BOUN_IN = 3
BOUN_OUT = 4

from graphnics_ii import *

class Branch(SubDomain):
    def __init__(self, branch_ix, mesh, G, **kwargs):
        self.bix = branch_ix
        self.mesh = mesh
        self.G = G
        super().__init__(**kwargs)

    def inside(self, x, on_boundary):
        
        # FIXME: This is wrong
        tree = BoundingBoxTree()
        tree.build(self.mesh)
        cell = tree.compute_first_entity_collision(Point(x))
        print(x, cell)
        if cell>len(self.G.mf.array()):
            return False
        else:
            in_branch = self.bix == self.G.mf[cell]
        return in_branch



class AllBranches(SubDomain):
    def inside(self, x, on_boundary):
        return True

class MultiFenicsGraph(FenicsGraph_ii):
    """
    Make fenics mesh from networkx directed graph
    Attributes:
        global_mesh (df.mesh): mesh for the entire graph
        edges[i].mesh (df.mesh): submesh for edge i
        global_tangent (df.function): tangent vector for the global mesh, points along edge
        dds (function): derivative d/ds along edge
        mf (df.function): 1d meshfunction that maps cell->edge number
        vf (df.function): 0d meshfunction on  edges[i].mesh that stores bifurcation and boundary point data
        num_edges (int): number of edges in graph
    """

    def __init__(self):
        nx.DiGraph.__init__(self)


    def make_submeshes(self):
        """
        Generates and stores submeshes for each edge
        """

        # Make and store one submesh for each edge
        
        mf_test = MeshFunction("size_t", self.global_mesh, 1, 0)
        
        for i, (u, v) in enumerate(self.edges):            
            branch = Branch(i, self.global_mesh, self)
            self.edges[u, v]["restriction"] = MeshRestriction(self.global_mesh, branch)
            
            branch.mark(mf_test, i)
            
        File('branch.pvd') << mf_test
        
            
        # Initialize meshfunction marking bifurcation nodes
        for e in self.edges():
            self.edges()[e]['vf'] = MeshFunction("size_t", self.global_mesh, 0, 0)
        
        # Give each edge a Meshfunction that marks the vertex if its a boundary node
        # or a bifurcation node
        # A bifurcation node is tagged BIF_IN if the edge points into it or BIF_OUT if the edge points out of it
        # A boundary node is tagged BOUN_IN if the edge points into it or BOUN_OUT if the edge points out of it

        # Mark the meshfunction belonging to each edge with
        # - BIF_IN/BIF_OUT if the vertex belongs to an edge that goes in/out of the bif node
        self.mark_edge_vfs(self.bifurcation_ixs, BIF_IN, BIF_OUT)
        # - BOUN_OUT/BOUN_IN if the vertex is an inlet or outlet point
        self.mark_edge_vfs(self.boundary_ixs, BOUN_OUT, BOUN_IN)
        # Note: We need to reverse our notation here so that BOUN_OUT is applied to
        # edges going into a boundary node (which is actually an outlet)

    def mark_edge_vfs(self, node_ixs, tag_in, tag_out):
        """ Mark the meshfunction belonging to the edges submesh 
        with tag_in if the edge goes into a node in node_ixs
        and tag_out if the edge goes out of a node in node_ixs
        Args:
            node_ixs (list): _description_
            edges (func): _description_
            tag (int): corresponding tag 
        """
        for n in node_ixs:

            for e in self.in_edges(n):
                vf = self.edges[e]["vf"]

                node_ix_in_submesh = np.where(
                    (self.global_mesh.coordinates() == self.nodes[n]["pos"]).all(axis=1)
                )[0]
                if len(node_ix_in_submesh) > 0:
                    vf.array()[node_ix_in_submesh[0]] = tag_in


            for e in self.out_edges(n):
                vf = self.edges[e]["vf"]

                node_ix_in_submesh = np.where(
                    (self.global_mesh.coordinates() == self.nodes[n]["pos"]).all(axis=1)
                )[0]
                if len(node_ix_in_submesh) > 0:
                    vf.array()[node_ix_in_submesh[0]] = tag_out


        



def copy_from_nx_graph(G_nx):
    """
    Return deep copy of nx.Graph as FenicsGraph
    Args:
        G_nx (nx.Graph): graph to be coped
    Returns:
        G (FenicsGraph): fenics type graph with nodes and egdes from G_nx
    """

    G = MultiFenicsGraph()
    G.graph.update(G_nx.graph)
    
    # copy nodes and edges
    G.add_nodes_from((n, d.copy()) for n, d in G_nx._node.items())
    for u, v in G_nx.edges():
        G.add_edge(u, v)
    
    # copy node attributes
    for n in G.nodes():
        node_dict = G_nx.nodes()[n]
        for key in node_dict:
            G.nodes()[n][key] = G_nx.nodes()[n][key]
            
    # copy edge attributes
    for e in G.edges():
        edge_dict = G_nx.edges()[e]
        for key in edge_dict:
            G.edges()[e][key] = G_nx.edges()[e][key]
        

    return G


# Network Raviart-Thomas type space, appropriate for the dual mixed hydraulic model
RT = {
    "flux_space": "CG",
    "flux_degree": 1,
    "pressure_space": "DG",
    "pressure_degree": 0,
}


class MultiMixedHydraulicNetwork(MixedHydraulicNetwork_ii):
    """
    Bilinear forms a and L for the dual mixed form of the hydraulic equations
            Res*q + d/ds p = g
            d/ds q = f
    on graph G, with bifurcation conditions q_in = q_out and continuous
    normal stress
    Args:
        G (FenicsGraph): Network domain
        Res (dict): dictionary with edge->resistance
        f (df.expr): source term
        p_bc (df.expr): pressure boundary condition
    """

    def __init__(self, G, f=Constant(0), g = Constant(0), p_bc=Constant(0), space=RT):
        """
        Set up function spaces and store model parameters f and ns
        """

        # Graph on which the model lives
        self.G = G

        # Model parameters

        self.f = f
        self.g = g
        self.p_bc = p_bc

        # Setup function spaces:
        # - flux space on each segment, ordered by the edge list
        # - pressure space on the full mesh
        # - real space on each bifurcation

        restrictions = nx.get_edge_attributes(G, 'restriction').values()
        all_branch = AllBranches()
        restrictions = list(restrictions) + [ MeshRestriction(G.global_mesh, all_branch) for j in range(G.num_bifurcations+1)]
        
        P2 = FunctionSpace(G.global_mesh, 'CG', 1)
        P1 = FunctionSpace(G.global_mesh, 'DG', 0)  
        LM = FunctionSpace(G.global_mesh, 'R', 0)
        spaces = [P2 for i in range(G.num_edges)] + [P1] + [LM for j in range(G.num_bifurcations)]
        
        W = BlockFunctionSpace(spaces, restrict=restrictions)

        self.W = W

        # Trial and test functions
        self.qp = BlockTrialFunction(W)
        self.vphi =BlockTestFunction(W)



    def a_form(self):
        """
        Add edge contributions to the bilinear form
        """

        ## Init a as list of lists
        a = [[0 for i in range(0, len(self.qp))] for j in range(0, len(self.qp))]

        qp = block_split(self.qp)
        vphi = block_split(self.vphi)
        
        G = self.G
        
        # split out the components
        n_edges = G.num_edges
        
        qs, vs = qp[:n_edges], vphi[:n_edges]
        p, phi = qp[n_edges], vphi[n_edges]
        lams, xis = qp[n_edges+1:], vphi[n_edges+1:]

        dx = Measure("dx")(subdomain_data=G.mf)
            
        # add edge terms
        for i, e in enumerate(G.edges):
            
            a[i][i] += qs[i] * vs[i] * dx(i)
            a[n_edges][i] += +G.dds_i(qs[i], i) * phi * dx(i)
            a[i][n_edges] += -p * G.dds_i(vs[i], i) * dx(i)
            
            
        # add vertex terms
        edge_list = list(G.edges.keys())
        
        for b_ix, b in enumerate(G.bifurcation_ixs):

            # make info dict of edges connected to this bifurcation
            conn_edges = {
                **{e: (1, BIF_IN, edge_list.index(e)) for e in G.in_edges(b)},
                **{e: (-1, BIF_OUT, edge_list.index(e)) for e in G.out_edges(b)},
            }

            for e in conn_edges:
                vf = G.edges()[e]['vf']
                dS = Measure("dS")(subdomain_data=vf)
                
                sign, tag, e_ix = conn_edges[e]

                a[e_ix][n_edges + 1 + b_ix] += (
                    sign * vs[e_ix]('+') * lams[b_ix]('+') * dS(tag)
                )
                a[e_ix][n_edges + 1 + b_ix] += (
                    sign * vs[e_ix]('+') * lams[b_ix]('+') * dS(tag)
                )
                
                a[n_edges + 1 + b_ix][e_ix] = (
                    sign * qs[e_ix]('-') * xis[b_ix]('+') * dS(tag)
                )
                a[n_edges + 1 + b_ix][e_ix] = (
                    sign * qs[e_ix]('-') * xis[b_ix]('+') * dS(tag)
                )
                
                a[n_edges + 1 + b_ix][e_ix] = (
                    sign * qs[e_ix]('+') * xis[b_ix]('-') * dS(tag)
                )
                a[n_edges + 1 + b_ix][e_ix] = (
                    sign * qs[e_ix]('+') * xis[b_ix]('-') * dS(tag)
                )
                
                
                a[n_edges + 1 + b_ix][e_ix] = (
                    sign * qs[e_ix]('-') * xis[b_ix]('-') * dS(tag)
                )
                a[n_edges + 1 + b_ix][e_ix] = (
                    sign * qs[e_ix]('-') * xis[b_ix]('-') * dS(tag)
                )

        return a

    
    def L_form(self):
        """
        The right-hand side linear form
        """

        L = [0 for i in range(len(self.vphi))]

        vphi = block_split(self.vphi)

        # split out the components
        n_edges = self.G.num_edges
        vs, phi, xis = vphi[0:n_edges], vphi[n_edges], vphi[n_edges + 1 :]

        dx = Measure("dx")(subdomain_data=self.G.mf)

        # Assemble edge contributions to a and L
        for i, e in enumerate(self.G.edges):
            ds_edge = Measure("ds")(subdomain_data=self.G.edges[e]["vf"])
            
            L[i] -= self.p_bc*vs[i]*ds_edge(BOUN_OUT) - self.p_bc*vs[i]*ds_edge(BOUN_IN)
            L[i] += self.g * vs[i] * dx(i)
            
        L[n_edges] += self.f * phi * dx
            
        for i in range(0, len(self.G.bifurcation_ixs)):
            L[n_edges + 1 + i] += Constant(0) * xis[i] * dx

        return L
    
    def solve(self):
        
        W = self.W
        a = self.a_form()
        L = self.L_form()

        A = block_assemble(a)
        b = block_assemble(L)
        
        U = BlockFunction(W)
        block_solve(A, U.block_vector(), b)

        return U
    
    

class TangentFunction(UserExpression):
    """
    Tangent expression for graph G, which is
    constructed from G.tangents
    """

    def __init__(self, G, degree, **kwargs):
        """
        Args:
            G (nx.graph): Network graph
            degree (int): degree of resulting expression
        """
        self.G = G
        self.degree = degree
        super().__init__(**kwargs)

    def eval_cell(self, values, x, cell):
        edge = self.G.mf[cell.index]
        values[0] = self.G.tangents[edge][1][0]
        values[1] = self.G.tangents[edge][1][1]
        if self.G.geom_dim == 3:
            values[2] = self.G.tangents[edge][1][2]

    def value_shape(self):
        return (self.G.geom_dim,)


if __name__ == "__main__":
    
    G = nx.hexagonal_lattice_graph(2,2)
    G = nx.convert_node_labels_to_integers(G)
    G = copy_from_nx_graph(G)
    G.make_mesh()
    G.make_submeshes()
    model = MultiMixedHydraulicNetwork(G)
    qp = model.solve()
    
    
    
    