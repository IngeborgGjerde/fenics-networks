import networkx as nx
import numpy as np
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


class FenicsGraph_ii(nx.DiGraph):
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


    def get_mesh(self, n=1):
        """
        Returns a fenics mesh on the graph with 2^n cells on each edge
        Args:
            n (int): number of refinements
            make_edgemeshes (bool): make submesh for each edge 
            
        Returns:
            mesh (df.mesh): the global mesh
            mf (df.MeshFunction): edge marker function 
        """
        
        # Make list of vertex coordinates and the cells connecting them
        vertex_coords = np.asarray([self.nodes[v]["pos"] for v in self.nodes()])
        cells_array = np.asarray([[u, v] for u, v in self.edges()])

        # We first make a mesh with 1 cell per edge
        mesh = Mesh()
        editor = MeshEditor()
        
        some_node_ix = list(self.nodes())[0]
        geom_dim = len(self.nodes()[some_node_ix]["pos"])
        
        editor.open(mesh, "interval", 1, geom_dim)
        editor.init_vertices(len(vertex_coords))
        editor.init_cells(len(cells_array))

        [editor.add_vertex(i, xi) for i, xi in enumerate(vertex_coords)]
        [editor.add_cell(i, cell.tolist()) for i, cell in enumerate(cells_array)]

        editor.close()

        # Make meshfunction containing edge ixs
        mf = MeshFunction("size_t", mesh, 1)
        mf.array()[:] = range(0, len(self.edges()))

        # Refine global mesh until desired resolution
        for i in range(0, n):
            mesh = refine(mesh)
            mf = adapt(mf, mesh)

        return mesh, mf

    
    def make_mesh(self, n=1):
        """
        Generates and stores a fenics mesh on the graph with 2^n cells on each edge
        Args:
            n (int): number of refinements
            make_edgemeshes (bool): make submesh for each edge 
            
        Returns:
            mesh (df.mesh): the global mesh
        """

        # Store the coordinate dimensions
        some_node_ix = list(self.nodes())[0]
        geom_dim = len(self.nodes()[some_node_ix]["pos"])
        self.geom_dim = geom_dim
        self.num_edges = len(self.edges)
        
                
        # Store refined global mesh and refined mesh function marking branches
        mesh, mf = self.get_mesh(n)
        self.global_mesh = mesh
        self.mf = mf
        
        # Store lists with bifurcation and boundary nodes
        self.record_bifurcation_and_boundary_nodes()
        
        # Compute tangent vectors
        self.assign_tangents()
        
    
    def make_submeshes(self):
        """
        Generates and stores submeshes for each edge
        """

        # Make and store one submesh for each edge
        for i, (u, v) in enumerate(self.edges):
            self.edges[u, v]["submesh"] = EmbeddedMesh(self.mf, i)
            
        # Initialize meshfunction for each edge
        for e in self.edges():
            msh = self.edges[e]["submesh"]
            vf = MeshFunction("size_t", msh, 0, 0)
            self.edges[e]["vf"] = vf
        
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

        
    def record_bifurcation_and_boundary_nodes(self):
        '''
        Identify bifurcation nodes (connected to two or more edges)
        and boundary nodes (connected to one edge) and store these as 
        class attributes
       '''
        
        bifurcation_ixs = []
        boundary_ixs = []
        for v in self.nodes():

            num_conn_edges = len(self.in_edges(v)) + len(self.out_edges(v))

            if num_conn_edges == 1:
                boundary_ixs.append(v)
            elif num_conn_edges > 1:
                bifurcation_ixs.append(v)
            elif num_conn_edges == 0:
                print(f"Node {v} in G is lonely (i.e. unconnected)")
        
        # Store these as global variables
        self.bifurcation_ixs = bifurcation_ixs
        self.num_bifurcations = len(bifurcation_ixs)
        self.boundary_ixs = boundary_ixs


    def compute_edge_lengths(self):
        """
        Compute and store the length of each edge
        """

        for e in self.edges():
            v1, v2 = e
            dist = np.linalg.norm(
                np.asarray(self.nodes()[v2]["pos"])
                - np.asarray(self.nodes()[v1]["pos"])
            )
            self.edges()[e]["length"] = dist
            
    def compute_vertex_degrees(self):
        """
        Compute and store the min and max weighted vertex degrees
        """
        
        # Check that edge lengths have been computed
        try:
            e = list(self.edges())[0]
            length = self.edges()[e]["length"]
        except KeyError:
            # if not we compute them now
            self.compute_edge_lengths()
        
        # vertex degree = sum_i L_i/2 for all edges i connected to the vertex
        for v in self.nodes():
            l_v = 0
            for e in self.in_edges(v):
                l_v += self.edges()[e]["length"]
                
            for e in self.out_edges(v):
                l_v += self.edges()[e]["length"]

            self.nodes()[v]['degree'] = l_v/2
            
        degrees = nx.get_node_attributes(self, 'degree')
        self.degree_min = min(degrees.values())
        self.degree_max = max(degrees.values())

            
        
    def assign_tangents(self):
        """
        Assign a tangent vector list to each edge in the graph
        The tangent vector lists are stored
            * for each edge in G.edges[i]['tangent']
            * as a lookup dictionary in G.tangents
            * as a fenics function in self.global_tangent
        """

        for u, v in self.edges():
            tangent = np.asarray(self.nodes[v]["pos"]) - np.asarray(
                self.nodes[u]["pos"], dtype=np.float64
            )
            tangent_norm = np.linalg.norm(tangent)
            tangent_norm_inv = 1.0 / tangent_norm
            tangent *= tangent_norm_inv
            self.edges[u, v]["tangent"] = tangent
        self.tangents = list(nx.get_edge_attributes(self, "tangent").items())

        tangent = TangentFunction(self, degree=1)
        tangent_i = interpolate(
            tangent, VectorFunctionSpace(self.global_mesh, "DG", 0, self.geom_dim)
        )
        self.global_tangent = tangent_i

    def dds(self, f):
        """
        function for derivative df/ds along graph
        """
        return dot(grad(f), self.global_tangent)

    def dds_i(self, f, i):
        """
        function for derivative df/ds along graph on branch i
        """
        tangent = self.tangents[i][1]
        return dot(grad(f), Constant(tangent))

    def get_num_inlets_outlets(self):
        num_inlets, num_outlets = 0, 0

        for e in self.edges():
            vf_vals = self.edges[e]["vf"].array()

            num_inlets += len(list(np.where(vf_vals == BOUN_IN)[0]))
            num_outlets += len(list(np.where(vf_vals == BOUN_OUT)[0]))

        return num_inlets, num_outlets


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
                msh = self.edges[e]["submesh"]
                vf = self.edges[e]["vf"]

                node_ix_in_submesh = np.where(
                    (msh.coordinates() == self.nodes[n]["pos"]).all(axis=1)
                )[0]
                if len(node_ix_in_submesh) > 0:
                    vf.array()[node_ix_in_submesh[0]] = tag_in


            for e in self.out_edges(n):
                msh = self.edges[e]["submesh"]
                vf = self.edges[e]["vf"]

                node_ix_in_submesh = np.where(
                    (msh.coordinates() == self.nodes[n]["pos"]).all(axis=1)
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

    G = FenicsGraph_ii()
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


class MixedHydraulicNetwork:
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
        submeshes = list(nx.get_edge_attributes(G, "submesh").values())

        P2s = [
            FunctionSpace(msh, space["flux_space"], space["flux_degree"])
            for msh in submeshes
        ]
        P1s = [
            FunctionSpace(
                G.global_mesh, space["pressure_space"], space["pressure_degree"]
            )
        ]
        LMs = [FunctionSpace(G.global_mesh, "R", 0) for b in G.bifurcation_ixs]

        ### Function spaces
        W = P2s + P1s + LMs
        self.W = W

        self.meshes = submeshes + [G.global_mesh] * (
            G.num_bifurcations + 1
        )  # associated meshes

        # Trial and test functions
        self.qp = list(map(TrialFunction, W))
        self.vphi = list(map(TestFunction, W))



    def a_form(self, a=None):
        """
        Add edge contributions to the bilinear form
        """

        if not a:
            a = self.init_a_form()

        qp, vphi = self.qp, self.vphi
        G = self.G
        

        # split out the components
        n_edges = G.num_edges
        
        qs, vs = qp[:n_edges], vphi[:n_edges]
        p, phi = qp[n_edges], vphi[n_edges]
        lams, xis = qp[n_edges+1:], vphi[n_edges+1:]

        # get submeshes and restriction of p to edge
        submeshes = list(nx.get_edge_attributes(G, "submesh").values())
        ps = [Restriction(p, msh) for msh in submeshes]
        phis = [Restriction(phi, msh) for msh in submeshes]

        # add edge terms
        for i, e in enumerate(G.edges):
            dx_edge = Measure("dx", domain=G.edges[e]["submesh"])
            
            try: 
                res = self.G.edges()[e]["Res"]
            except KeyError:
                res = Constant(1)

            a[i][i] += res * qs[i] * vs[i] * dx_edge
            
            a[n_edges][i] += +G.dds_i(qs[i], i) * phis[i] * dx_edge
            a[i][n_edges] += -ps[i] * G.dds_i(vs[i], i) * dx_edge
            
            
        # add vertex terms
        edge_list = list(G.edges.keys())
        
        for b_ix, b in enumerate(G.bifurcation_ixs):

            # make info dict of edges connected to this bifurcation
            conn_edges = {
                **{e: (1, BIF_IN, edge_list.index(e)) for e in G.in_edges(b)},
                **{e: (-1, BIF_OUT, edge_list.index(e)) for e in G.out_edges(b)},
            }

            for e in conn_edges:
                ds_edge = Measure(
                    "ds", domain=G.edges[e]["submesh"], subdomain_data=G.edges[e]["vf"]
                )

                sign, tag, e_ix = conn_edges[e]

                a[e_ix][n_edges + 1 + b_ix] += (
                    sign * vs[e_ix] * lams[b_ix] * ds_edge(tag)
                )
                a[n_edges + 1 + b_ix][e_ix] += (
                    sign * qs[e_ix] * xis[b_ix] * ds_edge(tag)
                )

        return a

    def init_a_form(self):
        """
        Init a
        """

        ## Init a as list of lists
        a = [[0 for i in range(0, len(self.qp))] for j in range(0, len(self.qp))]

        # Init zero diagonal elements (for shape info)
        for i, msh in enumerate(self.meshes):
            a[i][i] += (
                Constant(0) * self.qp[i] * self.vphi[i] * Measure("dx", domain=msh)
            )

        return a


    def init_L_form(self):
        """
        Init L
        """

        L = [0 for i in range(0, len(self.vphi))]

        # Init zero diagonal elements (for shape info)
        for i, msh in enumerate(self.meshes):
            dx = Measure("dx", domain=msh)
            L[i] += Constant(0) * self.vphi[i] * dx
            

        return L
    

    def L_form(self):
        """
        The right-hand side linear form
        """

        L = self.init_L_form()

        vphi = self.vphi

        # split out the components
        n_edges = self.G.num_edges
        vs, phi, xis = vphi[0:n_edges], vphi[n_edges], vphi[n_edges + 1 :]

        submeshes = list(nx.get_edge_attributes(self.G, "submesh").values())
        phis = [Restriction(phi, msh) for msh in submeshes]
        
        fs = [project(self.f, FunctionSpace(msh, 'CG', 1)) for msh in submeshes]
        gs = [project(self.g, FunctionSpace(msh, 'CG', 1)) for msh in submeshes]
        p_bcs = [project(self.p_bc, FunctionSpace(msh, 'CG', 1)) for msh in submeshes]
        
        # Assemble edge contributions to a and L
        for i, e in enumerate(self.G.edges):
            ds_edge = Measure("ds", domain=self.G.edges[e]["submesh"], subdomain_data=self.G.edges[e]["vf"])
            dx_edge = Measure("dx", domain=self.G.edges[e]["submesh"])

            L[i] -= p_bcs[i]*vs[i]*ds_edge(BOUN_OUT) - self.p_bc*vs[i]*ds_edge(BOUN_IN)
            L[i] += gs[i] * vs[i] * dx_edge
            
            L[n_edges] += fs[i] * phis[i] * dx_edge
            
        for i in range(0, len(self.G.bifurcation_ixs)):
            L[n_edges + 1 + i] += Constant(0) * xis[i] * dx

        return L
    
    def get_bc(self):
        return [[], []]
    
    def solve(self):
        
        W = self.W
        a = self.a_form()
        L = self.L_form()

        A, b = map(ii_assemble, (a, L))
        A, b = map(ii_convert, (A, b))

        qp = ii_Function(W)
        solver = LUSolver(A, "mumps")
        solver.solve(qp.vector(), b)
        
        return qp
    
    

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
    model = MixedHydraulicNetwork(G)
    qp = model.solve()
    
    
    
    