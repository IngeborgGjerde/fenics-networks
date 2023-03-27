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

from graphnics_ii import *

class FenicsGraph_mixed_dim(FenicsGraph_ii):
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

    def make_submeshes(self):
        """
        Generates and stores submeshes for each edge
        """

        # Make and store one submesh for each edge
        for i, (u, v) in enumerate(self.edges):
            self.edges[u, v]["submesh"] = MeshView.create(self.mf, i)
            
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

   
    def ip_jump_lm(self, qs, xi, i):
        '''
        Returns the inner product between the jump of edge fluxes [qs]
        and the lagrange multiplier xi over bifurcation i
        
        Args:
            qs (list): list of edge flux functions
            xi (df.function): lagrange multiplier on bifurcation i
            i (int): bifurcation index
        Returns:
            ip 
        '''

        edge_list = list(self.edges.keys())

        ip = 0
        for e in self.in_edges(i):
            ds_edge = Measure('ds', domain=self.edges[e]['submesh'], subdomain_data=self.edges[e]['vf'])
            edge_ix = edge_list.index(e)
            ip += qs[edge_ix]*xi*ds_edge(BIF_IN)

        for e in self.out_edges(i):
            ds_edge = Measure('ds', domain=self.edges[e]['submesh'], subdomain_data=self.edges[e]['vf'])
            edge_ix = edge_list.index(e)
            ip -= qs[edge_ix]*xi*ds_edge(BIF_OUT)

        return ip

# Network Raviart-Thomas type space, appropriate for the dual mixed hydraulic model
RT = {
    "flux_space": "CG",
    "flux_degree": 1,
    "pressure_space": "DG",
    "pressure_degree": 0,
}


class MixedHydraulicNetwork_mixed_dim:
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
        spaces = P2s + LMs + P1s
        self.W = MixedFunctionSpace(*spaces) 

        self.meshes = submeshes + [G.global_mesh] * (
            G.num_bifurcations + 1
        )  # associated meshes

        

    def solve(self):
        
        G = self.G
        f = self.f
        p_bc = self.p_bc
        
        # Trial and test functions
        vphi = TestFunctions(self.W)
        qp = TrialFunctions(self.W)

        # split out the components
        qs = qp[0:G.num_edges]
        lams = qp[G.num_edges:-1]
        p = qp[-1]

        vs = vphi[0:G.num_edges]
        xis = vphi[G.num_edges:-1]
        phi = vphi[-1]

        
        ## Assemble variational formulation 

        # Initialize a and L to be zero
        dx = Measure('dx', domain=G.global_mesh)
        a = Constant(0)*p*phi*dx
        L = f*phi*dx

        # Assemble edge contributions to a and L
        for i, e in enumerate(G.edges):
            
            msh = G.edges[e]['submesh']
            vf = G.edges[e]['vf']
            
            dx_edge = Measure("dx", domain = msh)
            ds_edge = Measure('ds', domain=msh, subdomain_data=vf)

            # Add variational terms defined on edge
            a += qs[i]*vs[i]*dx_edge        
            a -= p*G.dds(vs[i])*dx_edge
            a += phi*G.dds(qs[i])*dx_edge

            # Add boundary condition for inflow/outflow boundary node
            L += p_bc*vs[i]*ds_edge(BOUN_IN)
            L -= p_bc*vs[i]*ds_edge(BOUN_OUT)

        # Assemble vertex contribution to a, i.e. the bifurcation condition
        for i, b in enumerate(G.bifurcation_ixs):
            a += G.ip_jump_lm(qs, xis[i], b) + G.ip_jump_lm(vs, lams[i], b)
        
        
        # Solve
        qp0 = Function(self.W)
        solve(a==L, qp0)
        
        return qp0


def copy_from_nx_graph(G_nx):
    """
    Return deep copy of nx.Graph as FenicsGraph
    Args:
        G_nx (nx.Graph): graph to be coped
    Returns:
        G (FenicsGraph): fenics type graph with nodes and egdes from G_nx
    """

    G = FenicsGraph_mixed_dim()
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


if __name__ == "__main__":
    
    G = nx.hexagonal_lattice_graph(1,1)
    G = nx.convert_node_labels_to_integers(G)
    G = copy_from_nx_graph(G)
    G.make_mesh()
    G.make_submeshes()
    model = MixedHydraulicNetwork_mixed_dim(G)
    qp = model.solve()
    
    
    
    