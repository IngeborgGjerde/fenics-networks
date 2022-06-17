import networkx as nx
import numpy as np
from fenics import *

'''
The FenicsGraph class constructs fenics meshes from networkx directed graphs.
'''


# Marker tags for inward/outward pointing bifurcation nodes and boundary nodes
BIF_IN = 1
BIF_OUT = 2
BOUN_IN = 3
BOUN_OUT = 4


class FenicsGraph(nx.DiGraph):
    '''
    Make fenics mesh from networkx directed graph

    Attributes:
        global_mesh (df.mesh): mesh for the entire graph
        edges[i].mesh (df.mesh): submesh for edge i
        mf (df.function): 1d meshfunction that maps cell->edge number
        vf (df.function): 0d meshfunction on  edges[i].mesh that stores bifurcation and boundary point data
        global_tangent (df.function): tangent vector for the global mesh, points along edge
    '''


    def __init__(self):
        nx.DiGraph.__init__(self)

        self.global_mesh = None # global mesh for all the edges in the graph

    def make_mesh(self, n=1): 
        '''
        Makes a fenics mesh on the graph with n cells on each edge
        The full mesh is stored in self.global_mesh and a submesh is stored
        for each edge
        '''


        # Store the coordinate dimensions
        geom_dim = len(self.nodes[1]['pos']) 
        self.geom_dim = geom_dim

        # Make list of vertex coordinates and the cells connecting them
        vertex_coords = np.asarray( [self.nodes[v]['pos'] for v in self.nodes()  ] )
        cells_array = np.asarray( [ [u, v] for u,v in self.edges() ] )


        # We first make a mesh with 1 cell per edge
        mesh = Mesh()
        editor = MeshEditor()
        editor.open(mesh, "interval", 1, geom_dim)
        editor.init_vertices(len(vertex_coords))
        editor.init_cells(len(cells_array))

        [editor.add_vertex(i, xi) for i, xi in enumerate(vertex_coords)]
        [editor.add_cell(i, cell.tolist()) for i, cell in enumerate(cells_array)]

        editor.close()


        # Make meshfunction containing edge ixs
        mf = MeshFunction('size_t', mesh, 1)
        mf.array()[:]=range(0,len(self.edges()))
        self.mf = mf

        
        # Refine global mesh until desired resolution
        for i in range(0, n):
            mesh = refine(mesh)
            mf = adapt(mf, mesh)

        # Store refined global mesh and refined mesh function marking branches
        self.global_mesh = mesh
        self.mf = mf


        # Make and store one submesh for each edge
        for i, (u,v) in enumerate(self.edges):
            self.edges[u,v]['submesh']=MeshView.create(self.mf, i)

        # Compute tangent vectors
        self.assign_tangents()



        # Give each edge a Meshfunction that marks the vertex if its a boundary node 
        # or a bifurcation node
        # A bifurcation node is tagged BIF_IN if the edge points into it or BIF_OUT if the edge points out of it
        # A boundary node is tagged BOUN_IN if the edge points into it or BOUN_OUT if the edge points out of it

        # Initialize meshfunction for each edge
        for e in self.edges():
            msh = self.edges[e]['submesh']
            vf = MeshFunction('size_t', msh, 0, 0)
            self.edges[e]['vf'] = vf

        # Make list of bifurcation nodes (connected to three or more edges)
        # and boundary nodes (connected to one edge)
        D = nx.adj_matrix(self).todense()
        conns1 = np.asarray([np.sum(D[i,:]) for i in range(0, D.shape[0])]) # positive connection?
        conns2 = np.asarray([np.sum(D[:,i]) for i in range(0, D.shape[0])]) # negative connection?
        conns = conns1+conns2

        # Store these as global variables
        self.bifurcation_ixs = list(np.where(conns==3)[0])
        self.boundary_ixs = list(np.where(conns==1)[0])


        # Loop through all bifurcation ixs and mark the vfs
        for b in self.bifurcation_ixs:

            for e in self.in_edges(b):
                msh = self.edges[e]['submesh']
                vf = self.edges[e]['vf']

                bif_ix_in_submesh = np.where((msh.coordinates() == self.nodes[b]['pos']).all(axis=1))[0][0]
                vf.array()[bif_ix_in_submesh]=BIF_IN

            for e in self.out_edges(b):
                msh = self.edges[e]['submesh']
                vf = self.edges[e]['vf']

                bif_ix_in_submesh = np.where((msh.coordinates() == self.nodes[b]['pos']).all(axis=1))[0][0]
                vf.array()[bif_ix_in_submesh]=BIF_OUT 

        for b in self.boundary_ixs:
            for e in self.in_edges(b):
                msh = self.edges[e]['submesh']
                vf = self.edges[e]['vf']

                bif_ix_in_submesh = np.where((msh.coordinates() == self.nodes[b]['pos']).all(axis=1))[0][0]
                vf.array()[bif_ix_in_submesh]=BOUN_IN

            for e in self.out_edges(b):
                msh = self.edges[e]['submesh']
                vf = self.edges[e]['vf']

                bif_ix_in_submesh = np.where((msh.coordinates() == self.nodes[b]['pos']).all(axis=1))[0][0]
                vf.array()[bif_ix_in_submesh]=BOUN_OUT
            


    def assign_tangents(self):
        '''
        Assign a tangent vector list to each edge in the graph
        The tangent vector lists are stored 
            * for each edge in G.edges[i]['tangent']
            * as a lookup dictionary in G.tangents
            * as a fenics function in self.global_tangent
        '''
        
        for u,v in self.edges():
            tangent = np.asarray(self.nodes[v]['pos'])-np.asarray(self.nodes[u]['pos'])
            tangent *= 1/np.linalg.norm(tangent)
            self.edges[u,v]['tangent'] = tangent
        self.tangents = list( nx.get_edge_attributes(self, 'tangent').items() )

        tangent = TangentFunction(self, degree=1)
        tangent_i = interpolate(tangent, VectorFunctionSpace(self.global_mesh, 'DG', 0, self.geom_dim) )
        self.global_tangent = tangent_i
        

class GlobalFlux(UserExpression):
    '''
    Evaluated P2 flux on each edge
    '''
    def __init__(self, G, qs, **kwargs):
        '''
        Args:
            G (nx.graph): Network graph
            qs (list): list of fluxes on each edge in the branch
        '''
        
        self.G=G
        self.qs = qs
        super().__init__(**kwargs)

    def eval_cell(self, values, x, cell):
        edge = self.G.mf[cell.index]
        tangent = self.G.tangents[edge][1]
        values[0] = self.qs[edge](x)*tangent[0]
        values[1] = self.qs[edge](x)*tangent[1]
        if self.G.geom_dim == 3: values[2] = self.qs[edge](x)*self.G.tangents[2]

    def value_shape(self):
        return (self.G.geom_dim,)

class TangentFunction(UserExpression):
    '''
    Tangent expression for graph G, which is 
    constructed from G.tangents
    '''
    def __init__(self, G, degree, **kwargs):
        '''
        Args:
            G (nx.graph): Network graph
            degree (int): degree of resulting expression
        '''
        self.G=G
        self.degree=degree
        super().__init__(**kwargs)

    def eval_cell(self, values, x, cell):
        edge = self.G.mf[cell.index]
        values[0] = self.G.tangents[edge][1][0]
        values[1] = self.G.tangents[edge][1][1]
        if self.G.geom_dim==3: values[2] = self.G.tangents[edge][1][2]
    def value_shape(self):
        return (self.G.geom_dim,)



def copy_from_nx_graph(G_nx):
    '''
    Return deep copy of nx.Graph as FenicsGraph
    Args:
        G_nx (nx.Graph): graph to be coped
    Returns:
        G (FenicsGraph): fenics type graph with nodes and egdes from G_nx
    '''
    
    G = FenicsGraph()
    G.graph.update(G_nx.graph)
    G.add_nodes_from((n, d.copy()) for n, d in G_nx._node.items())
    for u,v in G_nx.edges():
        G.add_edge(u,v)
    return G


def test_fenics_graph():
    # Make simple y-bifurcation


    G = FenicsGraph()
    from graph_examples import make_Y_bifurcation
    G = make_Y_bifurcation()
    G.make_mesh()
    mesh = G.global_mesh

    # Check that all the mesh coordinates are also
    # vertex coordinates in the graph
    mesh_c = mesh.coordinates()
    for n in G.nodes:
        vertex_c = G.nodes[n]['pos']
        vertex_ix = np.where((mesh_c == vertex_c).all(axis=1))[0]
        assert len(vertex_ix)==1, 'vertex coordinate is not a mesh coordinate'
        
        




if __name__ == '__main__':
    
    test_fenics_graph()

    
