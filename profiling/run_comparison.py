
from fenics import *
import networkx as nx

import graphnics_ii as gn_ii
import graphnics_mixed_dim as gn_md
import multigraphnics as gn_mp
import numpy as np

set_log_level(50)

import time


# We profile using a honeycomb network with lots of edges 
# NOTE: This is a modified version of the function in graphnics/graphnics/generators.py
def honeycomb(n, m):
    """
    Make honeycomb mesh with inlet

    Args:
        m (int): honeycomb rows
        n (int): honeycomb cols
    """

    # Make hexagonal mesh
    G = nx.hexagonal_lattice_graph(n, m)
    G = nx.convert_node_labels_to_integers(G)

    G.add_node(len(G.nodes))
    G.nodes[len(G.nodes) - 1]["pos"] = [0, -1]
    G.add_edge(len(G.nodes) - 1, 0)

    inlet_node = len(G.nodes) - 1

    # Add outlet
    # We want it positioned at the top right of the mesh
    pos = nx.get_node_attributes(G, "pos")
    all_coords = np.asarray(list(pos.values()))
    all_node_dist_from_origin = np.linalg.norm(all_coords, axis=1)
    furthest_node_ix = np.argmax(all_node_dist_from_origin, axis=0)
    coord_furthest_node = all_coords[furthest_node_ix, :]

    # Add new node a bit above the furthest one
    G.add_node(len(G.nodes))
    G.nodes[len(G.nodes) - 1]["pos"] = coord_furthest_node + np.asarray([0.7, 1])
    G.add_edge(len(G.nodes) - 1, furthest_node_ix)

    # Usually the inlet edge is oriented outwards, we want it inwards
    if (0, inlet_node) in G.edges():
        G.remove_edge(0, inlet_node)
        G.add_edge(inlet_node, 0)

    return G



# Runtime comparisons

np.set_printoptions(precision=1)
from multiphenics import *

if __name__ == "__main__":
    
    fenics_ii_times = []
    fenics_mixed_dim_times = []
    multiphenics_times = []
    
    # clean cache
    import os
    os.system('dijitso clean')
    
    p_bc = Expression('x[0]+x[1]', degree=1)
    
    print('n_edges   fenics_ii    multiphenics  fenics_mixed_dim')
    # repeat N=0 so the cache can be initialized for both
    # versions of the code before we start timing
    for N in [0,0,1,2,3,4,5,6,7]:
        
        
        if N ==0:
            g = nx.Graph()
            g.add_node(0)
            g.add_node(1)
            g.add_node(2)
            g.add_edge(0,1)
            g.add_edge(1,2)
            g.nodes()[0]['pos']=[0,0]
            g.nodes()[1]['pos']=[1,0]
            g.nodes()[2]['pos']=[2,0]
            
        else:
            g = honeycomb(N, N)
            g = nx.convert_node_labels_to_integers(g)
        
        # Run and time fenics_ii version
        G = gn_ii.copy_from_nx_graph(g)
        G.make_mesh(n=1)
        G.make_submeshes()
        
        start_time = time.time()
        model_ii = gn_ii.MixedHydraulicNetwork_ii(G, p_bc=p_bc)
        qp_ii = model_ii.solve()
        
        fenics_ii_times.append(time.time() - start_time)
        
        print(norm(qp_ii[G.num_edges]))
        
        # record number of dofs
        num_dofs = sum([w.dim() for w in model_ii.W])
        
        
        # Run and time multiphenics version
        G = gn_mp.copy_from_nx_graph(g)
        G.make_mesh(n=1)
        G.make_submeshes()
        
        start_time = time.time()
        model = gn_mp.MultiMixedHydraulicNetwork(G, p_bc=p_bc)
        qp = model.solve()
        
        A = block_assemble(model.a_form())
        b = block_assemble(model.L_form())
        
        from IPython import embed; embed()
        
        multiphenics_times.append(time.time() - start_time)
        
        print(norm(qp[G.num_edges]))
        
        File('pressure.pvd') << qp[G.num_edges]
        
        # Run and time fenics_mixed_dim version
        G = gn_md.copy_from_nx_graph(g)
        G.make_mesh(n=1)
        G.make_submeshes()
        
        start_time = time.time()
        #model = gn_md.MixedHydraulicNetwork_mixed_dim(G, p_bc=p_bc)
        #qp = model.solve()
        
        fenics_mixed_dim_times.append(time.time() - start_time)
        

        
        print(f'  {G.num_edges:2.0f}  | {num_dofs} |  {fenics_ii_times[-1]:1.1f}s  |  {multiphenics_times[-1]:1.1f}s  | {fenics_mixed_dim_times[-1]:1.1f}s |')
        
        
        


