
from fenics import *
import networkx as nx

import graphnics_ii as gn_ii
import graphnics_mixed_dim as gn_md

set_log_level(50)

import time

if __name__ == "__main__":
    
    fenics_ii_times = []
    fenics_mixed_dim_times = []
    
    # clean cache
    import os
    os.system('dijitso clean')
    
    p_bc = Expression('x[0]', degree=1)
    
    print('n_edges   fenics_ii   fenics_mixed_dim')
    # repeat N=0 so the cache can be initialized for both
    # versions of the code before we start timing
    for N in [0,0,1,2,3,4,5,6,7]:
        
        
        if N ==0:
            g = nx.Graph()
            g.add_node(0)
            g.add_node(1)
            g.add_edge(0,1)
            g.nodes()[0]['pos']=[0,0]
            g.nodes()[1]['pos']=[1,0]
            
        else:
            g = nx.hexagonal_lattice_graph(N,N)
            g = nx.convert_node_labels_to_integers(g)
        
        # Run and time fenics_ii version
        G = gn_ii.copy_from_nx_graph(g)
        G.make_mesh()
        G.make_submeshes()
        
        start_time = time.time()
        model = gn_ii.MixedHydraulicNetwork(G, p_bc=p_bc)
        qp = model.solve()
        
        fenics_ii_times.append(time.time() - start_time)
        
        
        if N != 0:
            File('p.pvd') << qp[G.num_edges]
        
        # Run and time fenics_mixed_dim version
        G = gn_md.copy_from_nx_graph(g)
        G.make_mesh()
        G.make_submeshes()
        
        start_time = time.time()
        model = gn_md.MixedHydraulicNetwork_mixed_dim(G, p_bc=p_bc)
        qp = model.solve()
        
        fenics_mixed_dim_times.append(time.time() - start_time)
        
        print(f'  {G.num_edges:2.0f}     {fenics_ii_times[-1]:1.1f}s    {fenics_mixed_dim_times[-1]:1.1f}s')
        
        
        