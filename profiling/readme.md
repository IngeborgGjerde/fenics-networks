# Time profiling of mixed-dimensional assembly in FEniCS

The simulations performed in this paper uncovered significant issues with respect to the assembly time of mixed-dimensional systems in FEniCS. The assembly time was found to depend majorly on the choice of mixed-dimensional library. In particular, we found that using [fenics_ii](https://github.com/MiroK/fenics_ii) yields order-of-magnitude lower assembly times compared [fenics-mixed-dim](https://dl.acm.org/doi/abs/10.1145/3471138).

In these notes we extend [an early benchmarking experiment](https://github.com/IngeborgGjerde/graphnics/blob/afc2a866b32e254a0382fea87cb6222737ce00af/demo/Tree%20profiling.ipynb) to compare these two libraries. 

## Model setup
The hydraulic network model solves for a cross-section flux $q$ and pressure $p$ solving

$$\mathcal{R} q + \nabla p = 0 \text{ on } G, $$

$$ \nabla \cdot q = f \text{ on }G, $$

where $\nabla \cdot q$ is the graph divergence.

This model can be discretized in primal-mixed or dual-mixed forms, where the latter yields improved approximations of the flux $q$ but can be computationally expensive. To be more precise, the dual-mixed flux is discretized using the network equivalent of the Raviart-Thomas element,

$$ q_h = \sum_{e \in E} q_e \phi_e, $$

where $\phi_e$ is a piecewise linear function (P1) on the edge meshes $e$. The pressure takes values globally and on graph nodes, and is discretized using discontinuous Galerkin elements over the global mesh and real-functions $\mathbb{R}$ at each bifurcation. 

The assembly of this mixed-dimensional system requires constructing submeshes, generating function spaces on these submeshes, and assembling the local matrices and vectors. 

## Profiling results

Comparing implementations of this model in [fenics_ii](https://github.com/MiroK/fenics_ii), [multiphenics](https://github.com/multiphenics/multiphenics) and [fenics-mixed-dim](https://dl.acm.org/doi/abs/10.1145/3471138) shows **vast** differences in the assembly time:

| $n_{\text{edges}}$  |  $n_{\text{dofs}}$ | $T_{\text{fenics}_\text{ii}}$  |  $T_{\text{multifenics}}$  |  $T_{\text{fenics-mixed-dim}}$  
|---|---|---|---|---|
|   6  | 36 |  0.2s  |  7.2s | 24.5s |
|  19  | 111 |  0.9s  |  24.6s | 78.0s |
|  38  | 220 |  3.3s  |  45.0s | 205.5s |
|  63  | 363 |  7.6s  |  68.9s | 533.7s |
|  94  | 540 |  16.6s  |  102.7s |   1144.6s |


Here, $n_{\text{edges}}$ denotes the number of edges in the network, $n_{\text{dofs}}$ denotes the total number of finite element degrees of freedom, and $T_{{\text{fenics}}_\text{ii}}$ and $T_{\text{fenics-mixed-dim}}$ are the assembly times in seconds.

The profiling test can be run using the following command:

```python
    python3 run_comparison.py
```

## Discussion
We see that the assembly times are substantial, especially considering that the finite element system is quite small (<1000 degrees of freedom). We see that for (dual mixed) network simulations the runtime is fully dominated by the performance of the mixed-dimensional library.  

We caution that the runtime for multiphenics should be taken with a grain of salt, as we have not spent time optimizing the implementation. Contrarily, we have spent considerable time investigating the fenics_ii and fenics-mixed-dim implementations. Multiphenics was included mostly for the sake of comparison.


## Conclusion
The main takeaway from this comparison is that fenics-mixed-dim underperforms compared to other mixed-dimensional libraries. It is unclear at this point what causes this runtime disparity. **This issue should be investigated further.**

