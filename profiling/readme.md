# Time profiling mixed-dimensional assembly in FEniCS

The assembly of mixed-dimensional graph PDEs in FEniCS has uncovered significant performance issues, depending on the choice of mixed-dimensional library. In particular, we find that using [fenics_ii](https://github.com/MiroK/fenics_ii) yields order-of-magnitude lower assembly times than [fenics-mixed-dim](https://dl.acm.org/doi/abs/10.1145/3471138).

## Model setup
The hydraulic network model solves for a cross-section flux $q$ and pressure $p$ solving

$$\mathcal{R} q + \nabla p = 0 \text{ on } G, $$

$$ \nabla \cdot q = f \text{ on }G, $$

where $\nabla \cdot q$ is the graph divergence.

This model can be discretized in primal-mixed or dual-mixed forms, where the latter yields improved approximations of the flux $q$ but can be computationally expensive. To be more precise, the dual-mixed flux is discretized using the network equivalent of the Raviart-Thomas element,

$$ q_h = \sum_{e \in E} q_e \phi_e, $$

where $\phi_e$ is a piecewise linear function (P1) on the edge meshes $e$. The pressure takes value globally and on graph nodes, and is discretized using discontinuous Galerkin elements over the global mesh and real-functions $\mathbb{R}$ at each bifurcation. 

The assembly of this mixed-dimensional system requires constructing submeshes, generating function spaces on these submeshes, and assembling the local matrices and vectors. 

## Profiling

Comparing implementations of this model in [fenics_ii](https://github.com/MiroK/fenics_ii) and [fenics-mixed-dim](https://dl.acm.org/doi/abs/10.1145/3471138) shows **vast** differences in the assembly time:

| $n_{\text{edges}}$  |  $T_{\text{fenics}_\text{ii}}$  |  $T_{\text{fenics-mixed-dim}}$  
|---|---|---|
| 1  | 0.0s  | 0.0s     |
| 6  | 0.6s  | 47.2    |
| 16  | 1.2s  | 65.1s   |
| 21  | 1.4s  | 92.9s   |
| 26  | 2.0s  | 118.6s  |
| 31  | 2.8s  | 148.3s  |
| 38 | 4.4s  | 217.3s  |

Here, $n_{\text{edges}}$ denotes the number of edges in the network and $T_{\text{fenics}_\text{ii}}$ and $T_{\text{fenics-mixed-dim}}$ are the assembly times in seconds.


The profiling test can be run using the following command:

```python
    python3 run_comparison.py
```