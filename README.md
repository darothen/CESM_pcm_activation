
Analyses for setting up chaos expansion activation parameterization in CESM.

## Listings

\* *denotes work in progress*

1. [`fixed_accom`](./fixed_accom/fixed_accom.md) - Comparison between simulations when the accomodation coefficient, $\alpha_c$ is either fixed or allowed to vary such that $0.1 \leq \alpha_c \leq 1.0$.
2. \*[`diagnostic_Nkn`](./diagnostic_Nkn/diagnostic_Nkn.md) - Given a detailed activation calculation and spectrum, can we back out an $S_\text{max}$ such that the diagnosed droplet concentration (using Köhler theory) will match the kinetically limited estimate?
3. [`aerosol_dists`](./aerosol_dists/) - Plot aerosol size distribution statistics from a reference model simulation.
4. \*[`global_greedy_activation`](./global_greedy_activation/README.md) - Iterative, greedy calculation to analyze which aerosol modes dominate activation dynamics. 
    - calculations performed on legion, then moved locally
