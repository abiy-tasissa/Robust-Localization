# robust_node_localization

The sets of scripts here consider the problem of localization, where the objective is to determine the positions of target nodes given their distances to m 
anchor nodes with known positions. We assume that the last anchor node serves as the central node. In our setting, a random subset (alpha) of the target nodes is highly 
corrupted. This means that, for these nodes, a random subset (k) of their distance measurements to the anchors is highly corrupted. The corruption model is as follows:
(1 +  U(beta_min,beta_max)) * d_{i,j}^2
where U(beta_min,beta_max) is uniform random variable in the interval [beta_min,beta_max], and  d_{i,j}^2 represents the squared distance. 

Our goal is to:
(a) Identify the highly corrupted nodes.
(b) Estimate their positions in a robust manner.

For details of the experimental setup, please refer to our manuscript. 

This set of codes were written by Abiy Tasissa. 

## MATLAB files description
`demo_resilient_localization.m`: This compares our proposed algorithms against structured robust PCA for two different experimental setups. It also quantifies
the quality of the position estimates using four evaluation criterion.

`demo_demo_resilient_localization_comparison.m`: This compares our algorithm to the method proposed in the work below:
>  Zaeemzadeh, Alireza, et al. "Robust target localization based on squared range iterative reweighted least squares." 2017 IEEE 14th International Conference on Mobile Ad Hoc and Sensor Systems (MASS). IEEE, 2017.

`robust_alg.m`: This is the main algorithm to estimate the position of a target node given sparse corrupted distances to anchors. 

# Dependencies

* MATLAB: 2022a
* CVX: 2.2
* MOSEK: 10.2.1

## Instructions

To test the algorithm, you can start with either of the scripts `demo_resilient_localization` 


## Feedback

If you have any questions about the code, email <a href="mailto:abiy19@gmail.com">Abiy Tasissa</a>.

