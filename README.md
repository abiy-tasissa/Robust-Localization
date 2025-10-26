# robust_node_localization

This is a compressive sensing based approach for localization, where the objective is to determine the positions of target nodes given their distances to $m$ anchor nodes with known positions. We assume that the last anchor node serves as the central node. In our setting, a random subset ($\alpha$) of the target nodes is highly 
corrupted. This means that, for these nodes, a random subset ($k$) of their distance measurements to the anchors is highly corrupted. The corruption model is as follows:

$$\bar{d}_{i,j}^2 = (1 +  U(\beta_{\min},\beta_{max}) d_{i,j}^2$$

where $U(\beta_{\min},\beta_{\max})$ is uniform random variable in the interval $[\beta_{\min},\beta_{\max}]$,   $d_{i,j}^2$ represents the squared distance and $\bar{d}_{i,j}^2$ represents the corrupted distance.

Our goals are to:

1. Identify the highly corrupted nodes
2. Estimate their positions in a robust manner

For details of the experimental setup, please refer to our manuscript. These codes were developed as part of a research project on the robust localization problem by Abiy Tasissa and Waltenegus Dargie. 

## MATLAB files description
`demo_resilient_localization.m`: This compares our proposed algorithms against structured robust PCA for two different experimental setups. It also quantifies
the quality of the position estimates using four evaluation criterion.

`demo_demo_resilient_localization_comparison.m`: This compares our algorithm to the method proposed in the work below:
>  Zaeemzadeh, Alireza, et al. "Robust target localization based on squared range iterative reweighted least squares." 2017 IEEE 14th International Conference on Mobile Ad Hoc and Sensor Systems (MASS). IEEE, 2017.

`robust_alg.m`: This is the main algorithm to estimate the position of a target node given sparse corrupted distances to anchors. 

`generate_near_far_nodes.m`: This generates a set of anchor nodes and target nodes in a rectangular region. A set of the nodes are close to the anchors (referred to be in the near region), and another set of the anchors are far from the anchors (referred to be in the far region). 

`PUpositionSRHybrid.m`: An implementation of a baseline algorithm used for comparison (see reference above).

# Dependencies

* MATLAB: 2022a
* CVX: 2.2
* MOSEK: 10.2.1

## Instructions

To test the algorithm, you can start with the scripts `demo_resilient_localization.

## Citation

If you use our code or find our paper useful and relevant, we would appreciate if you cite our paper. 
> Tasissa, Abiy, and Waltenegus Dargie. "Robust node localization for rough and extreme deployment environments." arXiv preprint arXiv:2507.03856 (2025).
[arXiv link](https://arxiv.org/abs/2507.03856). 

## Feedback

If you have any questions about the code, email <a href="mailto:abiy19@gmail.com">Abiy Tasissa</a>.

