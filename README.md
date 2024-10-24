# ControlledCBO
Numerical implementation for a variant of consensus-based optimization (CBO) algorithms, conrolled-CBO, which introduces a feedback control term to improve convergence towards global minimizers of non-convex functions in multiple dimensions.

The code also includes the implementation of polynomial approximation methods for solving high-dimensional Hamilton-Jacobi-Bellman equations, see more details in following references.

@article{kalise2018polynomial,
  title={Polynomial approximation of high-dimensional Hamilton--Jacobi--Bellman equations and applications to feedback control of semilinear parabolic PDEs},
  author={Kalise, Dante and Kunisch, Karl},
  journal={SIAM Journal on Scientific Computing},
  volume={40},
  number={2},
  pages={A629--A652},
  year={2018},
  publisher={SIAM}
}

@article{dolgov2021tensor,
  title={Tensor decomposition methods for high-dimensional Hamilton--Jacobi--Bellman equations},
  author={Dolgov, Sergey and Kalise, Dante and Kunisch, Karl K},
  journal={SIAM Journal on Scientific Computing},
  volume={43},
  number={3},
  pages={A1625--A1650},
  year={2021},
  publisher={SIAM}
}
