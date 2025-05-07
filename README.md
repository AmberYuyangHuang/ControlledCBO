# ControlledCBO
This repository contains the numerical implementation of Controlled Consensus-Based Optimization (Controlled-CBO), a variant of consensus-based optimization algorithms. Controlled-CBO introduces a feedback control term to improve convergence towards global minimizers of non-convex functions in multiple dimensions. See the following preprint:

Huang, Y., Herty, M., Kalise, D., & Kantas, N. (2024). Fast and robust consensus-based optimization via optimal feedback control. arXiv preprint arXiv:2411.03051.

The code also includes polynomial approximation methods for solving high-dimensional Hamilton-Jacobi-Bellman (HJB) equations. These methods are detailed in the following references:

Kalise, D., & Kunisch, K. (2018). Polynomial approximation of high-dimensional Hamilton–Jacobi–Bellman equations and applications to feedback control of semilinear parabolic PDEs. SIAM Journal on Scientific Computing, 40(2), A629–A652.

Dolgov, S., Kalise, D., & Kunisch, K. (2021). Tensor decomposition methods for high-dimensional Hamilton–Jacobi–Bellman equations. SIAM Journal on Scientific Computing, 43(3), A1625–A1650.

For ease of navigation, here is a brief description of the main scripts:

main_hd: the main code for testing the accuracy of Controlled CBO in high-dimensional case for benchmark functions Rastrigin and Ackley function.

illustration: Contains the code to illustrate the trajectories of particle system and evolution of the variance and 2-Wasserstein distance.

approximation_error:  This script plots the value function in a one-dimensional case to illustrate approximation error.

Note: When using a hyperbolic cross basis, the maximum degree of the basis is J = 2^gradmax. For a full basis truncated by total degree, the maximum degree is M = gradmax, where gradmax is the parameter specified in the code.
