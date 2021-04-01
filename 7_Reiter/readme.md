 # Solving the Krusell-Smith's economy with Reiter's Method
 
  Written by Seungyoon Jeong (jeong299@umn.edu)

  08/01/2020					      

 

- initparams.m : A function for initialize main parameters.

- Main.m : main file to run overall code. 

- Compute_trans.m : A function to compute the entire transition matrix with respect to all possible states.

- discretize_policy.m : A function that discretize policy function over the the finer grid over distribution to construct Q.

- eulerres.m     : functions that returns residual from Euler equation over grids.
  eulerres_gen.m 

- solve_pol.m    : functions that finds an optimal policy function using the policy function iteration.
  solve_pol_gen.m 

- find_dist_ss.m : a function finds steady-state distribution using Q.

- find_ststate.m : a function that find steady-state of Xt, Xss.

- intr.m, cwage.m : a function that finds prices using aggregate variables.

- res_equation.m : a function that takes Xt and Xt_1 and return residual of overall system.

- linearize.m : This function linearize func with respect to two input vectors around the Steady-State level, and return Jacobian matrices.

- simult.m : a function for simulation.

- make_IRF.m : a function that generates IRF.

- makeknotsd.m : a function that generates knots over distribution.




