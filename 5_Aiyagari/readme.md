# Solving a model economcy by Aiyagari, "Uninsured Idiosyncratic Risk and Aggregate Saving (1994)"

  Written by Seungyoon Jeong (jeong299@umn.edu)      

  07/23/2020					       

- EDG_main.m : a code to solve the household problem with Endogenoue Grid Method (EGM) and find a equilibrium price.

- EDG_grid.m : a code for solving the model with Endogenous Grid Method and find stationary distribution with Monte-Carlo
    (It has a dependency to a function search.m and moveNext.m)

- search.m : a function that for given asset level on asset grid, find adjacent grids from endogenous grid. It returns each adjacent grid and indices.

- moveNext.m : a function that construct consumption policy in next period, using budget constraint when it binds, otherwise using interpolation

- collocation_8.m : a code for solving the model with Finite Element Method with 8 collocation points and save the result.
    (It has a dependency to a function expected.m) 
    
- expected.m : a function that evaluate RHS of Euler equation.

- EDGvsFEM.m : a code that comparing results from FEM and EGM
