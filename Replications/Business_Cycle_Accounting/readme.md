
#qfqf

  Written by Seungyoon Jeong (jeong299@umn.edu)      
  				               
  08/15/2020					       


- Main.m : a code for solving the model and implement the accounting procedure.

- BCA_findSS.m : Given Sbar and parameters, this function finds Steady-State level of variables.

- BCA_solve.m : With Steady-state level of variables, this function solves equilibrium.
	 	
- eulerres.m : this function returns residuals from Euler equation given inputs of variables.

- BCA_data.m : this function imports data from data.csv file, cleans up them, and return log variables for model inputs.

- BCA_extract.m : this function extract wedges from observables.

- sol_fidexp.m : this function solves equilibrium controlling wedges and expectation.

- eulerres2.m : this function returns residuals from Euler equation given inputs of variables under the control of expectation

- find_Sb.m : this function numerically find Sbar (logzt, Tlt, Txt) to match model's main steady-state variables with long-run mean of data

