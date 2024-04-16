The code currently is tested with only 1 dataset (v1.1 + support functions with known constants),
If the code throws a crash for any reason please do leave a comment and I will update the validation 
inside of the code.

Other than that it works by taking mesh grids of starting values +- a percentage amount of each starting
value and has some recursive elements (meaning that max array ~900 dimension due to the stack overflow error
and basic settings of most IDEs, which is ~900 variables). These values are then fitted in the first functiton
called get_best_chi_fit and the second function that should be ran is get_best_chi_fit_error. The code finds uncertainty 
within 1 standard deviation of the chi squared value and also outputs the chi squared value as well as the fitted
parameters.

The code has been developed over the course of a week but I have spent a couple months on and off thinking
and trying to start the project so it is possible for there to be some bugs.

I will be testing it throughout my years at university and make any changes and update it as I need to or
notice anything that needs to be upgraded and have time.


For get_best_chi_fit_error, for n>2 it becomes extremely sensitive to resolution, pt (percentage of target chi squared)
and percentage_of_variable where you often need to set see_chi_mesh to true and play around with the 3 settings to 
locate a chi_mesh (as the time complexity of this code is quite huge reasonable resolutions for 3d on my computer
are about 500 for 2mins, 4d is about 120 for >20mins) this means that you have to get quite precise values for pt
and percentage_of_variable for the uncertainties to be located.


My computer ASUS TUF GAMING A15
num dimensions | resolution   | time
2              | <4000        | <1min
3              | 400          | ~1min
4              | 50           | <1min
4              | 120          | ~20min
