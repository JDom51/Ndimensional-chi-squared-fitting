The code currently is tested with only 1 dataset and doesnt support functions that have known constants,
to go around the problem of known constants you can multiply the x numpy array by them
Also the code doesn`t have any input validation yet and it is possible or the code to crash.

Other than that it works by taking mesh grids of starting values +- a percentage amount of each starting
value and has some recursive elements (meaning that max array ~900 dimension due to the stack overflow error
and basic settings of most IDEs, which is ~900 variables). These values are then fitted in the first functiton
called get_best_chi_fit and the second function that should be ran is get_best_chi_fit_error.

The code has been developed over the course of a week but I have spent a couple months on and off thinking
and trying to start the project so it is possible for there to be some bugs.

I will be testing it throughout my years at university and make any changes and update it as I need to or
notice anything that needs to be upgraded and have time.