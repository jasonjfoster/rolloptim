# rolloptim

## Version 1.0.1

* Added input checks for the number of variables, the `sigma` and `mu` arguments of the `roll_min_var` and `roll_max_utility` functions, and the `xx` and `xy` arguments of the `roll_min_rss` function

* Fixed an issue for solutions on the bounds in the `roll_min_var`, `roll_max_mean`, `roll_max_utility`, and `roll_min_rss` functions

* Fixed an issue for missing values in the `roll_max_utility` and `roll_min_rss` functions

* Added tests for missing values, infeasible targets, univariate `sigma`, `xx`, and `xy` arguments, and input checks

## Version 1.0.0

* Initial CRAN release
