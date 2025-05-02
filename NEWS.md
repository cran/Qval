# Qval 1.2.2

-   Added - The `summary` method has been added for the `CDM`, `simData`, and `validation` classes.
-   Added - The `plot` method has been added for the `CDM` and `simData` classes.
-   Added   The `updata` method has been added for the `CDM`, `simData`, and `validation` classes.
-   Change - Further standardized the hypothesis testing procedure for `Wald`.
-   Change - Improved the question names for `Beta`, `Priority`, `PVAF`, and `R2` to be consistent with the original Q matrix.
-   Fixed - Resolved the issue occurring when `method = 'beta'`.
-   Fixed - Improved the code standards for S3 methods

# Qval 1.2.1

-   Fixed - Resolved the issue occurring when `eps = 'logit'`.
-   Fixed - Modified the manual description.

# Qval 1.2.0

-   Added - Separately provided beta metrics.
-   Added - Appropriate S3 methods have been provided for the `CDM`, `simData` and `validation` classes defined in the package to offer better interaction for users.
-   Change - Optimized the efficiency of all the code, making it more in line with the characteristics of the R language.
-   Change - The usage of the `validation` function has been changed when not using the iterative process.
-   Fixed - Resolved the issue occurring when `method = 'Wald'`.
-   Fixed - Resolved the issue occurring when `method = 'beta'`.
-   Fixed - Modified the manual description.

# Qval 1.1.1

-   Added - A search algorithm for 'test.att' has been added for `GDI`, `Hull`, and `beta`.

# Qval 1.1.0

-   Added - Added the URL of the online manual to the `DESCRIPTION` field.
-   Added - A new Q-matrix validation method, the `beta` (β) method, has been added.
-   Added - `Wald` now includes the `SSA` search method and the `item.level` iteration level.
-   Added - Added a PVAF cut-off point for `GDI` predicted by logistic regression (Najera et al., 2019).
-   Changed - Some of the important code is rewritten in C++ to accelerate the Q-matrix validation.
-   Fixed - Resolved the issue occurring when `iter.level = 'item'`.

# Qval 1.0.0

-   Added - Provided a more comprehensive manual description for the `validation` function.
-   Added - Provided a more comprehensive manual description for the `sim.Q` function.
-   Added - Added a function `get.Rmatrix` for calculating the restriction matrix.
-   Added - Added a function `get.priority` for calculating the priority of attribute.
-   Added - Added a flexible function `Wald.test` for calculating the Wald test.
-   Added - Added `SSA` search for the `Hull` method.
-   Added - Added the plotting function `plot.Hull` for the Hull plot.
-   Changed - Rewrote the `Wald` method with the following updates: 1. If the search method is `stepwise` or `forward`, it will call the `Qval` function from the `GDINA` package. 2. If the search method is `PAA`, the search will follow the `PAA`. 3. The information matrix used is the full information matrix, implemented by the internal function `inverse_crossprod` from the `GDINA` package.
-   Changed - Renamed the functions `getQRR`, `getVRR`, `getTPR`, `getTNR`, `getUSR`, and `getOSR` to `zQRR`, `zVRR`, `zTPR`, `zTNR`, `zUSR`, and `zOSR`, respectively.

# Qval 0.1.7

-   Fixed - Corrected the manual description for the `validation` function.

# Qval 0.1.6

-   Fixed - Corrected the manual description for the `validation` function.

# Qval 0.1.5

-   Fixed - Corrected the manual description for the `validation` function.

# Qval 0.1.4

-   Initial release
