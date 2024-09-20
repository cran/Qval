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
