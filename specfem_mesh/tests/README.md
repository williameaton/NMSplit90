# Test suite


A current list of tests in this directory, listed by subdirectory: 

### `integration`: 
-   `test_1_integration`: tests integration of scalar, real functions $f = 1$ and $f = r $ over the volume of the inner core
-   `test_ylm_integration`: tests integration of complex function $f = Y_{lm}^{*} Y_{l'm'}$ over the 
    volume of the inner core. The integration of this function over the solid angle is an orthogonality 
    condition  $$ \int_V Y_{lm}^{*}\, Y_{l'm'} \mathrm{d}\Omega = \delta_{ll'} \delta_{mm'} $$ 
    and so this function drops to the integral of 1 over the radius for the volume integral. This test
    therefore tests both the orthogonality of the spherical harmonics, as well as the integration of a 
    complex scalar function.


### `lagrange`: 
-   None


### `legendre`: 
-   None


### `plm`: 
-   `test_plm`: tests the `Plm` function for a number of l and m values


### `spline`: 
-   Currently a plotting test - needs updating to run from cmd line.


### `ylm`: 
-   `test_ylm`: tests the `ylm_complex` function for a number of l and m values against scipy values. Note that scipy does not include the (-1)^m. 

