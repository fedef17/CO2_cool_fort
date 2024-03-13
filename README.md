# CO2_cool - Fortran version
New parametrization of CO2 heating rate in 15um band in non-LTE conditions.

## Stand-alone version
A main program is also provided in `source/main.f90` to test the parametrization on individual profiles.

### To compile:
- Open the Makefile and change the Fortran compiler to your preferred choice (gfotran/ifort).
- from this folder, run: `make`

### To test:
- `cp input_test.dat input.dat`
- `./run_cool`
- Check that the result in `output.dat` is consistent with `output_test.dat`

### Inputs
- The input file `input.dat` is in a fixed format. Do not change the number of commented lines!
- First input at line 9: Number of levels in the atmosphere (n_lev), bottom level (should be 1 or n_lev), surface temperature
- Starting from line 12, 6 atmospheric profiles are read: pressure, temperature, VMRs of CO2, O, O2, N2. 
Profiles can either go from ground to top or reverse (determined by which is the bottom level), temperature in K, pressure in hPa, vmrs in absolute fraction (not ppm). 

- Profiles are needed up to 0.001 hPa for the parametrization to work. Bottom level should be at the ground for an accurate calculation.


## Library version

The code is organized in a library (in directory source/modules/) that can be included in a more sophisticated model. 
The main subroutine is CO2_NLTE_COOL, inside co2cool.f90. 

Collisional rates can be specified in the constants.f90 module. 

### To modify the collisional rates:
- Rates are defined in the form: z = a*np.sqrt(T) + b * np.exp(-g * T**(-1./3)). Name of the coefficients are as follows: 
    - for CO2-O: a_zo, b_zo, g_zo  (default: 3.5e-13, 2.32e-9, 76.75)
    - for CO2-O2: a_zo2, b_zo2, g_zo2  (default: 7.0e-17, 1.0e-9, 83.8)
    - for CO2-N2: a_zn2, b_zn2, g_zn2  (default: 7.0e-17, 6.7e-10, 83.8)
