# CO2_cool - Fortran version
A new parameterization of the CO2 15 µm cooling in non-LTE conditions.

## Library 

The code is organized in a library (in directory source/modules/) that can be included in a more sophisticated model. 
The main subroutine is CO2_NLTE_COOL, inside module file co2cool.f90. 

## Inputs

The following inputs are required by CO2_NLTE_COOL:
- n_lev: Number of levels in the atmosphere (n_lev);
- ilev0: index of the lower atmospheric level (maximum  pressure level) to be considered. Parametrization will only be activated above the selected level (at lower pressures);
- T_surf: surface temperature (if set to a negative value, the temperature of the first level from the surface is used);
- 6 atmospheric profiles: pressure, temperature, VMRs of CO2, O, O2, N2 
- temperature in K, pressure in hPa, vmrs in mol/mol (not ppm);
- the whole vertical range is needed, from the surface.
- Input profiles can either go from ground to top or reverse;

## Output

The output is expressed as heating rate in units of K/day, on the same input grid.

## To compile:
- Open the Makefile and change the Fortran compiler to your preferred choice (gfortran/ifort).
- From this folder, run: `make`

The compilation produces a test program `run_cool` and a module library file `lib/libco2_cool.a`.

## Test program
A main program is also provided in `source/main.f90` to test the parametrization on individual profiles.

### To test:
- `cp input_test.dat input.dat`
- `./run_cool`
- Check that the result in `output.dat` is consistent with `output_test.dat`

### Input file
- The input file `input.dat` is in a fixed format. Do not change the number of commented lines!

- First input at line 9: n_lev, ilev0, T_surf

- Starting from line 12:
    - 6 atmospheric profiles are read (n_lev rows are expected)

### Output file

Output is written in the `output.dat` file.

## To modify the collisional rates:
Collisional rates can be specified in the constants.f90 module. 

- Rates are defined in the form: z = a*np.sqrt(T) + b * np.exp(-g * T**(-1./3)). Name of the coefficients are as follows: 
- for CO2-O: a_zo, b_zo, g_zo  (default: 3.5e-13, 2.32e-9, 76.75)
- for CO2-O2: a_zo2, b_zo2, g_zo2  (default: 7.0e-17, 1.0e-9, 83.8)
- for CO2-N2: a_zn2, b_zn2, g_zn2  (default: 7.0e-17, 6.7e-10, 83.8)