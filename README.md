# CO2_cool - Fortran version
A new parameterization of the CO2 15 µm cooling in non-LTE conditions for the Earth's atmosphere.

## Library 

The code is organized in a library (in directory source/modules/) that can be included in a more complex GCM model. 
The subroutine to be called is `CO2_NLTE_COOL`, inside module file `co2cool.f90`. 

## Inputs

The following inputs are required (in order) by `CO2_NLTE_COOL`:
- 6 atmospheric profiles: temperature, pressure, and 4 VMRs of CO2, O, O2, and N2; 
- `lev0`: index of the pressure specified above corresponding to the maximum pressure level (lower boundary) to be considered for calculating the heating rate. Heating rates will be calculate from that pressure up to the minimum pressure specified in the pressure array. E.g, if `p` is specified in the range of 1e3 to 1e-6 hPa (or 1e-6 to 1e3 hPa) and `lev0=index(p(1e0))`, the heating rate will be calculated in the range of 1e0 to 1e-6 hPa.
- `surf_temp`: surface temperature (if set to a negative value, the temperature of the maximum pressure level will be used);
- `hr`: heating rate. This is an input/out array with the same dimension of pressure. It will be calculated ONLY at pressures in the range of `p(lev0)` (max. pressure considered) to the minimum specified pressure (minimum(pressure)). 
- Units: Temperature in K, pressure in hPa, vmrs in mol/mol (not ppm), heating rate in K/day;
- Input profiles can run either from ground to top of the atmosphere (decreasing pressures) or reverse (top to ground with increasing pressures). The pressure grid can be irregular. 
- Important notes:
  1) Pressure levels should include the surface pressure (near 1e3 hPa), even if the 15 µm cooling is to be calculated only at lower pressure levels (higher altitudes), i.e.,   p(lev0) << 1e3 hPa.
  2) If 15 µm cooling shall be calculated only in the non-LTE regime, it is recommended to set up the lower boundary, `p(lev0)`, close to the limit of the LTE/non-LTE transition, e.g. 1 hPa. In this way, more time-consuming calculations in the LTE region will be avoided.

## Output

The output is expressed as heating rate in units of K/day, on the given input grid in the range of p(lev0) to the minimum specified pressure.

## To compile:
- Edit the `Makefile` and change the Fortran compiler to your preferred choice (e.g., gfortran/ifort).
- From this folder, run: `make`

The compilation produces a test program `run_cool` (see below) and a module library file `lib/libco2_cool.a`.

## Test the parametrization on individual profiles
A test program is also provided in `source/main.f90` to test the parametrization on individual profiles. 
### Input file
- The input file `input.dat` is in a fixed format. Do not change the number of commented lines!

- First input at line 9: `n_lev`, `lev0`, `T_surf`.

- Starting from line 12:
    - 6 atmospheric profiles are read (`n_lev` rows are expected). 

### Output file
Output is written in the `output.dat` file.

### To test `main.f90`:
Two input/output files (examples) are provided: input_test.dat and input_test2.dat. The first computes the heating in the full pressure range provided. The second only at pressures smaller than ~1hPa. Follow these steps: 
- `cp input_test.dat input.dat`
- `./run_cool`
- Check that the results in `output.dat` are consistent with `output_test.dat`
- The same procedure can be done for test2.

## To modify the collisional rates:
Collisional rates can be specified in the `constants.f90` module. Default values are as in Funke et al., JQSRT, 2012.
- Rates are defined in the form: z = a * sqrt(T) + b * exp(- g * T^(-1./3)). Name of the coefficients are as follows: 
- for CO2-O: a_zo, b_zo, g_zo  (default: 3.5e-13, 2.32e-9, 76.75)
- for CO2-O2: a_zo2, b_zo2, g_zo2  (default: 7.0e-17, 1.0e-9, 83.8)
- for CO2-N2: a_zn2, b_zn2, g_zn2  (default: 7.0e-17, 6.7e-10, 83.8)
