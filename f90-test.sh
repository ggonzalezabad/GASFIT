# Compilation for development:
#gfortran -g -c -fbacktrace -fcheck=all -Wall gasfit_module.f90 main.f90
#gfortran -g -fbacktrace -fcheck=all -Wall -o gasbeta.exe main.o gasfit_module.o
# Compilation for production
gfortran -O3 -c gasfit_module.f90 main.f90
gfortran -O3 -o gasbeta.exe main.o gasfit_module.o

# Running the program
./gasbeta.exe << EOF
gasfit.inp
gasfit_ace_io_output.txt
gasfit_ace_geolocation_for_fitting_results.txt
gasfit_ace_solar_fitting_results.txt
gasfit_ace_radiance_fitting_results.txt
gasfit_ace_spectrum_fit_residual_output.txt
spec_ace_zn_vis_2016_12_28.txt
sao2010_solref_air_380_480.dat
EOF

# Removing compilation by-products
rm -f *.o *.mod gasbeta.exe
