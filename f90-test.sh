gfortran -g -c -fbacktrace -fcheck=all -Wall gasfit_module.f90 main.f90
gfortran -g -fbacktrace -fcheck=all -Wall -o gasbeta.exe main.o gasfit_module.o

#gfortran -O3 -c gasfit_module.f90 main.f90
#gfortran -O3 -o gasbeta.exe main.o gasfit_module.o


./gasbeta.exe << EOF
gasfit_ace_io.inp
gasfit_ace_io_output.txt
gasfit_ace_geolocation_for_fitting_results.txt
gasfit_ace_solar_fitting_results.txt
gasfit_ace_radiance_fitting_results.txt
spec_ace_zn_vis_2016_12_19.txt
test_spectrum_output.txt
sao2010_solref_air_380_480.dat
res.bro228_10cm.boas
  1.d0, f, 1.d0
EOF

rm -f *.o *.mod gasbeta.exe
#diff 81003031.bro.part 81003031.bro.part.sav
#
# irrad: shift, 1 sigma =  -5.467635E-03  4.499953E-04
#  nvar_sun =            10
# irrad: hw1e, 1 sigma =   9.288516E-02  7.137651E-04
#   nfirst =            21
# rad: shift, 1 sigma =  -4.698179E-03  4.346268E-04
#  finished with iterate_sun
#                        Avg rms =   5.04303E-04
#                       Avg dgas =   2.27183E+13
#   Avg relative gas uncertainty =   1.43233E-01
#  Avg radiance shift uncertainty =  3.89668E-05
