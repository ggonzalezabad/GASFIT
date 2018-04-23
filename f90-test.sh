gfortran -g -c -fbacktrace -fcheck=all -Wall gasfit_module.f90 main.f90
gfortran -g -fbacktrace -fcheck=all -Wall -o gasbeta.exe main.o gasfit_module.o


./gasbeta.exe << EOF
gasfit_ace_io.inp
gasfit_ace_io_output.txt
ace_zn_vis_2017_01_16_gasbeta.txt
test_spectrum_output.txt
sao2010_solref_air_380_480.dat
sao2010_solref_air_380_480.dat
newkpno_afgl.ring.bro
  1.d12, f, 1.d0
o3_221.vac.new.bro
  1.d21, f, 1.d0
o3_241.vac.new.bro
  1.d21, f, 1.d0
no2_220_iasb.vac.bro
  1.d18, f, 1.d0
bro228.10cm.bro
  1.d17, f, 1.d0
oclolt.dat.bro
  1.d18, f, 1.d0
O4_294K_BISA.BrO
  1.d5, f, 1.d0
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
