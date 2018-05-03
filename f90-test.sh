# Compilation for development:
#gfortran -g -c -fbacktrace -fcheck=all -Wall gasfit_module.f90 main.f90
#gfortran -g -fbacktrace -fcheck=all -Wall -o gasbeta.exe main.o gasfit_module.o
# Compilation for production
gfortran -O3 -c gasfit_module.f90 main.f90
gfortran -O3 -o gasbeta.exe main.o gasfit_module.o

# Running the program
./gasbeta.exe << EOF
gasfit.inp
Output_file1.txt
Output_file2.txt
Output_file3.txt
Output_file4.txt
Output_file5.txt
Input_file1.txt
sao2010_solref_air_380_480.dat
EOF

# Removing compilation by-products
rm -f *.o *.mod gasbeta.exe
