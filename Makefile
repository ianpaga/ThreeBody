# This is an commentary line in a makefile

# Start of the makefile

Three_Body: Global.o Three_Body.o Func.o

	gfortran -o Three_Body *.o

Global.mod: Global.o Global.f90

	gfortran -c Global.f90

Global.o: Global.f90

	gfortran -c Global.f90

Newton.mod: Newton.o Newton.f90

	gfortran -c Newton.f90

Newton.o: Newton.f90

	gfortran -c Newton.f90

Func.mod: Func.o Func.f90

	gfortran -c Global.f90

Func.o: Global.mod Func.f90

	gfortran -c Func.f90

Init.mod: Init.o Init.f90

	gfortran -c Init.f90

Init.o: Global.o Func.o Newton.o Init.f90

	gfortran -c Init.f90

Three_Body.o: Global.mod Newton.mod Func.mod Init.mod Three_Body.f90

	gfortran -c Three_Body.f90

clean:

	rm *.mod *.o Three_Body

# End of the makefile
