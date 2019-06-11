inexactproxgradient: inexactproxgradient.f90
	gfortran -c TVproximal.f90
	gfortran -O3 inexactproxgradient.f90 TVproximal.o -llapack
