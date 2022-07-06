all:
	rm -f pcomcot
	mpif90 -I/opt/local/include -L/opt/local/lib -lnetcdff VariableDefination.f90 okada.f pcomcotNetCDFlib.f90 pcomcotLIB.f90 solvEQNs.f90 pcomcot.f90  -o pcomcot
	rm -f VariableDefination.mod variabledefination.mod
nocdf:
	rm -f pcomcot
	mpif90 -FC=ifort VariableDefination.f90 okada.f pcomcotNetCDFlibEmpty.f90 pcomcotLIB.f90 solvEQNs.f90 pcomcot.f90 -O3 -o pcomcot
	rm -f VariableDefination.mod variabledefination.mod
checkbounds:
	rm -f pcomcot
	mpif90 -fcheck=bounds -I/opt/local/include -L/opt/local/lib -lnetcdff VariableDefination.f90 okada.f pcomcotNetCDFlib.f90 pcomcotLIB.f90 solvEQNs.f90 pcomcot.f90  -o pcomcot
	rm -f VariableDefination.mod variabledefination.mod
