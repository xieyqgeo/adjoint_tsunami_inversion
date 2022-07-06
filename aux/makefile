all:
	rm -f convertTsunamiDataTRInput conv_init_waterelev
	gfortran convertTsunamiDataTRInput.f90 -o convertTsunamiDataTRInput
	gfortran convertTimeReversalImageMatlabToInitalInput.f90 -o conv_init_waterelev
	gfortran convertdatares.f90 -o convert_data_res
