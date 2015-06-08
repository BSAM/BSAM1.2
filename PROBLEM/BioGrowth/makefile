#FOR = ifc -i_dynamic -O3 -tpp7 -IMODF -module MODF
#FOR = ifc -i_dynamic -CB -IMODF -module MODF
#FOR = /opt/intel/compiler70/ia32/bin/ifc -IMODF -module MODF
#FOR = ifort -fast -IMODF -module MODF -heap-arrays

#FOR = f95 -C=array -IMODF
#FOR = f95 -O3 -w -IMODF

#FOR = g95 -O3 -IMODF -fmod=MODF
#FOR = g95 -march=nocona -ffast-math -funroll-loops -O3 -IMODF -fmod=MODF
#FOR = g95 -fbounds-check -IMODF -fmod=MODF -ftrace=full

FOR = gfortran -IMODF -JMODF -O3 -ffast-math -funroll-loops
#FOR = gfortran -IMODF -JMODF -O3 -mtune=nocona -ffast-math -funroll-loops
#FOR = gfortran -IMODF -JMODF -O3 -mcpu=970 -mtune=970 -mpowerpc64 -mpowerpc-gpopt -force_cpusubtype_ALL -ffast-math -funroll-loops
#FOR = gfortran -IMODF -JMODF -fbounds-check

EXE = xbsam

LIB = ../../LIBRARY

OBJ = OF/nodeinfodef.o OF/problemdef.o OF/treeops.o  OF/gridutilities.o \
	OF/problem.o OF/bsamstorage.o OF/boundary.o OF/bsaminputoutput.o \
	OF/afasroutines.o OF/bsamroutines.o OF/bsamdriver.o

$(EXE): $(OBJ)
	$(FOR) $(OBJ) -o $(EXE)

OF/nodeinfodef.o: $(LIB)/nodeinfodef.f90
	$(FOR) -c $(LIB)/nodeinfodef.f90 -o OF/nodeinfodef.o

OF/problemdef.o: problemdef.f90
	$(FOR) -c problemdef.f90 -o OF/problemdef.o

OF/treeops.o: $(LIB)/treeops.f90
	$(FOR) -c $(LIB)/treeops.f90 -o OF/treeops.o

OF/gridutilities.o: $(LIB)/gridutilities.f90
	$(FOR) -c $(LIB)/gridutilities.f90 -o OF/gridutilities.o

OF/problem.o: problem.f90
	$(FOR) -c problem.f90 -o OF/problem.o

OF/bsamstorage.o: $(LIB)/bsamstorage.f90
	$(FOR) -c $(LIB)/bsamstorage.f90 -o OF/bsamstorage.o

OF/boundary.o: $(LIB)/boundary.f90
	$(FOR) -c $(LIB)/boundary.f90 -o OF/boundary.o

OF/bsaminputoutput.o: $(LIB)/bsaminputoutput.f90
	$(FOR) -c $(LIB)/bsaminputoutput.f90 -o OF/bsaminputoutput.o

OF/afasroutines.o: $(LIB)/afasroutines.f90
	$(FOR) -c $(LIB)/afasroutines.f90 -o OF/afasroutines.o

OF/bsamroutines.o: $(LIB)/bsamroutines.f90
	$(FOR) -c $(LIB)/bsamroutines.f90 -o OF/bsamroutines.o

OF/bsamdriver.o: $(LIB)/bsamdriver.f90
	$(FOR) -c $(LIB)/bsamdriver.f90 -o OF/bsamdriver.o

reset:
	rm $(EXE) MODF/* OF/*

remove:
	rm OUT/*.dat output.dat
