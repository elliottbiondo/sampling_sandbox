all:
	g++ /opt/moab/build/tools/measure.o -g main.cpp sampling.cpp -fpic -I/opt/moab/src/moab-4.6.2/tools/ -o samp -lMOAB
fortran:
	gfortran -c wrapper.F90
	g++ -c /opt/moab/build/tools/measure.o sampling.cpp -I/opt/moab/src/moab-4.6.2/tools/ -o sampling.o -lMOAB
	gfortran /opt/moab/build/tools/measure.o sampling.o wrapper.o -lstdc++ -lMOAB -o fort_samp
clean:
	rm -rf sampling.o wrapper.o fort_samp
