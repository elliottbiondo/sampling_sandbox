all:
	g++ /opt/moab/build/tools/measure.o -g sampling.cpp -I/opt/moab/src/moab-4.6.2/tools/ -o samp -lMOAB
