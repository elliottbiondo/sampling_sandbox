all:
	g++ /opt/moab/build/tools/measure.o -g main.cpp sampling.cpp -fpic -I/opt/moab/src/moab-4.6.2/tools/ -o samp -lMOAB

