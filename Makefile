all:
	g++ /home/elliott/Research/opt/MOAB/build/tools/measure.o -g main.cpp sampling.cpp -fpic -I/home/elliott/Research/opt/MOAB/moab-4.6.2/tools/ -o samp -lMOAB

