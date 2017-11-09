
#------------------------------------------------
#  makefile for CITUP+SCITE
#------------------------------------------------

# change the following two lines according to your system
BOOST_DIRECTORY=/home/andrew/install/include/boost
CPLEX_DIRECTORY=/home/andrew/install/ibm/ILOG/CPLEX_Studio125
CPLEX_BUILD=x86-64_sles10_4.1

CC=/home/smalkic/install/gcc-5.2.0/bin/g++

scite.exe: *.cpp
	$(CC) -ggdb *.cpp -o scite.exe -O3 -DIL_STD -DILOSTRICTPOD -I../External/ -I$(CPLEX_DIRECTORY)/concert/include -I$(CPLEX_DIRECTORY)/cplex/include -L$(CPLEX_DIRECTORY)/concert/lib/$(CPLEX_BUILD)/static_pic/ -lconcert -L$(CPLEX_DIRECTORY)/cplex/lib/$(CPLEX_BUILD)/static_pic -l ilocplex -l cplex -l m -l pthread -std=c++11


