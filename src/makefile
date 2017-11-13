
#------------------------------------------------
#  makefile for B-SCITE
#------------------------------------------------

# set the values of following variables according to your system (example of their values on our system are shown below)
#CPLEX_DIRECTORY=/home/andrew/install/ibm/ILOG/CPLEX_Studio125
#CPLEX_BUILD=x86-64_sles10_4.1
#CC=/home/smalkic/install/gcc-5.2.0/bin/g++

CPLEX_DIRECTORY=
CPLEX_BUILD=
CC=

bscite.exe: *.cpp
	$(CC) -ggdb *.cpp -o bscite.exe -O3 -DIL_STD -DILOSTRICTPOD -I../External/ -I$(CPLEX_DIRECTORY)/concert/include -I$(CPLEX_DIRECTORY)/cplex/include -L$(CPLEX_DIRECTORY)/concert/lib/$(CPLEX_BUILD)/static_pic/ -lconcert -L$(CPLEX_DIRECTORY)/cplex/lib/$(CPLEX_BUILD)/static_pic -l ilocplex -l cplex -l m -l pthread -std=c++11
