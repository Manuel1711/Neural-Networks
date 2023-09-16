CC = g++ -g -Wall #-I /home/users/tarante/armadillo-10.8.2/include/

all: boot.out

boot.out: bootstrap.cpp
	$(CC) bootstrap.cpp -oboot.out
	#./boot.out  -icorrelazioni_primofile.txt -otemp.dat
	#sort -k1 -k2 -n temp.dat > bootstrap_primofile.dat
