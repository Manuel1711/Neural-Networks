CC = g++ -g -Wall #-I /home/users/tarante/armadillo-10.8.2/include/

all: boot.out

boot.out: bootstrap.cpp
	$(CC) bootstrap.cpp -oboot.out
	#./boot.out  -icorrelazioni_primofile.txt -otemp.dat
	#sort -k1 -k2 -n temp.dat > bootstrap_primofile.dat

ntk.out: main.cpp include/
	$(CC) main.cpp -ontk.out -larmadillo -llapack -lblas -fopenmp

dec.out: dec_Chebys.cpp
	$(CC) dec_Chebys.cpp -odec.out -larmadillo -llapack -lblas -fopenmp

rec.out: rec_Chebys.cpp
	$(CC) rec_Chebys.cpp -orec.out -larmadillo -llapack -lblas -fopenmp

cr_data.out: Data_Func.C
	$(CC) Data_Func.C -ocr_data.out 
