CC=mpic++ 
CFLAGS=-std=c++11 -lmpi

a2: final.cpp
	$(CC) $(CFLAGS) final.cpp -o a2

run: final.cpp
	mpirun -np 8 ./a2 --taskid=1 --inputpath=test2_1/test-input-2.gra --headerpath=test2_1/test-header-2.dat --outputpath=output-2.txt --verbose=1 --startk=1 --endk=8 --p=4

clean:
	rm -rf *.o a2