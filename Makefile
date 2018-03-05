#CC=mpic++ -std=c++0x -Wno-deprecated -g -L/soft/perftools/mpiP/lib -lmpiP -liberty
#CC=mpixlcxx -std=gnu++0x -Wall -Wno-deprecated -O3 
CC=mpic++ -std=gnu++0x -O3 

all	:swap 
swap 	: kmerGraph.o sequence.o mympi.o
	$(CC) kmerGraph.o sequence.o mympi.o  -o swap
	rm -f kmerGraph.o sequence.o mympi.o	

kmerGraph.o : src/kmerGraph.cpp src/kmerGraph.h src/mympi.h src/sequence.h src/bloomFilter.h
	$(CC) -c  src/kmerGraph.cpp -o kmerGraph.o
sequence.o : src/sequence.cpp src/sequence.h src/mympi.h
	$(CC) -c  src/sequence.cpp -o sequence.o
mympi.o : src/mympi.cpp src/mympi.h
	$(CC) -c  src/mympi.cpp -o mympi.o

clean :
	rm -f stats swap graph.o kmerGraph.o sequence.o mympi.o
