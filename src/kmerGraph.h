#ifndef _KMERGRAPH_H_
#define _KMERGRAPH_H_

#include <map>
#include <string>
#include <ext/hash_map>
#include <unistd.h>
#include <unordered_map>
#include "sequence.h"
#include "mympi.h"
#include "bloomFilter.h"

using namespace std;

class kmerGraph
{
public:
	unordered_map<unsigned long long, unsigned long long> kmolecules;	
	unordered_map<unsigned long long, unsigned long long> tmp;	
	long long size;
	unsigned long long read_pos;
	double commtime, commtimetot;
	double cuttime, storagetime;
	MPI_Datatype commType;

	unsigned long long reverseComplement(unsigned long long kmerDescriptor, parameter *parameters);
    	static unsigned long long stringToLongLong(const char *buf, int start, int end, parameter *parameters);
	static string longLongToString(unsigned long long a, parameter *parameters);
	unsigned long long getProcsID(unsigned long long kmerID, int hashLength, MPIEnviroment *MPIcontrol);
	int arcPos(unsigned long long &A, unsigned long long &B, char directA, char directB, int hashLength);
	long long addKmers(const char *s, const char *end, unsigned long long *kmers, unsigned char *arcs, parameter *parameters);

	kmerGraph()
	{
		size     = 0;
		read_pos = 0;
		commtime = 0;
		commtimetot = 0;
		cuttime = 0;
		storagetime = 0;	
		
		MPI_Type_contiguous(KMER_COMM_TYPE_LEN, MPI_UNSIGNED_LONG_LONG, &commType);
		MPI_Type_commit(&commType);
	}
	
	~kmerGraph()
	{
	}

public:
	int  constructKmerGraph(parameter *parameters, MPIEnviroment *MPIcontrol, bloomFilter *bFilter);
	void printKmerCounting(parameter *parameters, MPIEnviroment *MPIcontrol);
	void distributeKmerGraph(parameter *parameters, MPIEnviroment *MPIcontrol);
	void printKmerFreq(parameter *parameters, MPIEnviroment *MPIcontrol);
	void extractKmerGraph(parameter *parameters, MPIEnviroment *MPIcontrol, bloomFilter *bFilter, char *readBuf, long long int localEnd, long long int readRound, clock_t & commTime, clock_t & checkTime);
	void dataShuttleBetweenGroupProcs(parameter *parameters, MPIEnviroment *MPIcontrol, char *& readBuf, long long int & localEnd, long long int readRound);
};

#endif
