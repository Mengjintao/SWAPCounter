#ifndef _MYMPI_H_
#define _MYMPI_H_

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <math.h>
#include <sys/resource.h>
#include <errno.h>
#include <fstream>
#include <iostream>
#include <sstream>

#include <ctype.h>
#include <cstdlib>
#include <set>
#include <map>
#include <queue> 
#include <deque>
#include <string>
#include <vector>
#include <algorithm>
using namespace std;
 
typedef vector<int> vi;
typedef vector<vi> vvi;
typedef vector<vvi> vvvi;
typedef vector<bool> vb;
typedef vector<vb> vvb;
typedef vector<string> vs;
typedef vector<vs> vvs;

//vs tokS(string a, string d){vs r; char b[10000000],*s; strcpy(b,a.c_str()); for(s=strtok(b,d.c_str());s;s=strtok(0,d.c_str())) r.push_back(s); return r;}


class MPIEnviroment
{
public:
	//process variables
        int rank, nprocs;
	char processor_name[256];
	char hostname[255];
	int  memusage;
        int namelen, rc;
	
	//file variables
        MPI_File cFile;
        MPI_Offset size;
	long long int start_pos, end_pos, datasize;
	long long int read_offset;
	clock_t            elapsedTime;	
	clock_t            locateTime, readTime;
	//process status
        MPI_Status status;
	int groupProc;

public:
        void init(int argc, char **argv);
        void finalize();
        void File_open(char *File_name);
	void File_locate();
        void File_close();
        unsigned long long File_read(char *buf, long long n);
	void print(const char *message);
	void print(const char *message, FILE *fp);
	void File_write(char *File_name, const char *buf,unsigned long long n);
/*
	void MPI_Send_Lock(unsigned long long srcID, unsigned long long dstID);
	int  MPI_Recv_Lock_Reply(unsigned long long srcID, unsigned long long dstID);

	void MPI_Send_Leftdata(node *pnode);
	void MPI_Send_Rightdata(node *pnode);
	void MPI_Send_Unlock(unsigned long long NodeID);

	int  MPI_Recv_Type(unsigned long long *srcdst);
	int  MPI_Send_Lock_Success(unsigned long long src, unsigned long long dst);
	int  MPI_Send_Lock_Failed(unsigned long long src, unsigned long long dst);
	string MPI_recv_Edge(unsigned long long src, unsigned long long dst);
*/
};

#endif
