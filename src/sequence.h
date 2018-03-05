#ifndef _SEQUENCE_H_
#define _SEQUENCE_H_

#include "mympi.h"
#include <sys/stat.h>
#include <sys/types.h>

#define DEBUG 0
#define BUF_SIZE         (1048576*32)
#define Filter_Threshold 0 
#define Edge_Threshold   6 
#define MersPrime        ((1<<31)-1)
#define READ_PROCESS_STEP 655360000
#define KMER_COMM_TYPE_LEN 10 
#define Contig_Length     100

class parameter
{
public:
	//file path
	char fastaPath[255];
	char outputPath[255];
	char contigsPath[255];
	char LogPath[255];
	char graphPath[255];
	char masterContigPath[255];
	char JungGraph_arc[255];
        char JungGraph_mul[255];
        char arcFrequency[255];
        char kmerGraph[255];
	char kmerCount[255];
	
  	//parameters
	int hashLength;
	unsigned long long MASK;
	int filterThreshold;
	int cutoffThreshold;
	int kmerGraphFlag;
	int JungGraphFlag;
	int distGraphFlag;
	int performanceFlag;
	int groupProc;
	
	//For Bloom Filter
	double refSize;
	int enableBF;
	int bitPerNode;

	//nucleotide variables
	char nucleotideValue[128];
	char nucleotideArray[4];
	char nucleotideReverse[128];

	void getopt(int argc, char **argv);
	void getParameters(int argc, char **argv, MPIEnviroment *MPIcontrol);
	
};
#endif
