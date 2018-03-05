#include "mympi.h"
#include "sequence.h"

void parameter::getopt(int argc, char **argv)
{
    this->fastaPath[0] = 0;
    this->outputPath[0] = 0;
    this->hashLength=23;
    this->cutoffThreshold = 5;
    this->groupProc = 1;
	
    this->kmerGraphFlag=0;
    this->JungGraphFlag=0;
    this->distGraphFlag=0;
    this->performanceFlag=0;
    this->refSize = 0;

    char comStr[100], localChar;
    int  ret;
    for(int i=1;i<argc;i++)
    {
	sscanf(argv[i], "%s", comStr);
 	if(comStr[0] != '-')  
	{
	  	printf("Usage: %s -k <kmerlength> -r <reference size> -c <cutoff Threshod> -g <grouped processes> -i <input Fasta data> -o <output data directory>\n",argv[0]);
//	    MPI_Abort(MPI_COMM_WORLD,1);		
	    exit(0);
	}	
	switch(comStr[1])
	{
	    case 'h':
	    	printf("Usage: %s -k <kmerlength> -r <reference size> -c <cutoff Threshod> -g <grouped processes> -i <input Fasta data> -o <output data directory>\n",argv[0]);
		exit(0);
	    case 'v':
		printf("SWAP-Assembler version 0.2\n");
		exit(0);
	    case 'k':
		if(i+1 == argc) 
		{
		    printf("Error: -k needs a value for the kmer size\n");
	    	    printf("Usage: %s -k <kmerlength> -r <reference size> -c <cutoff Threshod> -g <grouped processes> -i <input Fasta data> -o <output data directory>\n",argv[0]);
		    exit(0);
		}
		ret = sscanf(argv[i+1], "%d", &this->hashLength);
		i++;
		if(ret!=1 || this->hashLength<3 || this->hashLength>31 || this->hashLength%2!=1)	
		{
		    printf("Error: -k needs a value (odd number between 18 and 32) for the kmer size\n");
	    	    printf("Usage: %s -k <kmerlength> -r <reference size> -c <cutoff Threshod> -g <grouped processes> -i <input Fasta data> -o <output data directory>\n",argv[0]);
	            exit(0);
		}
		break;	
	    case 'c':
		if(i+1 == argc) 
		{
		    printf("Error: -c needs a value for the cutoff threshold\n");
	    	    printf("Usage: %s -k <kmerlength> -r <reference size> -c <cutoff Threshod> -g <grouped processes> -i <input Fasta data> -o <output data directory>\n",argv[0]);
		    exit(0);
		}
		ret = sscanf(argv[i+1], "%d", &this->cutoffThreshold);
		i++;
		if(ret!=1)	
		{
		    printf("Error: -c needs a value for the cutoff threshold (it is suggested to be 3~10 percent of the coverage of dataset)\n");
	    	    printf("Usage: %s -k <kmerlength> -r <reference size> -c <cutoff Threshod> -g <grouped processes> -i <input Fasta data> -o <output data directory>\n",argv[0]);
	            exit(0);
		}
		break;
            case 'i':
		if(i+1 == argc) 
		{
		    printf("Error: -i needs a path for the file in fasta format\n");
	    	    printf("Usage: %s -k <kmerlength> -r <reference size> -c <cutoff Threshod> -g <grouped processes> -i <input Fasta data> -o <output data directory>\n",argv[0]);
		    exit(0);
		}
		ret = sscanf(argv[i+1], "%s", this->fastaPath);
		i++;
		if(ret==0)
		{
		    printf("Error: -i needs a path for the file in fasta format\n");
	    	    printf("Usage: %s -k <kmerlength> -r <reference size> -c <cutoff Threshod> -g <grouped processes> -i <input Fasta data> -o <output data directory>\n",argv[0]);
		    exit(0);
		}
		break;
            case 'o':
		if(i+1 == argc) 
		{
		    printf("Error: -o needs a directory name used for the output data\n");
	    	    printf("Usage: %s -k <kmerlength> -r <reference size> -c <cutoff Threshod> -g <grouped processes> -i <input Fasta data> -o <output data directory>\n",argv[0]);
		    exit(0);
		}
		ret = sscanf(argv[i+1], "%s", this->outputPath);
		i++;
		if(ret==0)
		{
		    printf("Error: -o needs a directory path for the output data\n");
	    	    printf("Usage: %s -k <kmerlength> -r <reference size> -c <cutoff Threshod> -g <grouped processes> -i <input Fasta data> -o <output data directory>\n",argv[0]);
		    exit(0);
		}
		break;
	    case 'r':
		if(i+1 == argc) 
		{
		    printf("Error: -o needs a directory name used for the output data\n");
	    	    printf("Usage: %s -k <kmerlength> -r <reference size> -c <cutoff Threshod> -g <grouped processes> -i <input Fasta data> -o <output data directory>\n",argv[0]);
		    exit(0);
		}
		localChar = 0;
		ret = sscanf(argv[i+1], "%lf%c", &this->refSize, &localChar);
		if(localChar=='G')	this->refSize *= (1024*1024*1024);
		if(localChar=='M')	this->refSize *= (1024*1024);
		if(localChar=='K')	this->refSize *= 1024;
		i++;
		if(ret==0)
		{
		    printf("Error: -r needs a reference size for the Bloom filter\n");
	    	    printf("Usage: %s -k <kmerlength> -r <reference size> -c <cutoff Threshod> -g <grouped processes> -i <input Fasta data> -o <output data directory>\n",argv[0]);
	            exit(0);
		}
		break;
	    case 'g':
		if(i+1 == argc) 
		{
		    printf("Error: -g needs a value for the grouped processes number for reading input data (suggested value should be larger than 1 and less than the number of cores in computing nodes)\n");
	    	    printf("Usage: %s -k <kmerlength> -r <reference size> -c <cutoff Threshod> -g <grouped processes> -i <input Fasta data> -o <output data directory>\n",argv[0]);
		    exit(0);
		}
		ret = sscanf(argv[i+1], "%d", &this->groupProc);
		i++;
		if(ret!=1)	
		{
		    printf("Error: -g needs a value for the grouped processes number for reading input data (suggested value should be larger than 1 and less than the number of cores in computing nodes)\n");
	    	    printf("Usage: %s -k <kmerlength> -r <reference size> -c <cutoff Threshod> -g <grouped processes> -i <input Fasta data> -o <output data directory>\n",argv[0]);
	            exit(0);
		}
		break;
		this->groupProc = 4;
	    case 's':
		this->kmerGraphFlag = 1;
		break;
	    case 'j':
		this->JungGraphFlag = 1;
		break;
	    case 'd':
		this->distGraphFlag = 1;
		break;
	    case 'p':
		this->performanceFlag = 1;
		break;
		
	}	
    }
    if(this->fastaPath[0]==0 || this->outputPath[0]==0)
    {
        printf("Usage: %s -k <kmerlength> -r <reference size> -c <cutoff Threshod> -i <input Fasta data> -o <output data directory>\n",argv[0]);
        exit(0);
    }

}

void parameter::getParameters(int argc, char **argv, MPIEnviroment *MPIcontrol)
{
    
    getopt(argc, argv);

    if(MPIcontrol->rank==0)
    {    
    	int stat = mkdir(this->outputPath, S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);
    	if(stat==-1)	
    	{
		printf("Error: Creating directory %s\n", this->outputPath);
		exit(0);
    	} 
    }

    strcpy(this->contigsPath, this->outputPath);
    strcpy(this->LogPath, this->outputPath);
    strcpy(this->graphPath, this->outputPath);
    strcpy(this->masterContigPath,this->outputPath);
    strcpy(this->JungGraph_arc, this->outputPath);
    strcpy(this->JungGraph_mul, this->outputPath);
    strcpy(this->arcFrequency, this->outputPath);
    strcpy(this->kmerGraph, this->outputPath);
    strcpy(this->kmerCount, this->outputPath);
    
    strcat(this->contigsPath, "/noCEcontigs.fasta");
    strcat(this->LogPath, "/logtime.txt"); 
    strcat(this->graphPath, "/contigGraph.txt");
    strcat(this->masterContigPath,"/CEContig.fasta");
    strcat(this->JungGraph_arc, "/JungGraph_arc.txt");
    strcat(this->JungGraph_mul, "/JungGraph_mul.txt");
    strcat(this->arcFrequency, "/arcFrequency.txt");    
    strcat(this->kmerGraph, "/kmerGraph.txt");
    strcat(this->kmerCount, "/kmerCount.txt");

    MPIcontrol->groupProc = this->groupProc;

    if(MPIcontrol->rank==0)
    {    printf("Runing command: %s -k %d -c %d -i %s -o %s\n", argv[0], this->hashLength, this->cutoffThreshold, this->fastaPath, this->outputPath);
         printf("contigsPath=%s, logPath=%s, graphPath=%s\n", this->contigsPath, this->LogPath, this->graphPath);
    }
//    strcpy(this->LogPath, "logtime.txt");
//    strcpy(this->graphPath, "distGraph.txt");

    char ch0='A', ch1='C', ch2='G', ch3='T';
    this->nucleotideValue[(int)ch0] = 0;
    this->nucleotideValue[(int)ch1] = 1;
    this->nucleotideValue[(int)ch2] = 2;
    this->nucleotideValue[(int)ch3] = 3;

    this->nucleotideArray[0] = 'A';
    this->nucleotideArray[1] = 'C';
    this->nucleotideArray[2] = 'G';
    this->nucleotideArray[3] = 'T';

    this->nucleotideReverse[(int)ch0] = 'T';
    this->nucleotideReverse[(int)ch1] = 'G';
    this->nucleotideReverse[(int)ch2] = 'C';
    this->nucleotideReverse[(int)ch3] = 'A';

    this->MASK = ~(3ull<<(2*hashLength));
    this->filterThreshold = Filter_Threshold;
 //   this->arcThreshold    = 5;
}
