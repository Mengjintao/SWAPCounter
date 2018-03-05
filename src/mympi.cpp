#include "mpi.h"
#include "mympi.h"

#define ReadLength     100
#define PreReadLen    (ReadLength*3)

void MPIEnviroment::init(int argc, char *argv[])
{
	elapsedTime = 0;
	locateTime = 0;
	readTime = 0;
	groupProc=1;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Get_processor_name(processor_name, &namelen);
}

void MPIEnviroment::print(const char *message)
{
	int who = RUSAGE_SELF;
	struct rusage usage;
	int usageErrno;
	usageErrno = getrusage(who, &usage);
	if(usageErrno == EFAULT)  		printf("Error: EFAULT\n");
	else if(usageErrno == EINVAL) 		printf("Error: EINVAL\n");

	printf("Proc:%d (%s)[Mem=%ld KB] -> %s\n", rank, processor_name, usage.ru_maxrss, message); 
}

void MPIEnviroment::print(const char *message, FILE *fp)
{
	int who = RUSAGE_SELF;
	struct rusage usage;
	int usageErrno;
	usageErrno = getrusage(who, &usage);
	if(usageErrno == EFAULT)  		printf("Error: EFAULT\n");
	else if(usageErrno == EINVAL) 		printf("Error: EINVAL\n");

	fprintf(fp, "Proc:%d [Mem=%ld KB] ->%s\n", rank, usage.ru_maxrss, message);	
}

void MPIEnviroment::finalize()
{
	assert(MPI_Finalize()==MPI_SUCCESS);
}

void MPIEnviroment::File_open(char *File_name)
{
	if(rank%this->groupProc)	return;
	clock_t startTime, endTime;

	startTime = clock();	
        rc = MPI_File_open(MPI_COMM_SELF, File_name, MPI_MODE_RDONLY, MPI_INFO_NULL, &cFile);

        if(rc) {
                printf("Proc:%d Unable to open file \"%s\"\n", rank, File_name); fflush(stdout);
                exit(0);
        }
        else    {
                printf("Proc:%d file\"%s\" opened ...\n", rank, File_name);
        }
	endTime = clock();
	elapsedTime += endTime-startTime;
}

void MPIEnviroment::File_locate()
{
	clock_t startTime, endTime;
	startTime = clock();		

        size = -1;
        if(rank==0)	MPI_File_get_size(cFile, &size);
	MPI_Bcast(&size,1,MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);	
	printf("Proc:%d size of file is: %lld ...\n", rank,  size);
	fflush(stdout);

	if(rank%this->groupProc)	{
		start_pos = end_pos = 0;
		datasize = 0;
//        	printf("Proc:%d file start at %llu and end at %llu...\n", rank,  start_pos, end_pos);
		return;
	}
        unsigned long long step = (size-1)/(nprocs/groupProc)+1;
	int groupRank = rank/groupProc;
	//set start position and end position for file pointer
	start_pos = groupRank*step;
	start_pos = start_pos<size?start_pos:size;
	end_pos   = (groupRank+1)*step;
	end_pos   = end_pos<size?end_pos:size;

	//Find the start position for one read.
	char str[PreReadLen+2];
	unsigned long long  preLen = size-start_pos;
	if(preLen>PreReadLen)	preLen=PreReadLen;	
	MPI_File_read_at(cFile, start_pos, str, preLen, MPI_CHAR, &status);
	unsigned long long i=0;	
	while(str[i]!='>' && start_pos<size && i<preLen) {i++; start_pos++;}
	assert(str[i]=='>'||start_pos==size||i==preLen);
	
	//Find the end position for one read.
	preLen = min(size-end_pos, (long long int) PreReadLen);
	MPI_File_read_at(cFile,end_pos,str,preLen, MPI_CHAR, &status);
	i=0;
        while(str[i]!='>' && end_pos<size && i<preLen) {i++; end_pos++;}

	datasize = end_pos-start_pos;
//        printf("Proc:%d file start at %llu and end at %llu...\n", rank,  start_pos, end_pos);
	read_offset = start_pos;

	endTime = clock();
	elapsedTime += endTime-startTime;
	locateTime  += endTime-startTime;	
}

unsigned long long MPIEnviroment::File_read(char *buf, long long int n)
{
	if(rank%this->groupProc)	return 0;
	clock_t startTime, endTime;
	startTime = clock();
	long long int len = n<datasize? n : datasize;
	MPI_File_read_at_all(cFile, read_offset, buf, len, MPI_CHAR, &status);
	read_offset += len;	
	datasize -= len;
	endTime = clock();
	elapsedTime += endTime-startTime;
	readTime += endTime-startTime;
	return len;
}

void MPIEnviroment::File_write(char *filename, const char *buf, unsigned long long int n)
{
	clock_t startTime, endTime;
	startTime = clock();
	MPI_File writeFile;
	rc = MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &writeFile);
//        rc = MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE|MPI_MODE_APPEND, MPI_INFO_NULL, &writeFile);
        if(rc) {
                printf("Proc:%d Unable to open file \"%s\"\n", rank, filename); fflush(stdout);
                exit(0);
        }
        else    {
                printf("Proc:%d file\"%s\" was opened ...\n", rank, filename);
        }

	char *localBuf = const_cast<char*> (buf);
	unsigned long long  startPos;
	MPI_Scan(&n, &startPos, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
	startPos = startPos-n;

	MPI_Status writeStatus;
	MPI_File_write_at_all(writeFile, startPos, localBuf, n, MPI_CHAR, &writeStatus);
	MPI_File_close(&writeFile);

	endTime = clock();
	elapsedTime += endTime-startTime;
}

void MPIEnviroment::File_close()
{
	MPI_File_close(&cFile);
}

