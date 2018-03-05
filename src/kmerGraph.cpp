#include "sequence.h"
#include "mympi.h"
#include "kmerGraph.h"
#include <math.h>

unsigned long long kmerGraph::reverseComplement(unsigned long long a, parameter *parameters)
{
        unsigned long long rev = 0;
        for(int i=0; i<parameters->hashLength; i++){
                rev <<= 2;
                rev |= (3-(a&3));
                a >>= 2;
        }
        return rev;
}

unsigned long long kmerGraph::stringToLongLong(const char *buf, int start, int end, parameter *parameters)
{
        unsigned long long ret = 0;
        for(int i=start; i<end; i++){
                assert(buf[i]=='A'||buf[i]=='T'||buf[i]=='C'||buf[i]=='G');
		ret <<= 2;
                ret |= ((unsigned long long)parameters->nucleotideValue[(int)buf[i]]&3);
        }
        return ret;
}

string kmerGraph::longLongToString(unsigned long long a, parameter *parameters)
{
	string descriptor;
        descriptor.clear();
        for(int i=0;i<parameters->hashLength;i++)
        {
                descriptor += parameters->nucleotideArray[a%4];
                a = (a>>2);
        }
        reverse(descriptor.begin(), descriptor.end());
	return descriptor;
}

int kmerGraph::arcPos(unsigned long long &A, unsigned long long &B, char directA, char directB, int hashLength)
{
	int lastB;
	if(directB=='+')	lastB = B%4;
	else			lastB = 3-(B>>((hashLength-1)*2));
	
	if(directA=='-')	lastB += 4;
	assert(lastB>=0&&lastB<8);
	return lastB;
}

/*
unsigned long long kmerGraph::getProcsID(unsigned long long kmerID, int hashLength, MPIEnviroment *MPIcontrol)
{
	unsigned long long ret = 0, tmpID = kmerID;
	for(int i=0;i<hashLength;i++)
	{
		ret = ret * 4;
//		ret |= ((unsigned long long)3 - (tmpID&(unsigned long long)3));
		ret += (3 - tmpID%4);
		tmpID = tmpID / 4;	
	}
//	ret = ret ^ kmerID;
//	return ( (ret % (unsigned long long) MersPrime) % (unsigned long long) MPIcontrol->nprocs);	
//	return ( (ret ) % (unsigned long long) MPIcontrol->nprocs);	
	if(ret<kmerID)	ret = kmerID;
	
	double tmp = ((sqrt(5.0)-1)/2) * ret ;
	double rr  = tmp - floor(tmp);
	return floor(rr * MPIcontrol->nprocs);
}
*/

unsigned long long kmerGraph::getProcsID(unsigned long long kmerID, int hashLength, MPIEnviroment *MPIcontrol)
{
        unsigned long long revKmer = 0, tmpID = kmerID;
        for(int i=0;i<hashLength;i++)
        {
                revKmer = revKmer * 4;
                revKmer += ( 3 - tmpID%4 );
                tmpID = tmpID / 4;
        }

        if(revKmer>kmerID) revKmer = kmerID;

         unsigned int factor = 19;
         unsigned int numBytes = (hashLength + 3) / 4;

         unsigned long long sum = 0;
         for(unsigned i = 0; i < numBytes; i++){
                 sum = sum * factor + (revKmer & 0xFF);
                 revKmer >>= 8;
         }
         return sum % MPIcontrol->nprocs;
}

/*
long long kmerGraph::addKmers(const char *s, const char *end, unsigned long long *kmers, unsigned char *arcs, parameter *parameters)
{
	unsigned long long curNodeID=0,     primeNodeID=0,   	twinNodeID=0;
	unsigned long long preNodeID=0,     prePrimeNodeID=0; 
	unsigned long long kmer_index=0;
	
	curNodeID = this->stringToLongLong(s, 0, parameters->hashLength-1, parameters);
	int len = end - s;
        for(int j=parameters->hashLength-1;j<len; j++) 
	{
		if(!(s[j]=='A'||s[j]=='T'||s[j]=='C'||s[j]=='G'))	continue;	
	       	curNodeID <<= 2;
               	curNodeID |= (parameters->nucleotideValue[(int) s[j]]&(3ull));
               	curNodeID &= parameters->MASK;
		
               	twinNodeID  = this->reverseComplement(curNodeID, parameters);
               	primeNodeID = max(curNodeID, twinNodeID);
		kmers[kmer_index]  = primeNodeID;
                kmer_index++;
	}
	return (long long int) kmer_index;
}
*/

long long kmerGraph::addKmers(const char *s, const char *end, unsigned long long *kmers, unsigned char *arcs, parameter *parameters)
{
//	if(s>end || end-s <= parameters->hashLength)	return 0;
	unsigned long long curNodeID=0,     primeNodeID=0,   	twinNodeID=0;
	unsigned long long preNodeID=0,     prePrimeNodeID=0; 
	unsigned long long kmer_index=0;
//	for(int i=0; i<end-s; i++)	printf("%c", *(s+i));
//	printf("\n");
	
	curNodeID = this->stringToLongLong(s, 0, parameters->hashLength-1, parameters);
        for(int j=parameters->hashLength-1; s+j<end; j++) 
	{
		if(!(s[j]=='A'||s[j]=='T'||s[j]=='C'||s[j]=='G'))	continue;	
		assert(s[j]=='A'||s[j]=='T'||s[j]=='C'||s[j]=='G');
	       	curNodeID <<= 2;
               	curNodeID |= (parameters->nucleotideValue[(int) s[j]]&(3ull));
               	curNodeID &= parameters->MASK;
		
               	twinNodeID  = this->reverseComplement(curNodeID, parameters);
               	primeNodeID = max(curNodeID, twinNodeID);
//		cout<<longLongToString(curNodeID, parameters)<<" "<<longLongToString(primeNodeID, parameters)<<endl;
	
		if(j-parameters->hashLength+1==0)	{
			preNodeID = curNodeID;
			prePrimeNodeID = primeNodeID; 
			kmers[kmer_index] = primeNodeID;	arcs[kmer_index] = 0;
			kmer_index ++;
			continue;
		}
		kmers[kmer_index]  = primeNodeID;	arcs[kmer_index]  = 0;

		int pos_ret;
		//preNode -> curNode  add(A, B, +, +, last(B+))
		if( preNodeID == prePrimeNodeID && curNodeID == primeNodeID)
		{
			pos_ret = arcPos(prePrimeNodeID, primeNodeID, '+', '+', parameters->hashLength);	
			assert(kmer_index-1>=0);
			arcs[kmer_index-1] |= ((unsigned char) 1<<pos_ret);
				
			pos_ret = arcPos(primeNodeID, prePrimeNodeID, '-', '-', parameters->hashLength); 
			arcs[kmer_index] |=  ((unsigned char) 1<<pos_ret);
		}
		//preNode -> ~curNode
		else if(preNodeID == prePrimeNodeID && curNodeID != primeNodeID)
		{
			pos_ret = arcPos(prePrimeNodeID, primeNodeID, '+', '-', parameters->hashLength);	
			assert(kmer_index-1>=0);
			arcs[kmer_index-1] |= ((unsigned char) 1<<pos_ret);
				
			pos_ret = arcPos(primeNodeID, prePrimeNodeID, '+', '-', parameters->hashLength); 
			arcs[kmer_index] |= ((unsigned char) 1<<pos_ret);
		}
		//~preNode -> curNode
		else if(preNodeID != prePrimeNodeID && curNodeID == primeNodeID)
		{
			pos_ret = arcPos(prePrimeNodeID, primeNodeID, '-', '+', parameters->hashLength);	
			assert(kmer_index-1>=0);
			arcs[kmer_index-1] |= ((unsigned char) 1<<pos_ret);
				
			pos_ret = arcPos(primeNodeID, prePrimeNodeID, '-', '+', parameters->hashLength); 
			arcs[kmer_index] |= ((unsigned char) 1<<pos_ret);
		}			
		//~preNode -> ~curNode
		else if(preNodeID != prePrimeNodeID && curNodeID != primeNodeID)
		{
			pos_ret = arcPos(prePrimeNodeID, primeNodeID, '-', '-', parameters->hashLength);	
			assert(kmer_index-1>=0);
			arcs[kmer_index-1] |= ((unsigned char) 1<<pos_ret);
				
			pos_ret = arcPos(primeNodeID, prePrimeNodeID, '+', '+', parameters->hashLength); 
			arcs[kmer_index] |=  ((unsigned char) 1<<pos_ret);
		}
		preNodeID = curNodeID;
		prePrimeNodeID = primeNodeID;
		kmer_index++;
	}
	return (long long int) kmer_index;
}


void kmerGraph::extractKmerGraph(parameter *parameters, MPIEnviroment *MPIcontrol, bloomFilter *bFilter, char *readBuf, long long int localEnd, long long int readRound, clock_t & commTime, clock_t & checkTime)
{
		clock_t t0, t1;
		clock_t tt0, tt1;

		long long int tot_kmer_num = 0, localTag=0, localLen=0;
		for(long long int i=0; i<localEnd&&readBuf; i++)
		{
			if(readBuf[i]=='>')	
			{	
				localTag=0;	
				if(localLen <= parameters->hashLength+1)	continue;
				tot_kmer_num += localLen-parameters->hashLength;
				localLen = 0;
			}
			else if(localTag==0&&readBuf[i]=='\n')	localTag=1, localLen=0;
			else if(localTag==1)			localLen++;
		}	
		if(localLen>parameters->hashLength+1)	tot_kmer_num += localLen-parameters->hashLength;
		size  = tot_kmer_num;

		unsigned long long *kmers = new unsigned long long [size];
		unsigned char      *arcs  = new unsigned char      [size];
		assert(kmers && arcs);	
	
		tt0 = clock();		
		tot_kmer_num=0, localTag=0, localLen=0;
		for(long long int i=0; i<localEnd&&readBuf; i++)
		{
			if(readBuf[i]=='>')	
			{	
				localTag=0;	
				if(localLen <= parameters->hashLength+1)	continue;
				long long r = addKmers(readBuf+i-localLen, readBuf+i-1, kmers+tot_kmer_num, arcs+tot_kmer_num, parameters); 
			
				if(r!=localLen-parameters->hashLength)
				{
					printf("------r = %ld %ld %ld %ld-----------\n", r, localLen, localLen-parameters->hashLength, i);
					for(char *p=readBuf+i-localLen; p<=readBuf+i;p++)	printf("%c", *p);
					printf("\n");
					fflush(stdout);
				}
				assert(r==localLen-parameters->hashLength);
				localLen = 0;
				tot_kmer_num += r;
			}
			else if(localTag==0&&readBuf[i]=='\n')	localTag=1, localLen=0;
			else if(localTag==1)			
			{
				if(readBuf[i]=='A'||readBuf[i]=='T'||readBuf[i]=='C'||readBuf[i]=='G'|| readBuf[i]=='\n');
				else	readBuf[i]='A';				
				localLen++;
			}
		}	
	
		if(localLen>parameters->hashLength+1&&readBuf)	{
//			if(readBuf[localEnd]!='>')	localLen--;
			long long r = addKmers(readBuf+localEnd-localLen, readBuf+localEnd-1, kmers+tot_kmer_num, arcs+tot_kmer_num, parameters); 
			tot_kmer_num += r;
		//	assert(r==localLen-parameters->hashLength);
		}
		assert(size==tot_kmer_num);
		size  = tot_kmer_num;
//		readRound++;

		tt1 = clock();
//		MPI_Barrier(MPI_COMM_WORLD);
//		if(MPIcontrol->rank==0)	printf("Proc %d: OK in Barrier 2, size=%lld\n", MPIcontrol->rank, size);

		//arrange kmers and arcs before the communication step
		unsigned long long *kmers_send = new unsigned long long [size];
		assert(kmers_send != NULL);
		unsigned char      *arcs_send  = new unsigned char      [size];
		assert(arcs_send != NULL);
		for(long long i=0; i < size; i++)	kmers_send[i] = arcs_send[i] = 0;

		int  *send_size = new int [MPIcontrol->nprocs];
		assert(send_size!=NULL);
		int  *recv_size = new int [MPIcontrol->nprocs];
		assert(recv_size!=NULL);
		int  *send_pos  = new int [MPIcontrol->nprocs];
		assert(send_pos!=NULL);
		int  *recv_pos  = new int [MPIcontrol->nprocs];
		assert(recv_pos!=NULL);

		for(int i=0; i<MPIcontrol->nprocs; i++)	send_size[i] = recv_size[i] = 0;

//		if(MPIcontrol->rank==0)	printf("Proc %d: OK in Barrier 2 size=%llu\n", MPIcontrol->rank, size);
//		fflush(stdout);
//		MPI_Barrier(MPI_COMM_WORLD);

		unsigned long long send_sum=0,  recv_sum=0;
		unsigned long long ProcsID;	

		for(long long i=0; i<size; i++)
		{
			ProcsID = getProcsID(kmers[i], parameters->hashLength, MPIcontrol);
			assert(ProcsID>=0 && ProcsID<(unsigned long long) MPIcontrol->nprocs);
			send_size[ProcsID]++;
		}

		checkTime += tt1-tt0;
		t0 = clock();
//		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Alltoall(send_size, 1, MPI_INT, recv_size, 1, MPI_INT, MPI_COMM_WORLD);
		t1 = clock();
		commTime += (t1-t0);
		tt0 = clock();
		for(int i=0; i<MPIcontrol->nprocs; i++)	send_sum += send_size[i];
		for(int i=0; i<MPIcontrol->nprocs; i++)	recv_sum += recv_size[i];
		assert(send_sum == (unsigned long long) size);
	
//		MPI_Barrier(MPI_COMM_WORLD);
//		if(MPIcontrol->rank==0)	printf("Proc %d: OK in Barrier 3\n", MPIcontrol->rank);
//		printf("Proc %d: send_sum=%llu recv_sum=%llu\n", MPIcontrol->rank, send_sum, recv_sum);

		send_pos[0] = 0;
		for(int i=1; i<MPIcontrol->nprocs; i++)	send_pos[i] = send_pos[i-1] + send_size[i-1];
		recv_pos[0] = 0;
		for(int i=1; i<MPIcontrol->nprocs; i++)	recv_pos[i] = recv_pos[i-1] + recv_size[i-1];
//	        printf("proc%d: ", MPIcontrol->rank);
//        	for(int i=0;i<size;i++) printf("send_size[%d]=%d ", i, send_size[i]);
//        	printf("\t");
//		printf("proc%d: ", MPIcontrol->rank);
//       	for(int i=0;i<size;i++) printf("recv_size[%d]=%d ", i, recv_size[i]);

		for(long long i=0; i<size; i++)
		{
			ProcsID = getProcsID(kmers[i], parameters->hashLength, MPIcontrol);
			kmers_send[send_pos[ProcsID]]  = kmers[i]; 
			arcs_send[send_pos[ProcsID]++] = arcs[i];	
		}

		send_pos[0] = 0;
		for(int i=1; i<MPIcontrol->nprocs; i++)	send_pos[i] = send_pos[i-1] + send_size[i-1];

		delete [] kmers;
		delete [] arcs;
		kmers = NULL;  arcs=NULL;	
		if(recv_sum==0)	printf("Proc %d: recv_sum = %lld \n", MPIcontrol->rank, recv_sum);
		fflush(stdout);

//		MPI_Barrier(MPI_COMM_WORLD);
		kmers = new unsigned long long [recv_sum+1];
		assert(kmers!=NULL);
		arcs  = new unsigned char      [recv_sum+1];
		assert(arcs!=NULL);
		
//		MPI_Barrier(MPI_COMM_WORLD);
//		if(MPIcontrol->rank==0)	printf("Proc %d: OK in Barrier 4\n", MPIcontrol->rank);

		tt1 = clock();
//		checkTime += tt1-tt0;
		t0 = clock();
		MPI_Alltoallv(kmers_send,send_size,send_pos, MPI_LONG_LONG_INT, kmers,recv_size,recv_pos, MPI_LONG_LONG_INT, MPI_COMM_WORLD);	
		MPI_Alltoallv(arcs_send, send_size,send_pos, MPI_CHAR,          arcs, recv_size,recv_pos, MPI_CHAR,          MPI_COMM_WORLD);
	//	printf("Proc %d: AlltoAll Finished\n", MPIcontrol->rank);
		t1 = clock();
		commTime += (t1-t0);

		delete [] kmers_send;
		delete [] arcs_send;
		kmers_send = NULL;	arcs_send = NULL;
		size  = recv_sum;
	
		delete [] send_size;
		delete [] recv_size;
		delete [] send_pos;
		delete [] recv_pos;
		send_size = recv_size = send_pos = recv_pos = NULL;

		tt0 = clock();
		for(long long i=0;i<size;i++)
		{
			tmp[kmers[i]]++;
//			if(MPIcontrol->rank==39)	{printf("Proc0_kmer %lu\n", kmers[i]); fflush(stdout);}

			if(kmolecules.find(kmers[i])!=kmolecules.end())	 kmolecules[kmers[i]]++;
			else 
			{
			if(bFilter->bloomAdd(kmers[i]) < parameters->cutoffThreshold)	continue;	
//			int ProcsID = getProcsID(kmers[i], parameters->hashLength, MPIcontrol);
//			assert(ProcsID == MPIcontrol->rank);
//			unsigned char *p = (unsigned char *) &kmolecules[kmers[i]];
//			for(int j=0; j<8; j++) 	if(p[j]<250)	p[j] += ((arcs[i]&(1<<j))?1:0);
			else			kmolecules[kmers[i]]=parameters->cutoffThreshold;
			}
		}
	
		tt1 = clock();
		checkTime += tt1-tt0;
	
		long long int tot_kmol=0, cur_kmol=kmolecules.size();
		t0 = clock();
		MPI_Allreduce(&cur_kmol, &tot_kmol, 1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
//		if(MPIcontrol->rank==0)	printf("Proc %d: %llu Num of kmolecules is %llu(totNum is %llu)\n", MPIcontrol->rank, readRound, cur_kmol, tot_kmol);
		t1 = clock();
		commTime += (t1-t0);

		delete [] kmers;
		delete [] arcs;
		kmers=NULL;	arcs=NULL;
}

void  kmerGraph::dataShuttleBetweenGroupProcs(parameter *parameters, MPIEnviroment *MPIcontrol, char *& readBuf, long long int& localEnd, long long int readRound)
{
	int gProc = MPIcontrol->groupProc;
	//divide the sequences data for each group process, then distribute data to all its group processes. 
	if(MPIcontrol->rank%gProc==0)	
	{
		unsigned long long step;
		unsigned long long *start_pos = new unsigned long long [gProc+1];
        	step = (localEnd-1)/gProc + 1;
		for(int k=0;k<gProc;k++)
		{
			long long int i = step*k;
			while(i<localEnd && readBuf[i]!='>') i++; 
			assert(i==localEnd || readBuf[i]=='>');
			start_pos[k] = i;
		}
		start_pos[gProc]=localEnd;
		
		for(int k=1;k<gProc;k++)
		{
			unsigned long long int length = start_pos[k+1]-start_pos[k];
			MPI_Send(&length, 1, MPI_LONG_LONG_INT, MPIcontrol->rank+k, 0, MPI_COMM_WORLD);
			MPI_Send(readBuf+start_pos[k], length, MPI_CHAR, MPIcontrol->rank+k, 1, MPI_COMM_WORLD);
		}
		localEnd = start_pos[1];
		delete [] start_pos;
	}	
	else
	{
		unsigned long long length;
		MPI_Recv(&length, 1, MPI_LONG_LONG_INT, MPIcontrol->rank-MPIcontrol->rank%gProc, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		char *buffer = new char [length+2];	
		MPI_Recv(buffer, length, MPI_CHAR, MPIcontrol->rank-MPIcontrol->rank%gProc, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		readBuf = buffer;
		localEnd = length;
	}
}

int kmerGraph::constructKmerGraph(parameter *parameters, MPIEnviroment *MPIcontrol, bloomFilter *bFilter)
{
	clock_t IOTime=0, IOTimeRead=0, commTime=0, compTime=0, totTime=0, checkTime=0, t0, t1, t2;
	char message[1024];
	long long int readBases = min((long long int) MPIcontrol->nprocs/MPIcontrol->groupProc*1024*16, (long long int) BUF_SIZE);
	char          *readBuf = NULL; //new char [readBases*2+2];
	if(MPIcontrol->rank%MPIcontrol->groupProc==0)	readBuf = new char [readBases*2+2];	
	long long int readEnd=0, readStart=0, readLen=0, readRound=0;
	
	t0 = clock();
	MPIcontrol->File_open(parameters->fastaPath);
	
//	printf("Proc %d: %llu End\n", MPIcontrol->rank, readRound);
//	fflush(stdout);
	MPIcontrol->File_locate();
	t1 = clock();
	IOTime += t1-t0;
	totTime += t1-t0;
	int endTag = 1, totalTag = 0;
	
	do
	{
		t0 = clock();
		readLen = MPIcontrol->File_read(readBuf+readEnd, readBases);
		if(readBuf)	readBuf[readEnd+readLen] = 0;
		readEnd += readLen;		

		//losed the last sequence
//		if(MPIcontrol->datasize==0&&MPIcontrol->rank==MPIcontrol->nprocs-1)	{readBuf[readEnd]='>'; readBuf[++readEnd]=0;}	
		t1 = clock();
		IOTimeRead += t1-t0;
		if(DEBUG)
		{
			printf("Proc %d: %llu File_readLen=%llu, read_offset=%llu(%llu)\n", MPIcontrol->rank, readRound, readLen, MPIcontrol->read_offset-MPIcontrol->start_pos, MPIcontrol->datasize);
			fflush(stdout);
		}
		long long int localEnd = readEnd;
		
		if(readBuf)
		{
			while(localEnd>0 && readBuf[localEnd]!='>')	localEnd--;
			assert(readBuf[localEnd]=='>'||localEnd==0);
			readStart=localEnd;
		}
		//redistributed to groupProcs 	
		dataShuttleBetweenGroupProcs(parameters, MPIcontrol, readBuf, localEnd, readRound);

		extractKmerGraph(parameters, MPIcontrol, bFilter, readBuf, localEnd, readRound, commTime, checkTime);

		endTag=1, totalTag=0;
		if(MPIcontrol->datasize)	endTag=0;
		MPI_Allreduce(&endTag, &totalTag, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

		if(readBuf)	
		{
			for(long long int i=0; i<readEnd-readStart; i++)	readBuf[i] = readBuf[readStart+i];
			readBuf[readEnd-readStart] = 0;
			readEnd = readEnd-readStart;			
			readStart = 0;
		}

		t2 = clock();
		compTime += t2-t1;
		totTime  += t2-t0;
		readRound++;
		
		if(MPIcontrol->rank==0)
		{
			sprintf(message, "Proc %d: round %llu\n", MPIcontrol->rank, readRound);
			MPIcontrol->print(message);
		}
	}while(totalTag<MPIcontrol->nprocs);

	if(readBuf)
	{	
		for(long long int i=0; i<readEnd-readStart; i++)	readBuf[i] = readBuf[readStart+i];
		readBuf[readEnd-readStart] = 0;
	}
	readEnd = readEnd-readStart;			
	readStart = 0;
	long long int localEnd = readEnd;
	extractKmerGraph(parameters, MPIcontrol, bFilter, readBuf, localEnd, readRound, commTime, checkTime);
/*
	long long int tot_kmol=0, cur_kmol=0;
	
	for(unordered_map<unsigned long long, unsigned long long>::iterator I=kmolecules.begin(); I!=kmolecules.end(); I++)
//	for(auto I: kmolecules)
	{
		unsigned char *p = (unsigned char *) &(I->second);
		int t=0;
		for(int j=0;j<8;j++)	t += p[j];
		if(t>parameters->cutoffThreshold)	cur_kmol++;	
	}

	MPI_Allreduce(&cur_kmol, &tot_kmol, 1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
	printf("Proc %d: %llu Num of filtered kmolecules is %llu(totNum is %llu)\n", MPIcontrol->rank, readRound, cur_kmol, tot_kmol);
*/
	clock_t avgIOTime, avgIOTimeRead, avgCompTime, avgCommTime, avgTotTime, avgCheckTime;
	MPI_Allreduce(&IOTime,   &avgIOTime,   1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&IOTimeRead, &avgIOTimeRead, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&compTime, &avgCompTime, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&commTime, &avgCommTime, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&totTime,  &avgTotTime,  1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&checkTime,  &avgCheckTime,  1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
	avgIOTime   = avgIOTime/MPIcontrol->nprocs;
	avgIOTimeRead = avgIOTimeRead/MPIcontrol->nprocs;
	avgCommTime = avgCommTime/MPIcontrol->nprocs;
	avgCompTime = avgCompTime/MPIcontrol->nprocs;
	avgTotTime  = avgTotTime/MPIcontrol->nprocs;
	avgCheckTime = avgCheckTime/MPIcontrol->nprocs;
	
	long long int tot_kmol=0, cur_kmol=kmolecules.size();
	MPI_Allreduce(&cur_kmol, &tot_kmol, 1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);

	if(DEBUG)
	{	
		sprintf(message, "Proc %d: %llu Num of kmolecules is %llu(totNum is %llu)", MPIcontrol->rank, readRound, cur_kmol, tot_kmol);
		MPIcontrol->print(message);
	}

	if(MPIcontrol->rank==0)	
	{
		printf("IOTime=%.5fs IOTimeRead=%.5fs compTime=%.5fs commTime=%.5fs totTime=%.5fs checkTime=%.5fs\n", avgIOTime/(float)CLOCKS_PER_SEC, avgIOTimeRead/(float)CLOCKS_PER_SEC, avgCompTime/(float)CLOCKS_PER_SEC, avgCommTime/(float)CLOCKS_PER_SEC, avgTotTime/(float)CLOCKS_PER_SEC, avgCheckTime/(float)CLOCKS_PER_SEC);
	}

	MPIcontrol->File_close();
	if(readBuf)	delete [] readBuf;
	readBuf = NULL;
	return 1;
}

void kmerGraph::printKmerFreq(parameter *parameters, MPIEnviroment *MPIcontrol)
{
	long long int loss = 0;
	long long int missFreq = 0;
	long long int delta = 0;
	
	long long int kmernum=tmp.size(), validkmernum=0, bloomkmer=kmolecules.size();
	for(unordered_map<unsigned long long, unsigned long long>::iterator I=tmp.begin(); I!=tmp.end(); I++)
	{
		if(I->second>=parameters->cutoffThreshold)
		{
			if(kmolecules.find(I->first)==kmolecules.end())	loss++;
			else if(kmolecules[I->first]!=I->second)	missFreq++, delta += abs((int) (kmolecules[I->first]-I->second));
			validkmernum++;
		}	
	}

	printf("Proc %d: tmpkmernum = %lu, filterkmernum = %lu, bloomkmer=%lu\n", MPIcontrol->rank, kmernum, validkmernum, bloomkmer);
	printf("loss = %lu, miss = %lu, delta=%lu\n", loss, missFreq, delta/missFreq);
	long long int totloss=0, totmiss=0, totdelta=0, totkmernum=0, totbloomkmer=0, allkmers=0;	
	MPI_Reduce(&loss, &totloss, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);	
	MPI_Reduce(&missFreq, &totmiss, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);	
	MPI_Reduce(&delta, &totdelta, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);	
	MPI_Reduce(&kmernum, &allkmers, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);	
	MPI_Reduce(&validkmernum, &totkmernum, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);	
	MPI_Reduce(&bloomkmer, &totbloomkmer, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);	
	if(MPIcontrol->rank==0)	printf("totloss=%lu,  totmiss=%lu, totdelta=%.5f\n", totloss, totmiss, totdelta*1.0/totmiss);
	if(MPIcontrol->rank==0)	printf("allkmers=%lu, totfilterkmernum=%lu, totbloomkmernum=%lu\n", allkmers, totkmernum, totbloomkmer);


	ofstream fs(parameters->kmerCount);
      	if(!fs.is_open())  
       	{ 
		cout<<"Error opening file"<<endl; 
		exit(1); 
	}  
	int freqLen = 10000;
	int i=0;
	long long int a[freqLen];
	memset(a, 0, sizeof(long long int)*freqLen);
	for(unordered_map<unsigned long long, unsigned long long>::iterator I=kmolecules.begin(); I!=kmolecules.end(); I++)
	{
/*		unsigned char *p = (unsigned char *) &(I->second);
		int t=0;
		for(int j=0;j<8;j++)	t += p[j];
*/
		unsigned long long t = (unsigned long long) (I->second);
		if(t<freqLen)		a[t]++;	
		i++;	
	}
/*
	long long int sum=0;
        for(int i=1; i<freqLen; i++) sum+=a[i];
        a[1]=sum-a[1];
        for(int i=2;i<freqLen-1;i++) a[i]=a[i-1]-a[i];   //now the real value shift to left 
	for(int i=freqLen-1;i>1;i--) a[i]=a[i-1];
	a[1]=sum;
*/
	long long int b[freqLen];
	MPI_Reduce(&a, &b, freqLen, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);	
	if(MPIcontrol->rank==0)	{for(int i=1;i<freqLen;i++)	if(b[i])	fs<<i<<" "<<b[i]<<endl;}
	fs.close();
}


//Please Modify this function to print the frequence of all kmers.  jintaomeng
void kmerGraph::printKmerCounting(parameter *parameters, MPIEnviroment *MPIcontrol)
{
        ofstream fs(parameters->kmerGraph);
        if(!fs.is_open())
        {
                cout<<"Error opening file"<<endl;
                exit(1);
        }

	int lineLen = parameters->hashLength+1+4*8+5;
	int totLen = kmolecules.size()*lineLen;
	
	char *buf = new char [totLen];
	int i=0;

	for(unordered_map<unsigned long long, unsigned long long>::iterator I=kmolecules.begin(); I!=kmolecules.end(); I++)
	{
		sprintf(buf+i*lineLen, "%s ", longLongToString(I->first, parameters).c_str());
		unsigned char *p = (unsigned char *) &I->second;
		int t=0;
		/*for(int j=0;j<8;j++)	
		{
			t += p[j];
			sprintf(buf+i*lineLen+parameters->hashLength+1+j*4, "%3d ", (int) p[j]);
		}*/
		sprintf(buf+i*lineLen+parameters->hashLength+1+4*8, "%4d\n", t);	
		i++;	
	}
	buf[totLen] = 0;
	MPIcontrol->File_write(parameters->kmerGraph, buf, totLen);
}

int main(int argc, char *argv[])
{
	MPIEnviroment MPIcontrol;
	MPIcontrol.init(argc, argv);
	
	parameter parameters;
	parameters.getParameters(argc, argv, &MPIcontrol);	

	FILE *FP=NULL;
    	FP = fopen(parameters.LogPath, "a");
    	clock_t t0=0, t1=0;
    	if(MPIcontrol.rank==0)
	{
		t0 = clock();
    //		FP = fopen(parameters.LogPath, "w");
	}

	kmerGraph mygraph;
	bloomFilter bFilter((unsigned long long) parameters.refSize*4/MPIcontrol.nprocs, 3);
	if(DEBUG)   printf("Proc %d: MPIcontrol->groupProc %d\n", MPIcontrol.rank, MPIcontrol.groupProc);
	
	mygraph.constructKmerGraph(&parameters, &MPIcontrol, &bFilter);
	mygraph.printKmerFreq(&parameters, &MPIcontrol);
//	fflush(stdout);	
//	mygraph.printKmerCounting(&parameters, &MPIcontrol);	
    
    	if(MPIcontrol.rank==0)
   	{
		char message[1024];
		sprintf(message,"Time spend in I/O part in MPIcontrol %.5f\n", MPIcontrol.elapsedTime/(float)CLOCKS_PER_SEC);
 		MPIcontrol.print(message, FP);	
		sprintf(message, "Time spend in locate function %.5f\n", MPIcontrol.locateTime/(float)CLOCKS_PER_SEC);
		MPIcontrol.print(message,FP);
		sprintf(message,"Time spend in read function %.5f\n", MPIcontrol.readTime/(float)CLOCKS_PER_SEC);
		MPIcontrol.print(message,FP);
    	}

	if(MPIcontrol.rank==0)
	{
		t1 = clock();
		char message[1024];
		sprintf(message,"Time spend in I/O part in MPIcontrol %.5f\n", (t1-t0)/(float)CLOCKS_PER_SEC);
 		MPIcontrol.print(message, FP);	
		printf("Time spend in I/O part in MPIcontrol %.5f\n", (t1-t0)/(float)CLOCKS_PER_SEC);
		fclose(FP);
	}

	MPIcontrol.finalize();
	return 0;
}
