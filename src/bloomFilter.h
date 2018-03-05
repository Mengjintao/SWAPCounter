#include <stdio.h>
#include <assert.h>
#include <fcntl.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <iostream>
#include <vector>
#include <string>

using namespace std;

#define MDIVN 8
#define FunctionNum 4

class bloomFilter
{
private:
	unsigned long long *bitSet;
	unsigned long long size;
	int bits;
	int funNum;
	vector<unsigned long long> factorArray; // = {19,23,29,31, 37,39,257,263,269,271,277,281,283,293,307,311};

	bloomFilter()
	{
		bits=1;
		funNum=3;
		size=0;
		bitSet=NULL;
	}

	unsigned long long hashFun(unsigned long long a, int index)
	{
//		printf("hashFun a, index =%llu %d %d\n", a, index, factorArray.size());
//		fflush(stdout); 
         	unsigned long long factor = factorArray[index];
        	unsigned long long sum = 0;
         	for(unsigned i = 0; i<8; i++){
                 	sum = (sum * factor + (a & 0xFFFFFFFF))%size;
			if(a==0)	break;
			a = (a>>8);
         	}
         	return sum%size;
	}

public:	
	bloomFilter(unsigned long long sz=0, int bitsPerNode=1)
	{
		factorArray.push_back(19);
		factorArray.push_back(23);
		factorArray.push_back(29);
		factorArray.push_back(31);
		factorArray.push_back(37);
		factorArray.push_back(39);
		factorArray.push_back(257);
		factorArray.push_back(263);
		factorArray.push_back(269);
		factorArray.push_back(271);
		factorArray.push_back(277);
		factorArray.push_back(281);
		factorArray.push_back(283);
		factorArray.push_back(293);
		factorArray.push_back(307);
		factorArray.push_back(311);

		funNum = FunctionNum;
		size   = sz*MDIVN;
		bits   = bitsPerNode;
		bitSet = new unsigned long long [(bits*size)/64+1];
		printf("BloomFilter Set: funNum=%d size=%llu bitsPerNode=%d\n", funNum, size, bits); 
		for(unsigned long long i=0;i<(bits*size)/64+1;i++)	bitSet[i]=0;
	}

	~bloomFilter()
	{
		delete [] bitSet;
	}
	
	int bloomCheck(unsigned long long a)
	{
		if(size==0)	return 0;
		unsigned long long minC = (1<<30);
		for(int i=0;i<funNum;i++)
		{
			unsigned long long pos = hashFun(a, i);						
			pos *= bits;
			unsigned long long sum=0;
			for(int j=0;j<bits;j++)
			{
				unsigned long long lpos = pos+j;
				sum = sum*2 + ((bitSet[lpos/64]>>(lpos%64))&1);
			}	
			minC = min(minC, sum);
		}
		return (int) minC;
	}

	int bloomAdd(unsigned long long a)
	{
//		printf("bloomAdd: %llu\n", a);
//		fflush(stdout);	
		if(size==0)	return 0;
		unsigned long long minC = (1ull<<30);
		for(int i=0; i<funNum; i++)
		{
			unsigned long long pos = hashFun(a, i);		
			pos *= bits;
			unsigned long long sum=0;
			for(int j=0; j<bits; j++)
			{
				unsigned long long k = pos+j;
				sum = sum*2 + ((bitSet[k/64]>>(k%64))&1);
			}
			
			if(sum+1<(1ull<<bits))	sum++; 
			for(int j=0;j<bits;j++)
			{
				unsigned long long k = pos+j;
				bitSet[k/64] &= ~((unsigned long long) 1 << (k%64)); 
				bitSet[k/64] |= ( ((sum>>(bits-1-j))&1) << (k%64));
			}
			minC = min(minC, sum);
		}
		return (int) minC;
	}		

	void bloomPrint()
	{
		for(unsigned long long i=0;i<size;i++)		printf("%lld", i%10);
		printf("\n");
		for(unsigned long long i=0;i<size;i++)	
		{
			unsigned long long sum=0;
			for(int j=0;j<bits;j++)	
			{
				unsigned long long k = i*bits+j; 
				sum = sum*2+ ((bitSet[k/64]>>(k%64))&1);
			}
			printf("%llu", sum);
		}
		printf("\n");
	}
	
};

/*
int main()
{
	bloomFilter bFilter(10, 3);	
	unsigned long long a[10] = {5,3, 1, 1, 2, 6, 4, 4,5, 9};
	for(int i=0;i<10;i++)
	{
		bFilter.bloomPrint();
		cout<<a[i]<<" "<<bFilter.bloomAdd(a[i])<<endl;	
	}
	return 1;
}
*/
