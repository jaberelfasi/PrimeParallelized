#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <mpi.h>
#include <math.h>

int CountPrimes(unsigned long long int n0, unsigned long long int n1, unsigned long int *basePrimes, int nBasePrimes)
{
	int nPrimes = 0;

		int i;
		unsigned long long *working;
		working = (unsigned long long int*)malloc( (n1-n0)*sizeof(unsigned long long int) );
		memset( working, 0, (n1-n0)*sizeof(unsigned long long int) );

		for ( i = 0; i < nBasePrimes; i++ )
		{
			// find starting point
			long long int i0;
			if ( (n0 % basePrimes[i]) == 0 )
				i0 = n0;
			else
				i0 = (n0 / basePrimes[i] + 1)*basePrimes[i];
			while ( i0 < n1 )
			{
				working[i0 - n0] = 1;
				i0 += basePrimes[i];
			}
		}
		for ( i = 0; i < (n1-n0); i++ )
			if ( !working[i] )
				nPrimes++;

		free( working );

		return nPrimes;
}

int main( void )
{

	//printf("here starts the mpi part");

	unsigned long long int nPrimes=0;

	#define MAX_PRIME 100000 // we are counting primes up to 1e10, so we need a base primes list up to sqrt(1e10) = 1e5
	#define SQRT_MAX_PRIME 317 // the largest number we will need to cross out multiples of for our base prime list

	unsigned long int *working = (unsigned long int *)malloc( (MAX_PRIME+1)*sizeof(unsigned long int) );
	unsigned long int *primes = (unsigned long int *)malloc( MAX_PRIME*sizeof(unsigned long int) );

	int rank, size;


	//unsigned long long int ip;
	unsigned long long int delta =100000000/1000;

	unsigned long long int pmax =10000000000;

	//printf("here starts the mpi part");

	MPI_Init(NULL, NULL);
	double initialT = MPI_Wtime();
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);


	if(rank==0)
	{

			// construct the base prime list (primes up to 1e5) using the sieve of Eratosthenes
			int i, j;
			for ( i = 0; i < MAX_PRIME; i++ )
			{
				working[i] = 1;
			}
			working[0] = working[1] = 0;

			// sieve
			for ( i = 2; i < SQRT_MAX_PRIME; i++ )
			{
				if ( working[i] )
				{
					for ( j = 2*i; j < MAX_PRIME; j += i )
					{
						working[j] = 0;
					}
				}
			}
			nPrimes = 0;

			for ( i = 0; i < MAX_PRIME; i++ )
				if ( working[i] )
					primes[nPrimes++] = i;
			free( working );

	}


	//this line broad casts the nPrimes to all of the ranks
	MPI_Bcast(&nPrimes, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);

	//this line broad casts the prime array to all of the ranks
	MPI_Bcast(primes, nPrimes, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);



	int nBase = nPrimes;
	unsigned long long int getPrimes=0;

	int IPinitial = primes[nPrimes-1]+1;
	unsigned long long int IP;

			int i = rank;
			float x, y;
			if(rank==0)
			{
			x=ceil(pmax/delta); //gives us how many number of IPs will be in the computation.
			while(((int)x%size)!=0)
			{
				delta++;
				x=ceil(pmax/delta);
			}
			//y=ceil(x/size); //gives us how many iterations of the following IPs in each rank.
			}
			MPI_Bcast(&x, 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Bcast(&y, 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Bcast(&delta, 1, MPI_INT, 0, MPI_COMM_WORLD);

			IP = IPinitial+(i*delta);
			unsigned long long int IpElements=0;//this is the illustrate how many elements from IP to next or to pmax
			unsigned long long int ChunkCount=1;

			unsigned long long int next;
			while(IP<pmax)
			{
				next=IP+delta;
				if(next>pmax)
					next = pmax;
				printf("counting from %lld\n", IP);
				getPrimes +=CountPrimes( IP, next, primes, nBase );
				IpElements=next-IP;
				ChunkCount++;
				IP=IP+delta*size;
			}



	unsigned long long int sumPrimes=0;

	MPI_Reduce(&getPrimes, &sumPrimes, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

	if(rank==0)
	{
		printf("prime counting function(%lld) = %lld\n",pmax,sumPrimes+nPrimes);
		double endT = MPI_Wtime();
		double totalT=endT-initialT;
		printf("total time consumed is: %f\n",totalT);
	}


	// the code below is to illustrate the load balancing
	MPI_Send(&IpElements,1,MPI_UNSIGNED_LONG,0,0,MPI_COMM_WORLD);
	MPI_Send(&ChunkCount, 1, MPI_UNSIGNED_LONG, 0, 0, MPI_COMM_WORLD);

	free(primes);

	if(rank==0)
	{

		for(int i=0; i<size; i++)
		{
			MPI_Recv(&IpElements, 1, MPI_UNSIGNED_LONG, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Recv(&ChunkCount, 1, MPI_UNSIGNED_LONG, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			printf("in the last chunk of rank %d from IP to next or pmax there are %lld elements\n", i, IpElements);
			printf("rank %d processed %lld number of chunks\n", i, ChunkCount);
		}

	}


	MPI_Finalize();
	return 0;
}
