/* File:     mpi_prime.c
 * Purpose:  Find all primes less than or equal to an input value.
 *           This version doesn't bother checking even ints.
 *
 * Input:    n:  integer >= 2 (from command line)
 * Output:   Sorted list of primes between 2 and n,
 *
 * Compile:  gcc -g -Wall -o primes1 primes1.c -lm
 * Usage:    ./primes1 <n>
 *              n:  max int to test for primality
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mpi.h"

//#define DEBUG

void Merge(int** all_primes, int *all_count, int* received_primes,
		int received_count, int** temp);

/* Helper Functions */
int Is_prime(int i);
int Sum(int A[], int n);
int Get_n(int argc, char* argv[]);
int Smallest_power_two(int p);
void Usage(char prog[]);
void Print_vector(char* title, int y[], int m, int my_rank);

int main(int argc, char* argv[]) {
	int n, i, p, my_rank;
	int local_count, my_count, received_count, sum_count;
	int divisor, proc_diff, partner, done, i_send;
	int *local_primes, *received_primes, *all_primes, *all_counts, *temp;
	MPI_Comm comm;
	char title[100];

	MPI_Init(&argc, &argv);
	comm = MPI_COMM_WORLD;
	MPI_Comm_size(comm, &p);
	MPI_Comm_rank(comm, &my_rank);

	if (my_rank == 0)
		n = Get_n(argc, argv);

	MPI_Bcast(&n, 1, MPI_INT, 0, comm);

	if (n <= 1) {
		if (my_rank == 0)
			Usage(argv[0]);
		return 0;
	} else if (n == 2) {
		if (my_rank == 0)
			printf("The primes < 2 are: 2");
		return 0;
	}

	local_primes = malloc(n / (2 * p) + 2 * sizeof(int));

	// logic for getting 'my' primes
	for (i = 2 * my_rank + 3, local_count = 0; i < n; i += 2 * p) {
		if (Is_prime(i)) {
			local_primes[local_count] = i;
			local_count++;
		}
	}

	// get the amount of primes
	all_counts = malloc(p * sizeof(int));
	MPI_Allgather(&local_count, 1, MPI_INT, all_counts, 1, MPI_INT, comm);

	sum_count = Sum(all_counts, p);

	// make three arrays that are big enough to store all primes from all processes
	all_primes = malloc(sum_count * sizeof(int));
	received_primes = malloc(sum_count * sizeof(int));
	temp = malloc(sum_count * sizeof(int));

	// copy over my values to beginning of all_primes array
	for (i = 0; i < local_count; i++)
		all_primes[i] = local_primes[i];

#	ifdef DEBUG
	MPI_Barrier(comm);
	sprintf(title,"all_primes: only my primes");
	Print_vector(title, all_primes, local_count, my_rank);
	fflush(stdout);
#	endif

	divisor = 2;
	proc_diff = 1;
	done = 0;
	my_count = local_count;

	while (!done && divisor <= Smallest_power_two(p)) {
		i_send = my_rank % divisor;
		if (i_send) {
			partner = my_rank - proc_diff;
#			ifdef DEBUG
			sprintf(title,"my_count %d | proc_diff %d | snd all_primes to %d: ", my_count, proc_diff, partner);
			Print_vector(title, all_primes, my_count, my_rank);
			fflush(stdout);
#			endif
			// send all my primes over to my partner
			MPI_Send(all_primes, my_count, MPI_INT, partner, 0, comm);
			done = 1;
		} else { /* I'm receiving */
			partner = my_rank + proc_diff;
			if (partner < p) {
				// calculate how many primes i should be receiving from partner
				received_count = Sum(&(all_counts[partner]), (partner
						+ proc_diff < p) ? proc_diff : p - partner);
				// receive primes from my partner
				MPI_Recv(received_primes, received_count, MPI_INT, partner, 0,
						comm, MPI_STATUS_IGNORE);
#				ifdef DEBUG
				sprintf(title,"received_count %d | proc_diff %d | rcv all_primes from %d: ", received_count, proc_diff, partner);
				Print_vector(title, received_primes, received_count, my_rank);
				fflush(stdout);
#				endif

				// merge my primes with the ones i just received, updating all_primes
				Merge(&all_primes, &my_count, received_primes, received_count,
						&temp);
#				ifdef DEBUG
				sprintf(title,"my_count %d | proc_diff %d | merged: ", my_count, proc_diff);
				Print_vector(title, all_primes, my_count, my_rank);
				fflush(stdout);
#				endif
			}
			divisor *= 2;
			proc_diff *= 2;
		}
	}

	// process 0 will print out the final list
	if (my_rank == 0) {
		sprintf(title,"The primes < %d are: 2", n);
		Print_vector(title, all_primes, my_count, my_rank);
	}

	MPI_Finalize();
	free(local_primes);
	free(all_primes);
	free(received_primes);
	free(temp);
	return 0;
} /* main */

/*-------------------------------------------------------------------
 * Function:   Is_prime
 * Purpose:    Determine whether the argument is prime
 * Input arg:  i
 * Return val: true (nonzero) if arg is prime, false (zero) otherwise
 */
int Is_prime(int i) {
	int j;
	int sqrt_i;

	sqrt_i = sqrt(i);
	for (j = 2; j <= sqrt_i; j++)
		if (i % j == 0)
			return 0;
	return 1;
} /* Is_prime */

/*-------------------------------------------------------------------
 * Function:   Merge
 * Purpose:    Merge the contents of the all_primes with received_primes
 * Input args:
 * 	  received_primes: 	the array to be merged with all_primes
 *    received_count:  	the number of elements in received_primes
 * In/Out arg:
 *    all_primes:		array to which all values end up in after merge
 *    all_count:  		the number of elements in all_primes, updated to new size
 *    temp:				temporarily stores merged array. swapped with all_primes
 *    					to reduce malloc calls.
 */
void Merge(int** all_primes, int *all_count, int* received_primes,
		int received_count, int** temp) {
	int all_i, received_i, temp_i, all_count_t;
	int *all_p, *temp_p;

	// set temporary array pointers
	all_p = *all_primes;
	temp_p = *temp;

	// save all_primes count and set temp primes count to total
	all_count_t = *all_count;
	*all_count = all_count_t + received_count;

	all_i = received_i = temp_i = 0;

	// merge lists
	while (all_i < all_count_t && received_i < received_count) {
		if (all_p[all_i] <= received_primes[received_i]) {
			temp_p[temp_i] = all_p[all_i];
			all_i++;
		} else {
			temp_p[temp_i] = received_primes[received_i];
			received_i++;
		}
		temp_i++;
	}

	// fill in the tail end
	if (all_i >= all_count_t)
		for (; temp_i < *all_count; temp_i++, received_i++)
			temp_p[temp_i] = received_primes[received_i];
	else
		for (; temp_i < *all_count; temp_i++, all_i++)
			temp_p[temp_i] = all_p[all_i];

	// set all_primes to point to the answer array temp
	*all_primes = temp_p;
	// set temp to point to original all_primes (for re-use)
	*temp = all_p;

} /* Merge */

/*-------------------------------------------------------------------
 * Function:  Usage
 * Purpose:   Print a brief message explaining how the program is run.
 *            Then quit.
 * In arg:    prog:  name of the executable
 */
void Usage(char prog[]) {
	fprintf(stderr, "usage: %s <n>\n", prog);
	fprintf(stderr, "   n = max integer to test for primality\n");
} /* Usage */

/*-------------------------------------------------------------------
 * Function:    Get_n
 * Purpose:     Get the input value n
 * Input args:  argc:  number of command line args
 *              argv:  array of command line args
 */
int Get_n(int argc, char* argv[]) {
	int n;

	if (argc != 2)
		Usage(argv[0]);

	/* Convert string to int */
	n = strtol(argv[1], NULL, 10);

	return n;
} /* Get_n */

/*------------------------------------------------------------------
 * Function:    Smallest_power_two
 * Purpose:     Find the largest power of two smaller than an integer
 * In args:     integer n
 * Return:		the power of two
 */
int Smallest_power_two(int n) {
	int power_two = 1;

	while (power_two < n)
		power_two *= 2;

	return power_two;
} /* Smallest_power_two */

/*------------------------------------------------------------------
 * Function:    Sum
 * Purpose:     Add up all elements of integer vector
 * In args:     array A, size of array n
 * Return:		sum
 */
int Sum(int A[], int n) {
	int i, sum;

	for (i = 0, sum = 0; i < n; i++)
		sum += A[i];

	return sum;
} /* Sum */

/*------------------------------------------------------------------
 * Function:    Print_vector
 * Purpose:     Print a vector
 * In args:     title, A, n, my_rank
 */
void Print_vector(char* title, int A[], int n, int my_rank) {
	int i;
	char buffer[10000];

	sprintf(buffer,"%d> %s ", my_rank, title);
	for (i = 0; i < n; i++)
		sprintf(buffer, "%s %d ",buffer, A[i]);
	sprintf(buffer,"%s \n", buffer);
	printf("%s", buffer);
} /* Print_vector */
