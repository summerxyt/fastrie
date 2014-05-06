#include "miner.h"
#include "fermat.h"
#include <string.h>
#include <stdlib.h>
#include <primesieve.h>
#include <omp.h>
#include <sys/time.h>

//2*3*5*7
#define PRIMORAL (2*3*5*7)
#define SIEVESIZE (1024*512)
//#define SIEVESIZE (1024*512)
//#define MAX_SIEVE_PRIME (1000000000)
#define MAX_SIEVE_PRIME (1000000)

int main(int argc, char **argv)
{
    mpz_t min, rop, sievesize;
    unsigned long r;
    unsigned long loop_count = 0;
    FILE *fp = stdout;
    bool *sieve[128];
    int *primes, *primes_start, *sieve_loc;
    size_t prime_size, loc_size;

    mpz_init(min);
    mpz_init(rop);
    mpz_init(sievesize);

    mpz_set_str(min, argv[1], 0);
    mpz_set_ui(sievesize, SIEVESIZE);

    r = mpz_fdiv_ui(min, PRIMORAL);
    r = (97 - r + PRIMORAL) % PRIMORAL;

    mpz_add_ui(rop, min, r);

    primes_start =
	(int *) primesieve_generate_primes(0, MAX_SIEVE_PRIME, &prime_size,
					   INT_PRIMES);
    for (primes = primes_start; *primes <= 7; primes++);
    prime_size -= (primes - primes_start);
    //loc_size = 6 * prime_size;
    loc_size = 8 * prime_size;	//padding
    sieve_loc = (int *) malloc(sizeof(int) * loc_size);

    struct timeval start_t;
    struct timeval end_t;
    gettimeofday(&start_t, NULL);
    {
	mpz_t candidate_plus;
	mpz_init(candidate_plus);
	int tid = omp_get_thread_num();
	sieve[tid] = (bool *) malloc(sizeof(bool) * SIEVESIZE);
	if(sieve[tid]==NULL)
	  printf("%d nil\n", tid);
	memset(sieve[tid], true, SIEVESIZE);
	for (int i = 0; i < prime_size; i++) {
	    //printf("%d:%d\n",tid,i);
	    int p = primes[i];
	    int inv = modinv_i(PRIMORAL, p);
	    for (int j = 0; j < 6; j++) {
		mpz_add_ui(candidate_plus, rop, offsets[j]);
		unsigned long cmodp = mpz_fdiv_ui(candidate_plus, p);
		int index = ((p - cmodp) * inv) % p;
		sieve_loc[i * 8 + j] = index;
	    }
	}
	mpz_clear(candidate_plus);
    }
    gettimeofday(&end_t, NULL);
    float duration =
	(end_t.tv_sec - start_t.tv_sec) * 1E6 + (end_t.tv_usec -
						 start_t.tv_usec);
    printf("running time:%f\n", duration);

	    mpz_t tmp;
	    mpz_t candidate;
	    mpz_init(tmp);
	    mpz_init(candidate);
    while (1) {
	//step 1
	//printf("Sieving\n");

	{
	gettimeofday(&start_t, NULL);
	    int tid = omp_get_thread_num();
	    for (int i = 0; i < prime_size; i++) {
		    //printf("thread %d:prime%d\n", tid, primes[i]);
		for (int j = 0; j < 6; j++) {
		    //o = sieve_loc[sieve_loc_index];
		    unsigned int o = sieve_loc[8 * i + j];
		    while (o < SIEVESIZE) {
			sieve[tid][o] = false;
			o += primes[i];
		    }
		    sieve_loc[8 * i + j] = o - SIEVESIZE;
		}
	    }
    gettimeofday(&end_t, NULL);
    float duration =
	(end_t.tv_sec - start_t.tv_sec) * 1E6 + (end_t.tv_usec -
						 start_t.tv_usec);
    printf("running time:%f\n", duration);
	//}
	//printf("Testing Candidates\n");
	//#pragma omp parallel
	//{
	gettimeofday(&start_t, NULL);
	    int num = omp_get_num_threads();
	    for (int i = 0; i < SIEVESIZE; i++) {
		bool flag = true;
		for (int j = 0; j < num; j++) {
		    flag = flag && sieve[j][i];
		    sieve[j][i] = true;
		}
		if (flag) {
		    mpz_set_ui(tmp, i);
		    mpz_addmul_ui(tmp, sievesize, loop_count);
		    mpz_set(candidate, rop);
		    mpz_addmul_ui(candidate, tmp, PRIMORAL);
		    if (is_fermat_valid(candidate)) {
			mpz_out_str(fp, 16, candidate);
			exit(0);
		    }
		}
	    }
    gettimeofday(&end_t, NULL);
    duration =
	(end_t.tv_sec - start_t.tv_sec) * 1E6 + (end_t.tv_usec -
						 start_t.tv_usec);
    printf("running time:%f\n", duration);
	//}

	}
	loop_count++;
    }
	    mpz_clear(tmp);
	    mpz_clear(candidate);

    return 0;
}
