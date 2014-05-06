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
#define MAX_SIEVE_PRIME (1000000)

int main(int argc, char **argv)
{
    mpz_t min, rop, sievesize;
    unsigned long r;
    unsigned long loop_count = 0;
    FILE *fp = stdout;
    bool *sieve;
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

    sieve = (bool *) malloc(sizeof(bool) * SIEVESIZE);
    memset(sieve, true, SIEVESIZE);
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
#pragma omp parallel
    {
	mpz_t candidate_plus;
	mpz_init(candidate_plus);
#pragma omp for schedule(static, 16)
	for (int i = 0; i < prime_size; i++) {
	    int p = primes[i];
	    int inv = modinv_i(PRIMORAL, p);
	    //printf("%d %d %d %d\n", i, p, PRIMORAL, inv);
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

    while (1) {
	//step 1
	unsigned int sieve_loc_index = 0;
	//printf("Sieving\n");

	gettimeofday(&start_t, NULL);
#pragma omp parallel
	{
	    int tid = omp_get_thread_num();
	    int num = omp_get_num_threads();
	    int chunksize = (SIEVESIZE + num - 1) / num;
	    int start = chunksize * tid;
	    int end = start + chunksize;
	    if (end > SIEVESIZE)
		end = SIEVESIZE;

	    //printf("thread:%d start:%d end:%d\n", tid,start,end);
	    for (int i = 0; i < prime_size; i++) {
		for (int j = 0; j < 6; j++) {
		    //o = sieve_loc[sieve_loc_index];
		    int o = sieve_loc[8 * i + j];
		    int t_start = o;
		    int jump;
		    if (start > o) {
			jump = (start + primes[i] - 1 - o) / primes[i];
		    } else {
			jump = 0;
		    }

		    o += jump * primes[i];
		    //while(o<start) o+=primes[i];
		    while (o < end) {
			sieve[o] = false;
			o += primes[i];
		    }
		    //if (tid == num - 1)
		    //sieve_loc[8 * i + j] = o - SIEVESIZE;
		    //sieve_loc_index += 1;
		}
	    }
#pragma omp barriar
//update sieve_loc
#pragma omp for
	    for (int i = 0; i < prime_size; i++) {
		for (int j = 0; j < 6; j++) {
		    int o = sieve_loc[8 * i + j];
		    int t_start = o;
		    int jump;
		    if (SIEVESIZE > o) {
			jump = (SIEVESIZE + primes[i] - 1 - o) / primes[i];
		    } else {
			jump = 0;
		    }
		    o += jump * primes[i];
		    sieve_loc[8 * i + j] = o - SIEVESIZE;

		}
	    }
	}
	gettimeofday(&end_t, NULL);
	duration =
	    (end_t.tv_sec - start_t.tv_sec) * 1E6 + (end_t.tv_usec -
						     start_t.tv_usec);
	printf("running time1:%f\n", duration);

	gettimeofday(&start_t, NULL);
/*for(int i=0;i<SIEVESIZE;i++)
  if(sieve[i])
    printf("1\n");
  else
    printf("0\n");*/
/*for(int i=0;i<prime_size;i++){
  for(int j=0;j<6;j++)
    printf("%d:%d\n",primes[i], sieve_loc[8*i+j]);
}
break;*/
	//printf("Testing Candidates\n");
#pragma omp parallel
	{
	    mpz_t tmp;
	    mpz_t candidate;
	    mpz_init(tmp);
	    mpz_init(candidate);
#pragma omp for
	    for (int i = 0; i < SIEVESIZE; i++) {
		if (sieve[i]) {
		    mpz_set_ui(tmp, i);
		    mpz_addmul_ui(tmp, sievesize, loop_count);
		    mpz_set(candidate, rop);
		    mpz_addmul_ui(candidate, tmp, PRIMORAL);
		    if (is_fermat_valid(candidate)) {
			mpz_out_str(fp, 16, candidate);
			exit(0);
		    }
		}
		sieve[i] = true;
	    }
	    mpz_clear(tmp);
	    mpz_clear(candidate);
	}
	gettimeofday(&end_t, NULL);
	duration =
	    (end_t.tv_sec - start_t.tv_sec) * 1E6 + (end_t.tv_usec -
						     start_t.tv_usec);
	printf("running time2:%f\n", duration);

	loop_count++;
    }

    return 0;
}
