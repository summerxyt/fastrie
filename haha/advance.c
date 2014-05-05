#include "miner.h"
#include "fermat.h"
#include <string.h>
#include <stdlib.h>
#include <primesieve.h>

//2*3*5*7
#define PRIMORAL (2*3*5*7)
//#define SIEVESIZE (1024*1024)
#define SIEVESIZE (1024*512)
#define MAX_SIEVE_PRIME (1000000)
//#define NUM_PRIMES 78494
//#define NUM_LOCS (NUM_PRIMES * 6)

//static const unsigned int offsets[6] = {0, 4, 6, 10, 12, 16};

int main(int argc, char **argv){
  mpz_t min, max, rop, sievesize;
  unsigned long r;
  unsigned long loop_count=0;
  FILE* fp = stdout;
  bool *sieve;
  int *primes, *primes_start, *sieve_loc;
  size_t prime_size, loc_size;

  mpz_init(min);
  mpz_init(max);
  mpz_init(rop);
  mpz_init(sievesize);

  mpz_set_str(min, argv[1], 0);
  mpz_set_ui(sievesize, SIEVESIZE);
  //mpz_set_str(max, argv[2], 0);

  printf("here\n");

  r = mpz_fdiv_ui(min, PRIMORAL);
  r = (97 - r + PRIMORAL) % PRIMORAL;

  mpz_add_ui(rop, min, r);
  
  sieve = (bool*)malloc(sizeof(bool) * SIEVESIZE);
  memset(sieve, true, SIEVESIZE);
  primes_start = (int *) primesieve_generate_primes(0, MAX_SIEVE_PRIME, &prime_size, INT_PRIMES);
  for(primes=primes_start;*primes<=7;primes++);
  prime_size -= (primes - primes_start);
  loc_size = 6 * prime_size;
  sieve_loc = (int *)malloc(sizeof(int) * loc_size);

  mpz_t candidate_plus;
  mpz_init(candidate_plus);
  //#pragma omp parallel for
  for(int i=0;i<prime_size;i++){
    //if(primes[i] <= 7)
    //  continue;
    int p = primes[i];
    int inv = modinv_i(PRIMORAL, p);
    //printf("%d %d %d %d\n", i, p, PRIMORAL, inv);
    for(int j=0;j<6;j++){
      mpz_add_ui(candidate_plus, rop, offsets[j]);
      unsigned long cmodp = mpz_fdiv_ui(candidate_plus, p);
      int index = ((p - cmodp)*inv)%p;
      sieve_loc[i*6+j] = index;
    }
  }
      //mpz_t tmp;
      //mpz_t candidate;
      //mpz_inits(tmp, candidate);
      //mpz_init(tmp);
      //mpz_init(candidate);

  while(1) {
    //step 1
    unsigned int sieve_loc_index = 0;
    unsigned int i, j, o;
    //printf("Sieving\n");

    //#pragma omp parallel for
    for(i = 0; i < prime_size; i++) {
      for(j = 0; j < 6; j++){
	o = sieve_loc[sieve_loc_index];
	while(o < SIEVESIZE){
	  sieve[o] = false;
	  o += primes[i];
	}
	sieve_loc[sieve_loc_index] = o - SIEVESIZE;
	sieve_loc_index += 1;	
      } 
    }

    //printf("Testing Candidates\n");
      
  //    mpz_t tmp;
  //    mpz_t candidate;
      //mpz_inits(tmp, candidate);
  //    mpz_init(tmp);
  //    mpz_init(candidate);

    #pragma omp parallel
    { 
      mpz_t tmp;
      mpz_t candidate;
      //mpz_inits(tmp, candidate);
      mpz_init(tmp);
      mpz_init(candidate);
    #pragma omp for
    for(i = 0; i < SIEVESIZE; i++) {
      mpz_set_ui(tmp, i);
      mpz_addmul_ui(tmp, sievesize, loop_count); 
      mpz_set(candidate, rop);
      mpz_addmul_ui(candidate, tmp, PRIMORAL);

      if(sieve[i]){
	if(is_fermat_valid(candidate)){
	  mpz_out_str(fp, 10, candidate);
	  exit(0);
	}
      }
      sieve[i] = true;
    }
      mpz_clear(tmp);
      mpz_clear(candidate);
    }

    loop_count++;
  }

  return 0;
}
