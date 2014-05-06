#include "mpi.h"

#include <gmp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include "fermat.h"
#include "miner.h"

#define  MASTER  0
//2*3*5*7
#define PRIMORAL (2*3*5*7)
#define SIEVESIZE (1024*512)
#define MAX_SIEVE_PRIME (1000000)

int* primesieve_generate_primes(int prime_size, int *count){

    int *re = (int* )malloc(sizeof(int) * (prime_size+1)/2);
    bool *isprime = (bool* )malloc(sizeof(bool) * (prime_size+1));
    
    memset(isprime, true, prime_size+1);

    int sN = int(floor(sqrt(prime_size)));
    
    *count = 1;
    re[0] = 2;

    for(int i = 3;i < prime_size+1;i+=2){ 
        if (isprime[i]){
            re[*count] = i;
        *count = *count + 1;
            if (i < sN){
                int ni = 2*i;
                while (ni <= prime_size){
                    isprime[ni] = false;
                    ni += i;
                }
            }
        }
    }

    return re;
}


int main (int argc, char *argv[])
{
    int len;
    char hostname[MPI_MAX_PROCESSOR_NAME];
    
    mpz_t min;
    mpz_t max;
    mpz_t rop;
    mpz_t sievesize;

    char start[100] = "0x801A2F60588BF10BB614D6796A726025F88C7156E3FBDF68685FC0617F4266358C0000000000";
    unsigned long r;
    unsigned long loop_count;
    int  rank;
    int  np;        /* number of processes in job */
    bool *sieve;
    int *primes;
    int *primes_start;
    int *sieve_loc;
    int prime_size;
    int loc_size;
    struct timeval t_start;
    struct timeval t_end;

    mpz_init(min);
    mpz_init(max);
    mpz_init(rop);
    mpz_init(sievesize);

    // step 1 & 2: round up to nearest primorial
    mpz_set_str(min, start, 0);
    mpz_set_ui(sievesize, SIEVESIZE);

    r = mpz_fdiv_ui(min, PRIMORAL);
    r = (97 - r + PRIMORAL) % PRIMORAL;

    mpz_add_ui(rop, min, r);

    printf("init done\n");

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Get_processor_name(hostname, &len);

    sieve = (bool* )malloc(sizeof(bool) * SIEVESIZE * 4);

    printf ("Hello from task %d on %s!\n", rank, hostname);

    // prime.append(p)
    primes_start = primesieve_generate_primes(MAX_SIEVE_PRIME, &prime_size);

    // pass the first 7 primes
    for(primes = primes_start;*primes <= 7;primes++);

    prime_size -= (primes - primes_start);
    loc_size = 6 * prime_size;
    sieve_loc = (int *)malloc(sizeof(int) * loc_size);

    mpz_t candidate_plus;
    mpz_init(candidate_plus);

    gettimeofday(&t_start, NULL);

    for(int i = 0;i < prime_size;i++){
        // if(primes[i] <= 7)
        //   continue;
        int p = primes[i];
        int inv = modinv_i(PRIMORAL, p);

        for(int j = 0;j < 6;j++){
            mpz_add_ui(candidate_plus, rop, offsets[j]);
            unsigned long cmodp = mpz_fdiv_ui(candidate_plus, p);
            int index = ((p - cmodp) * inv) % p;
            // sieve_loc.append(index)
            sieve_loc[i*6+j] = index;
        }
    }

    // wait for all processes to generate sievesize
    MPI_Barrier(MPI_COMM_WORLD);
    
    gettimeofday(&t_end, NULL);

    double duration = t_end.tv_sec - t_start.tv_sec + (t_end.tv_usec - t_start.tv_usec) / 1000000.0;
    printf("Process %d use Time in creating prime index: %f\n", rank, duration);
    

    bool isStop[np];
    loop_count = 0;

    while(1) {
        //step 1
        memset(isStop, false, np);
        memset(sieve, true, SIEVESIZE);
        unsigned int sieve_loc_index = 0;
        unsigned int o;
        //printf("Sieving\n");

        gettimeofday(&t_start, NULL);

        // each process test its own set of prime
        int chunk = (prime_size + np - 1) / np;
        int index_start = rank * chunk;
        int index_end = index_start + chunk;

        if(index_end > prime_size){
            index_end = prime_size;
        }

        sieve_loc_index = 6 * index_start;

        for(int i = index_start; i < index_end; i++) {
            for(int j = 0; j < 6; j++){
                o = sieve_loc[sieve_loc_index];
                while(o < SIEVESIZE){
                    sieve[o] = false;
                    o += primes[i];
                }

                // prepare for the next round
                sieve_loc[sieve_loc_index] = o - SIEVESIZE;
                sieve_loc_index += 1;   
            } 
        }

        gettimeofday(&t_end, NULL);

        duration = t_end.tv_sec - t_start.tv_sec + (t_end.tv_usec - t_start.tv_usec) / 1000000.0;
        printf("Process %d use Time in creating sieve: %f\n", rank, duration);

        // wait each process finish its part of sieve
        MPI_Barrier(MPI_COMM_WORLD);

        gettimeofday(&t_start, NULL);

        MPI_Allgather(sieve, SIEVESIZE, MPI_C_BOOL, sieve, SIEVESIZE, MPI_C_BOOL, MPI_COMM_WORLD);
        //printf("Testing Candidates\n");
        for(int i = 0;i < SIEVESIZE;i++){
            for(int j = 0;j < np;j++){
                sieve[i] = sieve[i] && sieve[i+j*SIEVESIZE];
                if(!sieve[i]){
                    break;
                } 
            }
        }

        gettimeofday(&t_end, NULL);

        duration = t_end.tv_sec - t_start.tv_sec + (t_end.tv_usec - t_start.tv_usec) / 1000000.0;
        printf("Process %d use Time in combining sieve: %f\n", rank, duration);

        mpz_t tmp;
        mpz_t candidate;
        //mpz_inits(tmp, candidate);
        mpz_init(tmp);
        mpz_init(candidate);
        
        gettimeofday(&t_start, NULL);

        for(int i = index_start; i < index_end; i++) {

            if(sieve[i]){
                mpz_set_ui(tmp, i);
                mpz_addmul_ui(tmp, sievesize, loop_count); 
                mpz_set(candidate, rop);
                mpz_addmul_ui(candidate, tmp, PRIMORAL);

                if(is_fermat_valid(candidate)){
                    gmp_printf("An RIE found 0x%Zx\n", candidate);
                    isStop[0] = true;
                }
            }
            
            //sieve[i] = true;
        }

        gettimeofday(&t_end, NULL);

        duration = t_end.tv_sec - t_start.tv_sec + (t_end.tv_usec - t_start.tv_usec) / 1000000.0;
        printf("Process %d use Time in testing candidate: %f\n", rank, duration);

        MPI_Barrier(MPI_COMM_WORLD);            


        if(rank == MASTER){
            gmp_printf("MASTER: Number of test now is: 0x%Zx\n",candidate);
        }

        MPI_Allgather(isStop, 1, MPI_C_BOOL, isStop, 1, MPI_C_BOOL, MPI_COMM_WORLD);

        // decide if need to exit
        for(int j = 1;j < np;j++){
            isStop[0] = isStop[0] || isStop[j];
            if(isStop[0]){
                break;
            } 
        }

        if(isStop[0]){
            break;
        }

        loop_count++;
    }





    if(rank == MASTER){
            printf("MASTER: Number of MPI tasks is: %d\n",np);
    }
    
    MPI_Finalize();

    return 0;
}
