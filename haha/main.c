#include "miner.h"
#include "fermat.h"
#include <stdlib.h>
#include <gmp.h>

//2*3*5*7*11
#define PRIMORAL (2*3*5*7*11)

int main(int argc, char **argv)
{
    mpz_t min, max;
    unsigned int rv = 0;
    FILE *fp = stdout;
    char start[100] = "0x801A2F60588BF10BB614D6796A726025F88C7156E3FBDF68685FC0617F4266358C0000000000";

    mpz_init(min);
    mpz_init(max);

    mpz_set_str(min, start, 0);
    //mpz_set_str(max, argv[2], 0);

#define NUM_THREADS 4
    #pragma omp parallel num_threads(NUM_THREADS)
    {
	mpz_t rop, tmp;
	mpz_t rop1, rop2, rop3, rop4, rop5;
	mpz_init(rop);
	mpz_init(tmp);
	mpz_init(rop1);
	mpz_init(rop2);
	mpz_init(rop3);
	mpz_init(rop4);
	mpz_init(rop5);

	//int i = omp_get_thread_num();
	int i = omp_get_thread_num();
	mpz_add_ui(rop, min, i * PRIMORAL);
	unsigned long r = mpz_fdiv_ui(rop, PRIMORAL);
	r = (97 - r + PRIMORAL) % PRIMORAL;
	mpz_add_ui(rop1, rop, r);
	r = (937 - r + PRIMORAL) % PRIMORAL;
	mpz_add_ui(rop2, rop, r);
	r = (1147 - r + PRIMORAL) % PRIMORAL;
	mpz_add_ui(rop3, rop, r);
	r = (1357 - r + PRIMORAL) % PRIMORAL;
	mpz_add_ui(rop4, rop, r);
	r = (2197 - r + PRIMORAL) % PRIMORAL;
	mpz_add_ui(rop5, rop, r);


	while (1) {
	    if (is_fermat_valid(rop1)) {
		mpz_set(rop, rop1);
		break;
	    }
	    if (is_fermat_valid(rop2)) {
		mpz_set(rop, rop2);
		break;
	    }
	    if (is_fermat_valid(rop3)) {
		mpz_set(rop, rop3);
		break;
	    }
	    if (is_fermat_valid(rop4)) {
		mpz_set(rop, rop4);
		break;
	    }
	    if (is_fermat_valid(rop5)) {
		mpz_set(rop, rop5);
		break;
	    }
	    mpz_add_ui(tmp, rop1, NUM_THREADS * PRIMORAL);
	    mpz_set(rop1, tmp);
	    mpz_add_ui(tmp, rop2, NUM_THREADS * PRIMORAL);
	    mpz_set(rop2, tmp);
	    mpz_add_ui(tmp, rop3, NUM_THREADS * PRIMORAL);
	    mpz_set(rop3, tmp);
	    mpz_add_ui(tmp, rop4, NUM_THREADS * PRIMORAL);
	    mpz_set(rop4, tmp);
	    mpz_add_ui(tmp, rop5, NUM_THREADS * PRIMORAL);
	    mpz_set(rop5, tmp);
	}
	mpz_out_str(fp, 10, rop);
	exit(0);
    }

    //mpz_add_ui(rop, min, r);
    /*while(mpz_cmp(rop, max) < 0) {
       //mpz_add_ui(rop, rop, 210);
       if(is_fermat_valid(rop))
       rv++;
       mpz_add_ui(tmp, rop, PRIMORAL);
       mpz_set(rop, tmp);
       } */


    return 0;
}
