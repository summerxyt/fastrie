#include "miner.h"
#include "fermat.h"

//2*3*5*7
#define PRIMORAL (2*3*5*7*11)

int main(int argc, char **argv){
    mpz_t min, max;
    mpz_t rop, tmp;
    mpz_t rop1, rop2, rop3, rop4, rop5;
    unsigned long r;
    unsigned int rv=0;
    FILE* fp = stdout;
    int i=0;
    char start[100] = "0x801A2F60588BF10BB614D6796A726025F88C7156E3FBDF68685FC0617F4266358C0000000000";


    mpz_init(min);
    mpz_init(max);
    mpz_init(rop);
    mpz_init(tmp);

    mpz_init(rop1);
    mpz_init(rop2);
    mpz_init(rop3);
    mpz_init(rop4);
    mpz_init(rop5);

    mpz_set_str(min, start, 0);
    //mpz_set_str(max, argv[2], 0);

    r = mpz_fdiv_ui(min, PRIMORAL);
    r = (97 - r + PRIMORAL) % PRIMORAL;
    mpz_add_ui(rop1, min, r);
    r = (937 - r + PRIMORAL) % PRIMORAL;
    mpz_add_ui(rop2, min, r);
    r = (1147 - r + PRIMORAL) % PRIMORAL;
    mpz_add_ui(rop3, min, r);
    r = (1357 - r + PRIMORAL) % PRIMORAL;
    mpz_add_ui(rop4, min, r);
    r = (2197 - r + PRIMORAL) % PRIMORAL;
    mpz_add_ui(rop5, min, r);
    

    //printf("init done!\n");
    //mpz_add_ui(rop, min, r);
    /*while(mpz_cmp(rop, max) < 0) {
        //mpz_add_ui(rop, rop, 210);
        if(is_fermat_valid(rop))
            rv++;
        mpz_add_ui(tmp, rop, PRIMORAL);
        mpz_set(rop, tmp);
    }*/

    while(1) {
        if(is_fermat_valid(rop1)){
            mpz_set(rop, rop1);
            break;
        }
        if(is_fermat_valid(rop2)){
            mpz_set(rop, rop2);
            break;
        }
        if(is_fermat_valid(rop3)){
            mpz_set(rop, rop3);
            break;
        }
        if(is_fermat_valid(rop4)){
            mpz_set(rop, rop4);
            break;
        }
        if(is_fermat_valid(rop5)){
            mpz_set(rop, rop5);
            break;
        }

        mpz_add_ui(rop1, rop1, PRIMORAL);
        mpz_add_ui(rop2, rop2, PRIMORAL);
        mpz_add_ui(rop3, rop3, PRIMORAL);
        mpz_add_ui(rop4, rop4, PRIMORAL);
        mpz_add_ui(rop5, rop5, PRIMORAL);
        //i++;
        //printf("%d\n", i);
    }

    //printf("found!\n");
    mpz_out_str(fp, 10, rop);
    return 0;
}
