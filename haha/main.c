#include "miner.h"
#include "fermat.h"

//2*3*5*7
#define PRIMORAL (2*3*5*7)

int main(int argc, char **argv){
  mpz_t min, max;
  mpz_t rop, tmp;
  unsigned long r;
  unsigned int rv=0;
  FILE* fp = stdout;
  int i=0;

  mpz_init(min);
  //mpz_init(max);
  mpz_init(rop);
  mpz_init(tmp);

  mpz_set_str(min, argv[1], 0);
  //mpz_set_str(max, argv[2], 0);

  r = mpz_fdiv_ui(min, PRIMORAL);
  r = (97 - r + PRIMORAL) % PRIMORAL;

  mpz_add_ui(rop, min, r);
  /*while(mpz_cmp(rop, max) < 0) {
    //mpz_add_ui(rop, rop, 210);
    if(is_fermat_valid(rop))
      rv++;
    mpz_add_ui(tmp, rop, PRIMORAL);
    mpz_set(rop, tmp);
  }*/

  while(1) {
    if(is_fermat_valid(rop)){
      break;
    }
    mpz_add_ui(tmp, rop, PRIMORAL);
    mpz_set(rop, tmp);
    //i++;
    //printf("%d\n", i);
  }

  mpz_out_str(fp, 10, rop);
  return 0;
}
