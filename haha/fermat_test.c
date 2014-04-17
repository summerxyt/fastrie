#include<stdio.h>
#include<gmp.h>

int fermat_testi(unsigned int n) {
  int i, re;
  unsigned int pow = 1;

  if(n == 2)
    return 1;
  if(!(n & 0x0001))
    return 0;

  for(i = 0; i < n - 1; i++)
    pow *= 2;

  return (pow % n == 1) ? 1 : 0;
}

int fermat_testz(mpz_t n) {
  int i, re;
  mpz_t two, rop, n_minus_one;

  if(mpz_cmp_ui(n, 2) == 0)
    return 1;

  if(mpz_even_p(n))
    return 0;

  mpz_init(two);
  mpz_set_ui(two, 2);
  mpz_init(rop);
  mpz_init(n_minus_one);
  mpz_sub_ui(n_minus_one, n, 1);
  mpz_powm(rop, two, n_minus_one, n);
  if(mpz_cmp_ui(rop, 1) == 0)
    return 1;

  return 0;
}

int main(){
  int i=2;
  int rv;
  mpz_t n;

  for(;i<20;i++){
    printf("%d:%d\n", i, fermat_testi(2)); 
  }

  mpz_init(n);
  for(i=2;i<20;i++){
    mpz_set_ui(n, i);
    rv = fermat_testz(n);
    printf("%d:%d\n", i, rv); 
  }

  return 0;
}
