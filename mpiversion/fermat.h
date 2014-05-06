#ifndef FERMAT_H
#define FERMAT_H

#include<stdio.h>
#include<gmp.h>

static const unsigned int offsets[6] = {0, 4, 6, 10, 12, 16};

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
};

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
};

int is_fermat_valid(mpz_t n) {
  int i, rv;
  mpz_t temp;
  mpz_init(temp);
  mpz_set(temp, n);

  for(i = 0; i < 6; i++){
    mpz_add_ui(temp, n, offsets[i]);
    //serial
    if(fermat_testz(temp) == 0)
      return 0;
  }

  return 1;
}

#endif
