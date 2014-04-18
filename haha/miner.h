#include<stdio.h>
#include<stdlib.h>
#include<gmp.h>

int egcd_i(int a, int b, int *rx, int *ry) {
  int x, y, u, v, m, n;
  x = v = 0;
  y = u = 1;
  div_t d;

  while (a) {
    d = div(b, a);
    m = x - u * d.quot;
    n = y - v * d.quot;
    b = a;
    a = d.rem;
    x = u;
    y = v;
    u = m;
    v = n;
  }

  *rx = x;
  *ry = y;
 
  return b; 
}

int modinv_i(int a, int m) {
  int gcd, x, y;

  gcd = egcd_i(a, m, &x, &y);
  if (gcd != 1)
    return 0;
  else
    return (x+m) % m;
}
