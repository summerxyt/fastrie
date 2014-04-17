#include<stdio.h>
//#include<gmp.h>
//#include<omp.h>

#define PRIMORIAL_NUM 5
#define PRIMORIAL (2*3*5*7*11)
unsigned int primes[PRIMORIAL_NUM] = {2, 3, 5, 7, 11};
int offset[6] = {0, 4, 6, 10, 12, 16};

unsigned char is_candidate[PRIMORIAL] = {0};
//int is_candidate[PRIMORIAL] = {0};

void gen_candidate(int primorial) {
  int i,j,ni;

  // #pragma omp parallel for
  for(i = 0; i < PRIMORIAL_NUM; i++){
    ni = primes[i];
    while(ni < primorial + 16) {
      for(j = 0; j < 6; j++){
	if((ni - offset[j] > 0) && (ni - offset[j] < primorial)) {
	  is_candidate[ni-offset[j]] = 0x01;
	}
      }
      ni += primes[i];
    }
  }

  for(i = 3; i < primorial; i++) {
    if(!is_candidate[i])
      printf("%d\t", i);
  }
  printf("\n");
}

int main() {
  gen_candidate(PRIMORIAL);
  return 0;  
}
