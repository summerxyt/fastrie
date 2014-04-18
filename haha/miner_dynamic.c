#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <pthread.h>
#include "fermat.h"

// (2*3*5*7)
#define PRIMORAL 210

pthread_t finish_id = -1; 
mpz_t num;
pthread_mutex_t mutex;
pthread_mutex_t wheel;

void rieminer(void *arg) 
{
    int i = *(int* )arg;
    char start[100] = "0x801A2F60588BF10BB614D6796A726025F88C7156E3FBDF68685FC0617F4266358C0000000000";
    mpz_t min, max;
    mpz_t rop;
    mpz_t tmp;
    mpz_t cur;
    unsigned long r;
    unsigned int rv=0;
    //FILE* fp = stdout;
    //int i=0;

    mpz_init(min);
    mpz_init(cur);
    mpz_init(rop);
    mpz_init(tmp);

    mpz_set_str(min, start, 0);
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
        pthread_mutex_lock(&wheel);
        mpz_set(tmp,num);
        mpz_add_ui(num, num, 4*PRIMORAL);
        pthread_mutex_unlock(&wheel);
        for(int i = 0;i < 4;i++){
            mpz_add(cur, rop, tmp);
            if(is_fermat_valid(cur)){
                pthread_mutex_lock(&mutex);
                finish_id = pthread_self();
                printf("thread %d find rie prime!\n", finish_id);
                pthread_mutex_unlock(&mutex);
                gmp_printf ("An RIE found 0x%Zx\n", cur);  
                break;
            }
            mpz_add_ui(tmp, tmp, PRIMORAL);
        }
    }

    //mpz_out_str(fp, 10, rop);
    // pthread_mutex_lock(&mutex);
    // finish_id = pthread_self();
    // printf("thread %d find rie prime!\n", finish_id);
    // pthread_mutex_unlock(&mutex);
}

int main(int argc, char *argv[])
{
    pthread_t id[4];
    mpz_init(num);
    //pthread_mutex_init(&mutex, NULL);

    for(int i = 0;i < 4;i++){    
        pthread_create(&id[i], NULL, (void *)rieminer, &i);
        pthread_detach(id[i]);
    }

    while(1){
        pthread_mutex_lock(&mutex);
        if(finish_id != -1){
            //pthread_join(finish_id, NULL); 
            break;
        }
        pthread_mutex_unlock(&mutex);
    }
   
    pthread_mutex_destroy(&mutex);
    exit(0);
}



