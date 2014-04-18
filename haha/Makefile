sieve:sieve.c
	gcc -O3 -fopenmp sieve.c -o sieve
fermat:fermat_test.c
	gcc -O3 fermat_test.c -o fermat -lgmp
test:fermat.h main.c
	gcc -O3 main.c -o test -lgmp
miner:fermat.h miner_static.c
	gcc -O3 miner_static.c -o miner -lgmp