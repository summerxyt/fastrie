sieve:sieve.c
	gcc -O3 -fopenmp sieve.c -o sieve
fermat:fermat_test.c
	gcc -O3 fermat_test.c -o fermat -lgmp
test2:fermat.h main.c
	g++ -O3 -fopenmp main.c -o test2 -lgmp -lprimesieve
test:fermat.h test.c
	g++ -O3 -fopenmp test.c -o test -lgmp -lprimesieve
miner:fermat.h miner_static.c
	gcc -O3 miner_static.c -o miner -lgmp
