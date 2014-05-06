all:main why mains
main:fermat.h advance.c
	g++ -O3 -fopenmp advance.c -o main -lgmp -lprimesieve
mains:fermat.h adwhys.c
	g++ -O3 adwhys.c -fopenmp -o mains -lgmp -lprimesieve
why:fermat.h adwhy.c
	g++ -O3 -fopenmp adwhy.c -o why -lgmp -lprimesieve
clean:
	rm main why
