all:main why mains why2 why3
main:fermat.h advance.c
	g++ -O3 -fopenmp advance.c -o main -lgmp -lprimesieve
mains:fermat.h adwhys.c
	g++ -O3 adwhys.c -fopenmp -o mains -lgmp -lprimesieve
why:fermat.h adwhy.c
	g++ -O3 -fopenmp adwhy.c -o why -lgmp -lprimesieve
why2:fermat.h adwhy2.c
	g++ -O3 -fopenmp adwhy2.c -o why2 -lgmp -lprimesieve
why3:fermat.h adwhy3.c
	g++ -O3 -fopenmp adwhy3.c -o why3 -lgmp -lprimesieve
clean:
	rm main why why2 why3
