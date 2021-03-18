all: simple gmp

gyt-common.o:
	g++ -O4 -c -o gyt-common.o gyt-common.cpp

simple: gyt-common.o
	g++ -O4 -o ../gyt-2d     gyt-2d.cpp     gyt-common.o
	g++ -O4 -o ../gyt-2d-all gyt-2d-all.cpp gyt-common.o
	g++ -O4 -o ../gyt        gyt.cpp        gyt-common.o
	g++ -O4 -o ../gyt-pq     gyt-pq.cpp     gyt-common.o
	g++ -O4 -o ../gyt-rand   gyt-rand.cpp   gyt-common.o
	g++ -O4 -o ../gyt-proba  gyt-proba.cpp  gyt-common.o

gyt-common-gmp.o:
	g++ -O4 -c -o gyt-common-gmp.o gyt-common-gmp.cpp

gmp: gyt-common-gmp.o
	g++ -O4 -o ../gyt-2d-gmp      gyt-2d-gmp.cpp      gyt-common-gmp.o -lgmpxx -lgmp
	g++ -O4 -o ../gyt-2d-all-gmp  gyt-2d-all-gmp.cpp  gyt-common-gmp.o -lgmpxx -lgmp
	g++ -O4 -o ../gyt-gmp         gyt-gmp.cpp         gyt-common-gmp.o -lgmpxx -lgmp
	g++ -O4 -o ../gyt-pq-gmp      gyt-pq-gmp.cpp      gyt-common-gmp.o -lgmpxx -lgmp
	g++ -O4 -o ../gyt-rand-gmp    gyt-rand-gmp.cpp    gyt-common-gmp.o -lgmpxx -lgmp
	g++ -O4 -o ../gyt-proba-gmp   gyt-proba-gmp.cpp   gyt-common-gmp.o -lgmpxx -lgmp

.PHONY: clean scratch

clean:
	rm -f *.o

scratch: clean
	rm -f *~
