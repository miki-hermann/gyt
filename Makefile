all: simple gmp

simple:
	g++ -O4 -o young2d          young2d.cpp
	g++ -O4 -o young2d-all      young2d-all.cpp
	g++ -O4 -o gyt              gyt.cpp
	g++ -O4 -o gyt-pq           gyt-pq.cpp
	g++ -O4 -o gyt-rand         gyt-rand.cpp
	g++ -O4 -o gyt-proba        gyt-proba.cpp

gmp:
	g++ -O4 -o young2d-gmp      young2d-gmp.cpp      -lgmpxx -lgmp
	g++ -O4 -o young2d-all-gmp  young2d-all-gmp.cpp  -lgmpxx -lgmp
	g++ -O4 -o gyt-gmp          gyt-gmp.cpp          -lgmpxx -lgmp
	g++ -O4 -o gyt-pq-gmp       gyt-pq-gmp.cpp       -lgmpxx -lgmp
	g++ -O4 -o gyt-rand-gmp     gyt-rand-gmp.cpp     -lgmpxx -lgmp
	g++ -O4 -o gyt-proba-gmp    gyt-proba-gmp.cpp    -lgmpxx -lgmp

.PHONY: clean scratch

clean:
	rm -f young2d
	rm -f gyt gyt-pq gyt-rand gyt-proba
	rm -f *-gmp *-all

scratch: clean
	rm -f *~
