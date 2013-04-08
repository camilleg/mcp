test: aufeat
	./aufeat

aufeat: aufeat.cpp aquat.h aufeat.h optparse.h wavfile.h
	g++ -g -O0 -Wall -D__CHECK -I.. aufeat.cpp -lfftw3 -o $@
