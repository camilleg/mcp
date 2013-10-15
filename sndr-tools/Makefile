# Compiler
CXX= g++ # icpc

# Flags.  Don't define -ffast-math, because that breaks isnan().
CXXFLAGS = -O3 -I..
CXXFLAGS += -D__NO_RAND_SEED
ifeq ($(CXX),icpc)
	CXXFLAGS += -xHOST -ip -no-prec-div -mdynamic-no-pic -openmp -openmp-report0 -vec-report0 -fno-alias
else
	CXXFLAGS += -march=native -fopenmp -pedantic -Wall -D_GNU_SOURCE
	# _GNU_SOURCE defines NaN
endif
#CXXFLAGS += -g -O0
#CXXFLAGS += -D__CHECK
LDFLAGS = -lfftw3

# Use float precision (not recommended)
# CXXFLAGS += -D__FLOAT

# Make HMM speaker models
CXXFLAGS += -D__HMM_TRAIN

# Directories
BIN = bin
SRC = src

# Rules

all: $(BIN)/snd_track $(BIN)/sndr

$(BIN):
	mkdir $(BIN)

$(BIN)/snd_track: $(SRC)/snd_track.cpp $(SRC)/../../aufeat.h $(SRC)/logadd.h $(SRC)/gmm.h $(SRC)/hmm.h $(SRC)/../../array.h $(SRC)/intransp.h $(SRC)/../../optparse.h $(SRC)/../../wavfile.h Makefile
	$(CXX) $(SRC)/snd_track.cpp $(CXXFLAGS) $(LDFLAGS) -o $(BIN)/snd_track

$(BIN)/sndr: $(SRC)/sndr.cpp $(SRC)/sndr.h $(SRC)/../../aufeat.h $(SRC)/logadd.h $(SRC)/gmm.h $(SRC)/hmm.h $(SRC)/../../array.h $(SRC)/intransp.h $(SRC)/../../optparse.h $(SRC)/../../wavfile.h Makefile
	$(CXX) $(SRC)/sndr.cpp $(CXXFLAGS) $(LDFLAGS) -o $(BIN)/sndr

clean:
	-rm -rf $(BIN)/*
