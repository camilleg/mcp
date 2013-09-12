# Compiler
CXX= g++
#CXX = icpc

# Flags
CXXFLAGS = -g -O0 -I.. -pedantic -Wall -march=native
CXXFLAGS += -D__NO_RAND_SEED
ifeq ($(CXX),icpc)
	CXXFLAGS += -xHOST -O3 -ip -no-prec-div -mdynamic-no-pic -openmp -openmp-report0 -vec-report0 -fno-alias
else
	CXXFLAGS += -O3 -march=core2 -fopenmp
endif
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

$(BIN)/snd_track: $(SRC)/snd_track.cpp $(SRC)/../../aufeat.h $(SRC)/gmm.h $(SRC)/hmm.h $(SRC)/../../array.h $(SRC)/intransp.h $(SRC)/../../optparse.h $(SRC)/../../wavfile.h Makefile
	$(CXX) $(SRC)/snd_track.cpp $(CXXFLAGS) $(LDFLAGS) -o $(BIN)/snd_track

$(BIN)/sndr: $(SRC)/sndr.cpp $(SRC)/sndr.h $(SRC)/../../aufeat.h $(SRC)/gmm.h $(SRC)/hmm.h $(SRC)/../../array.h $(SRC)/intransp.h $(SRC)/../../optparse.h $(SRC)/../../wavfile.h Makefile
	$(CXX) $(SRC)/sndr.cpp $(CXXFLAGS) $(LDFLAGS) -o $(BIN)/sndr

clean:
	-rm -rf $(BIN)/*
