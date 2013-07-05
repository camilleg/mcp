# compiler
CXX= g++
#CXX = icpc

# flags
CXXFLAGS = -I.. -pedantic -Wall -D__NO_RAND_SEED
ifeq ($(CXX),icpc)
	CXXFLAGS += -xHOST -O3 -ip -no-prec-div -mdynamic-no-pic -openmp -openmp-report0 -vec-report0 -fno-alias
else
	CXXFLAGS += -O3 -march=core2 -fopenmp
endif
LDFLAGS = -lfftw3

# Use float precision (not recommended)
# CXXFLAGS += -D__FLOAT

# Use AquaTerm to do some plots (Mac only)
# CXXFLAGS += -DAQUA
# LDFLAGS += -laquaterm

# Make HMM speaker models
#CXXFLAGS += -D__HMM_TRAIN

# directories
BIN = bin
SRC = src

# rules

all: $(BIN)/snd_track $(BIN)/sndr

$(BIN):
	mkdir $(BIN)

$(BIN)/snd_track: $(SRC)/snd_track.cpp $(SRC)/../../aufeat.h $(SRC)/gmm.h $(SRC)/hmm.h $(SRC)/../../array.h $(SRC)/intransp.h $(SRC)/aquat.h $(SRC)/optparse.h $(SRC)/../../wavfile.h
	$(CXX) $(SRC)/snd_track.cpp $(CXXFLAGS) $(LDFLAGS) -o $(BIN)/snd_track

$(BIN)/sndr: $(SRC)/sndr.cpp $(SRC)/sndr.h $(SRC)/../../aufeat.h $(SRC)/gmm.h $(SRC)/hmm.h $(SRC)/../../array.h $(SRC)/intransp.h $(SRC)/aquat.h $(SRC)/optparse.h $(SRC)/../../wavfile.h
	$(CXX) $(SRC)/sndr.cpp $(CXXFLAGS) $(LDFLAGS) -o $(BIN)/sndr

clean:
	rm -r -f $(BIN)/*
