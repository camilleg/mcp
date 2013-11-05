CXX= g++ # icpc

CXXFLAGS = -O3 -I.. -std=c++0x
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
# Omit -ffast-math, because it breaks isnan().

# Use float precision (not recommended)
# CXXFLAGS += -D__FLOAT

# Make HMM not GMM speaker models
CXXFLAGS += -D__HMM_TRAIN


TST = tests
BIN = bin
SRC = src
SRC_H = $(SRC)/../..
H = $(SRC)/ann.h $(SRC)/gmm.h $(SRC)/hmm.h $(SRC)/intransp.h $(SRC)/logadd.h $(SRC)/real_t.h $(SRC)/state.h \
    $(SRC_H)/array.h $(SRC_H)/aufeat.h $(SRC_H)/optparse.h $(SRC_H)/wavfile.h

test: $(BIN)/sndr $(TST)/sines.rb
	 $(TST)/sines.rb
all: $(BIN)/snd_track $(BIN)/sndr
$(BIN):
	mkdir $(BIN)
$(BIN)/snd_track: $(BIN) $(SRC)/snd_track.cpp $(H) Makefile
	$(CXX) $(SRC)/snd_track.cpp $(CXXFLAGS) $(LDFLAGS) -o $(BIN)/snd_track
$(BIN)/sndr: $(BIN) $(SRC)/sndr.cpp $(SRC)/sndr.h $(H) Makefile
	$(CXX) $(SRC)/sndr.cpp $(CXXFLAGS) $(LDFLAGS) -o $(BIN)/sndr
clean:
	-rm -rf $(BIN)/*
