CFLAGS := -Wall -I.. -I/usr/local/include # or -I/usr/local/src/ffmpeg
CFLAGS += -D__CHECK
CFLAGS += -g -O0
#CFLAGS += -O3

# These libs were found by typing "make -n" in e.g /usr/local/src/ffmpeg/doc/examples .
LIBS_FFMPEG = -lavdevice -lpostproc -lavformat -lavcodec  -lswresample -lswscale -lavutil
LIBS_X = -lXfixes -lXext -lX11
LIBS_CODEC = -lSDL -lvorbisenc -lvorbis -ltheoraenc -ltheoradec -logg -lspeex -lrtmp -lgnutls -lopencore-amrwb -lopencore-amrnb -lmp3lame -lfdk-aac
LIBS_COMMON = -lfftw3 -lm -dl -lbz2 -lz -lrt

LIBS_CYGWIN = $(LIBS_X)                  $(LIBS_FFMPEG) $(LIBS_CODEC) -liconv                                        $(LIBS_COMMON)
LIBS_UBUNTU = $(LIBS_X) -L/usr/local/lib $(LIBS_FFMPEG) $(LIBS_CODEC) -lx264 -lvpx -lfaac -lass -lva -ljack -lasound $(LIBS_COMMON)

all: aufeat filt

clean: 
	rm -f aufeat filt *.o

test: aufeat
	./aufeat -i 80s.wav
	./aufeat -i /r/mp3/marvin/esoteria/audio_18.mp3

aufeat : aufeat.o auffmpeg.o
	g++ -o $@ $^ $(LIBS_UBUNTU)

filt : read-ffmpeg.o
	g++ -o $@ $^ $(LIBS_UBUNTU)

aufeat.o : aufeat.cpp ../array.h aquat.h aufeat.h optparse.h
	g++ $(CFLAGS) -c $<

auffmpeg.o : auffmpeg.cpp ../array.h
	g++ $(CFLAGS) -c $<

filt.o : read-ffmpeg.cpp
	g++ $(CFLAGS) -c $<
