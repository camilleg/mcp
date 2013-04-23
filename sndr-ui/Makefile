test: aufeat
	./aufeat -i 80s.wav
	./aufeat -i /r/mp3/marvin/esoteria/audio_18.mp3

aufeat: aufeat.cpp auffmpeg.cpp aquat.h aufeat.h optparse.h wavfile.h ../array.h
	g++ -g -O0 -Wall -D__CHECK -I.. aufeat.cpp auffmpeg.cpp -lfftw3 -o $@ \
	  -pthread -L/usr/local/lib -lavdevice -lpostproc -lavformat -lavcodec -ldl -lva -lXfixes -lXext -lX11 -ljack -lasound -lSDL -lx264 -lvpx -lvorbisenc -lvorbis -ltheoraenc -ltheoradec -logg -lspeex -lrtmp -lgnutls -lopencore-amrwb -lopencore-amrnb -lmp3lame -lfdk-aac -lfaac -lass -lbz2 -lz -lrt -lswresample -lswscale -lavutil -lm

filt: read-ffmpeg.cpp Makefile
	g++ -Wall -O3 -I/usr/local/include read-ffmpeg.cpp -pthread -L/usr/local/lib -lavdevice -lpostproc -lavformat -lavcodec -ldl -lva -lXfixes -lXext -lX11 -ljack -lasound -lSDL -lx264 -lvpx -lvorbisenc -lvorbis -ltheoraenc -ltheoradec -logg -lspeex -lrtmp -lgnutls -lopencore-amrwb -lopencore-amrnb -lmp3lame -lfdk-aac -lfaac -lass -lbz2 -lz -lrt -lswresample -lswscale -lavutil -lm -s -o $@
# These long commands came from typing "make -n" in e.g /usr/local/src/ffmpeg/doc/examples .
# -I/usr/local/include or -I/usr/local/src/ffmpeg
