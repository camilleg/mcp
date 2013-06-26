all: streamtest
	./streamtest 127.0.0.1 49990 localhost.sdp

streamtest: streamtest.cpp Makefile
	g++ -Wall -g streamtest.cpp -o streamtest \
	  -L/usr/local/lib -lavdevice -lpostproc -lavformat -lavcodec -ldl -lva -lXfixes -lXext -lX11 -ljack -lasound -lSDL -lx264 -lvpx -lvorbisenc -lvorbis -ltheoraenc -ltheoradec -logg -lspeex -lrtmp -lgnutls -lopencore-amrwb -lopencore-amrnb -lmp3lame -lfdk-aac -lfaac -lass -lbz2 -lz -lrt -lswresample -lswscale -lavutil -lm
