#/bin/bash

# Binaries
bin=bin

# Sounds
snd=../simpsons

# Parameters
#set pt = "-K 4 -H 4 -T .2 -t .5"
#set pc = "-f 10 -w 2 -p .9999"
set pt = "-K 8"
set pc = "-p .9999"

# Learn and find Homer
time $bin/snd_track -i $snd/simp.wav -t $snd/homer.wav -f $snd/homer-track.wav $pt $pc
