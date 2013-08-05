#/bin/bash

# Binaries
bin=bin

# Sounds
snd=../simpsons

# Parameters
set pt = "-K 8"
set pc = "-p .9999"

echo "Learn and find Homer."
$bin/snd_track -i $snd/simp.wav -t $snd/homer.wav -f $snd/homer-track.wav $pt $pc
