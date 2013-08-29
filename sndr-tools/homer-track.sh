#/bin/bash

# Binaries
bin=bin

# Sounds
snd=../simpsons

# Parameters
train="-K 8"
classify="-p .9999"

echo "Learn and find Homer."
$bin/snd_track -i $snd/simp.wav -t $snd/homer.wav -f $snd/homer-track.wav $train $classify
