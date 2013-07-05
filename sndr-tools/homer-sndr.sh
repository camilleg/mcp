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

# Learn models
# echo ...
time $bin/sndr $pt -i $snd/simp.wav -g $snd/homer.wav -M models/homer

# Find Homer
time $bin/sndr $pt -i $snd/simp.wav -g $snd/homer.wav -M models/homer

# Find not-Homer
time $bin/sndr $pc -i $snd/simp.wav -m models/homer-target models/homer-ubm -d homer
