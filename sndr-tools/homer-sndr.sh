#/bin/bash
set -e		 # Abort if any command fails.

# Binaries
bin=bin

# Sounds
snd=../simpsons

# Parameters
#set pt = "-K 4 -H 4 -T .2 -t .5"
#set pc = "-f 10 -w 2 -p .9999"
set pt = "-K 8"
set pc = "-p .9999"

echo "Learn model."
$bin/sndr $pt -i $snd/simp.wav -g $snd/homer.wav -M models/homer

echo
echo "Find Homer and not-Homer (into simp.wav.homer.?.wav)."
$bin/sndr $pc -i $snd/simp.wav -m models/homer-target models/homer-ubm -M models/homer-combined -d homer
