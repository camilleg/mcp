#/bin/bash
set -e		 # Abort if any command fails.

# Binaries
bin=bin

# Sounds
snd=../simpsons

# Parameters
train="-t .25 -H 8 -T .4 -l 45 -h 4500 -K 7 -e 120"
classify="-p .9999 -f 2"

echo "Learn model for Marge."
echo $bin/sndr $train -i $snd/simp.wav -g $snd/marge -M models/marge
$bin/sndr $train -i $snd/simp.wav -g $snd/marge.wav -M models/marge
echo
echo "Find Marge and not-Marge (into simp.wav.marge.*.wav)."
$bin/sndr $train $classify -i $snd/simp.wav -m models/marge-target models/marge-ubm -M models/marge-combined -d marge -D marge

echo "Learn model for Homer."
$bin/sndr $train -i $snd/simp.wav -g $snd/homer.wav -M models/homer
echo
echo "Find Homer and not-Homer (into simp.wav.homer.*.wav)."
$bin/sndr $train $classify -i $snd/simp.wav -m models/homer-target models/homer-ubm -M models/homer-combined -d homer -D homer
