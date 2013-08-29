#/bin/bash
set -e		 # Abort if any command fails.

# Binaries
bin=bin

# Sounds
snd=../simpsons

# Parameters
train="-t .25 -H 8 -T .4 -l 45 -h 4500 -K 7 -e 120"
classify="-p .9999 -f 2"

echo "Learn model."
$bin/sndr $train -i $snd/simp.wav -g $snd/homer.wav -M models/homer

echo
echo "Find Homer and not-Homer (into simp.wav.homer.*.wav)."
$bin/sndr $train $classify -i $snd/simp.wav -m models/homer-target models/homer-ubm -M models/homer-combined -d homer -D homer
