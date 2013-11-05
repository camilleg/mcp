#!/usr/bin/env ruby

# Two classes of sounds: 440 and 1000 Hz pure sinusoids.
# Splice them randomly, in a ratio from 5 to 95%.
# Classify the spliced recording.
# The test passes iff the ratio of the number of frames is about the same as the splicing ratio.

# I want a sequence of intervals a,b,a,b,a,b,...
# where each a and b has a distribution of durations, uniform for convenience,
# and where the ratio of sums of a- and b-durations is specified.

# TODO: generalize to n classes.
# Generalize fraction to a random set of weights that sum to 1:
# start with n weights, from [1/n,...,1/n] to [really skewed, even with some zeros],
# and then permute it randomly.

$fraction = 0.44 # of $Freqs[1]
$ratio = $fraction / (1.0 - $fraction)		# 0.001 -> 0.001, 0.5 -> 1.0, 0.999 -> 999.0.
$Durs = $ratio<1.0 ? [1.0/$ratio, 1.0] : [1.0, $ratio]
$Classnames = ["A440", "1kHz"]
$Freqs = [440, 1000]
$histogram = [0,0]

$bin = "bin" # ../bin if run from sndr-tools/tests/
$snd = "/run/shm" # ramdisk
$models = "#$snd/models"
`rm -rf #$models; mkdir -p #$models`

puts "Generate test."

$Infiles = ["c0.wav", "c1.wav"]
$Testfile = "testcase.wav"
$SecMaxSplice = 30.0
$splices = []
`cd #$snd; rm -f #$Testfile splice*.wav #{$Infiles.join ' '}`
2.times {|i| `sox -n #$snd/#{$Infiles[i]} synth #$SecMaxSplice sin #{$Freqs[i]} gain -3`}
# "gain -3" dB leaves headroom.

def makeSplice(iOut,iIn)
  wav = "splice#{iOut}.wav"
  $splices << wav
  # Compute dur as a function of iIn.
  durMean = 0.2 * $Durs[iIn]
  if (durMean*1.5 > $SecMaxSplice)
    puts "Failed to generate test.  Retry, after growing $SecMaxSplice from #$SecMaxSplice to at least #{durMean*1.5}."
    exit -1
  end
  dur = rand(durMean*0.5 ... durMean*1.5)
  dur -= dur.modulo(1.0 / $Freqs[iIn])		# Round to nearest full period, to reduce clicks.
  $histogram[iIn] += dur
  `cd #$snd; sox #{$Infiles[iIn]} #{wav} trim 0 #{dur}`
end

# After generalizing to n, todo: instead of alternating 0 and 1, randomly permute [0,1,2,3,4, 0,1,2,3,4, ... ];
# and then tweak that permutation to avoid repeats (to avoid clicks)
# if I haven't yet rounded dur to a multiple of sine period.
20.times {|i| makeSplice i, i%2 }

puts $histogram
puts "Fraction: target #$fraction, actual #{$histogram[1] / ($histogram[0] + $histogram[1])}"

`cd #$snd; sox #{$splices.join ' '} #$Testfile`
#`audacity #$snd/#$Testfile`
# Sox howto: http://billposer.org/Linguistics/Computation/SoxTutorial.html

puts "Learn model for #{$Classnames[0]}."
$train = "-t .1 -H 2 -l 200 -h 2000 -K 2 -e 15"
$classify="-p .9999 -f 2"
2.times {|i|
  puts "#$bin/sndr #$train -i #$snd/#$Testfile -g #$snd/#{$Infiles[i]} -M #$snd/models/#{$Infiles[i]}"
  puts `#$bin/sndr #$train -i #$snd/#$Testfile -g #$snd/#{$Infiles[i]} -M #$snd/models/#{$Infiles[i]}`
  puts
}
