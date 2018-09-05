#/bin/bash

# gcc 
# gcc 
gcc wavetrans.c fft.c wave.c -o wavetrans -lm

# test:
@echo "Making noise sample:"
./wavetrans -s < mosquito03.wav > noisy1.dat
@echo "Removing noise:"
./wavetrans -r noisy1.dat < mosquito03.wav > mosquito03-clean.wav

@echo "Equalizing:"
./wavetrans -e equalizer.txt < mosquito03.wav > mosquito03-equalized.wav
@echo "Done!"
@echo "Use playwave or your favourite player to play noisy1-clean.wav"
@echo "and sample1-equalized.wav"

#test 2:
@echo "Making noise sample:"
./wavetrans -s < daleng5.wav > noisy1.dat
@echo "Removing noise:"
./wavetrans -r noisy1.dat < daleng5.wav > daleng5-clean.wav

@echo "Equalizing:"
./wavetrans -e equalizer.txt < daleng5.wav > daleng5-equalized.wav

# show:
../waveform/waveform mosquito03.wav --png mosquito03.png --png-color-outer 00ff00ff --png-width 16284 --png-height 256 
../waveform/waveform mosquito03-clean.wav --png mosquito03-clean.png --png-color-outer 00ff00ff --png-width 16284 --png-height 256 
../waveform/waveform mosquito03-equalized.wav --png mosquito03-equalized.png --png-color-outer 00ff00ff --png-width 16284 --png-height 256 


../waveform/waveform daleng5.wav --png daleng5.png --png-color-outer 00ff00ff --png-width 16284 --png-height 512 
../waveform/waveform daleng5-clean.wav --png daleng5-clean.png --png-color-outer 00ff00ff --png-width 16284 --png-height 512 
../waveform/waveform daleng5-equalized.wav --png daleng5-equalized.png --png-color-outer 00ff00ff --png-width 16284 --png-height 512 