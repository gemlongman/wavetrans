# CFLAGS = -O2 -Wall
CFLAGS = -g -Wall
LDFLAGS = -lm -lefence
OBJS = fft.o wavetrans.o wave.o

.PHONY: all clean

all: wavetrans

wavetrans: $(OBJS)

fftexample: fft.o fftexample.o

test: all
	@echo "Making noise sample:"
	./wavetrans -s < noisy1.wav > noisy1.dat
	@echo "Removing noise:"
	./wavetrans -r noisy1.dat < noisy1.wav > noisy1-clean.wav
	@echo "Equalizing:"
	./wavetrans -e equalizer.txt < sample1.wav > sample1-equalized.wav
	@echo "Done!"
	@echo "Use playwave or your favourite player to play noisy1-clean.wav"
	@echo "and sample1-equalized.wav"

clean:
	rm -f *~ *.o wavetrans fftexample
