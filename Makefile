CFLAGS = -O2 -Wall
OBJS = fft.o wavetrans.o wave.o

.PHONY: all clean

all: wavetrans

wavetrans: $(OBJS)

clean:
	rm -f *~ *.o wavetrans
