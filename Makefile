# CFLAGS = -O2 -Wall
CFLAGS = -g -Wall
LDFLAGS = -lm -lefence
OBJS = fft.o wavetrans.o wave.o

.PHONY: all clean

all: wavetrans

wavetrans: $(OBJS)

test: fft.o test.o

clean:
	rm -f *~ *.o wavetrans
