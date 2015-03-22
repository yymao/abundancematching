CC=gcc
CFLAGS=-D_BSD_SOURCE -D_POSIX_SOURCE -D_POSIX_C_SOURCE=200809L -D_SVID_SOURCE -D_DARWIN_C_SOURCE -Wall -fno-math-errno -fPIC
LDFLAGS=-shared
OFLAGS=-lm -O3 -std=c99

OBJS=fiducial_deconvolute.so

all: $(OBJS)

fiducial_deconvolute.so: fiducial_deconvolute.c
	$(CC) $(CFLAGS) $(LDFLAGS) $< -o $@ $(OFLAGS)

clean:
	rm -f $(OBJS)

.PHONY: clean
