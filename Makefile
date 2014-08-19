all: crlprop

crlprop: crlprop.c crlprop.h
	gcc -o $@ $< -Wall -Wextra -Werror -O2 -lm -lfftw3

