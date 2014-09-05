all: crlprop

crlprop: crlprop.c copyparam.c functions.c crlprop.h
	gcc -o $@ crlprop.c copyparam.c functions.c -Wall -Wextra -Werror -O2 -lm -lfftw3

