CC=gcc
CFLAGS=-Wall -O3 -std=c99 -fPIC -I ~/libs/include -L ~/libs/lib
LDFLAGS=-lmol2 -ljansson -lm

all: minimize minimize_list

minimize: minimize.c
	$(CC) $(CFLAGS) $^ $(LDFLAGS) -o minimize

minimize_list: minimize_list.c
	$(CC) $(CFLAGS) $^ $(LDFLAGS) -o minimize_list
