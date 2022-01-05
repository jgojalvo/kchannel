PROG = kch1d
OBJECTS = dran.o
MAINOBJ = $(PROG).o
MAINSRC = $(PROG).c

CC = cc -I$(HOME)/include
CFLAGS = -O3

all:	$(MAINOBJ) $(OBJECTS)
	$(CC) $(MAINOBJ) $(OBJECTS) -o $(PROG) -lm

$(MAINOBJ):	$(MAINSRC) $(OBJECTS)
	$(CC) $(CFLAGS) -c $(MAINSRC)

%.o : %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm *.o
