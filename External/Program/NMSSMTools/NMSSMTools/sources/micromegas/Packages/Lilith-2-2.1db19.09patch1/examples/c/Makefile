# 'make'        build executable file 'lilith_compute'
# 'make clean'  removes all .o and executable files

# the path to the Lilith C/C++ API and to the main file should be given
LILITHCAPI = ../../lilith/c-api
MAIN = lilith_compute

# all lines below do not need to be modified
CC = gcc
CFLAGS += $(shell python2-config --cflags)
INCLUDES = -I$(LILITHCAPI)
LFLAGS += $(shell python2-config --ldflags)
SRCS = $(MAIN).c $(LILITHCAPI)/lilith.c
OBJS = $(SRCS:.c=.o)

.PHONY: clean

all: $(MAIN)

$(MAIN): $(OBJS) 
	$(CC) $(CFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS) $(LFLAGS)

.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $<  -o $@

clean:
	$(RM) *.o *~ $(MAIN) $(LILITHCAPI)/*.o

