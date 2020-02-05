#Unix makefile for c-files

# put the name of the target program here
TARGET = vor2fe
# the list of source files
SRC = vor2fe.c

# some definitions
RP_DIR = $(shell . /etc/environ.sh; use -e -r rappture; echo $${RAPPTURE_PATH})
CFLAGS=  -Wall -I$(RP_DIR)/include
#LDFLAGS = -L$(RP_DIR)/lib -lrappture -lm
LDFLAGS = -lm
CC = gcc

all:	$(TARGET)

$(TARGET): $(SRC:.c=.o)
	$(CC) -o $(TARGET) $(SRC:.c=.o) $(LDFLAGS)

install: all
	mv $(TARGET) ../bin

clean:	
	@rm -f $(SRC:.c=.o)

distclean: clean
	rm -f ../bin/$(TARGET)

.PHONY: all install clean distclean
