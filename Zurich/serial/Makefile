CXX=gcc
CFLAGS= -Wall -O2 -lm

SRCDIR := src
INCDIR := include

C_FILES = $(wildcard ./$(SRCDIR)/*.c)
H_FILES = $(wildcard ./$(INCDIR)/*.h)

EXEC = main

all: build

build: $(EXEC)

$(EXEC): $(C_FILES) $(H_FILES)
	$(CXX) 	$(CFLAGS) -I ./$(INCDIR) $(C_FILES) -lm $(EXEC).c -o $@

clean:
	rm -rf $(EXEC)