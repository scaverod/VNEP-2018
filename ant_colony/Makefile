CC=g++

boostdir = /usr/local/boost_1_53_0/ 
boostldir = /usr/local/boost_1_53_0/stage/lib
builddir = build/
testdir = test

CPPFLAGS = -W -Wall -pedantic -lm -std=c++0x -I $(boostdir) -L $(boostldir)
VPATH = src

OBJ = main.o util.o vmap.o
BOBJ = $(OBJ:%=$(builddir)%)

all:	$(builddir)vne.out

$(builddir)%.o: %.cpp %.h
	$(CC) $(CPPFLAGS) -o $@ -c $<

$(builddir)main.o: main.cpp
	$(CC) $(CPPFLAGS) -o $@ -c $<

$(builddir)vne.out:	$(BOBJ)
		$(CC) $(CPPFLOAGS) -o $@ $^

# $(builddir)main.o: main.cpp
#		$(CC) $(CPPFLAGS) -o $@ -c $<




clean:
	rm -rf build/*
	rm -rf *.out
