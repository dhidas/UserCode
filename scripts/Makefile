CC = g++
CFLAGS = `root-config --cflags`
LIBS = `root-config --libs` -Llib -lMinuit
INCLUDE = -Iinclude -Idict -I.
SCRIPTS = $(patsubst %.cc,%,$(wildcard *.cc))


all: $(SCRIPTS)

% : %.cc
	$(CC) -Wall $(CFLAGS) $(LIBS) $(INCLUDE) $< -o $@





clean:
	rm -f $(SCRIPTS)

