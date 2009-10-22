CC = g++
CFLAGS = `root-config --cflags`
LIBS = `root-config --libs` -lMinuit
SCRIPTS = $(patsubst %.cc,%,$(wildcard *.cc))


all: $(SCRIPTS)

% : %.cc
	$(CC) -Wall $(CFLAGS) $(LIBS) $< -o $@





clean:
	rm -f $(SCRIPTS)

