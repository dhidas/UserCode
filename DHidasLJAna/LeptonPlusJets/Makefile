CC = g++
LD = g++
CFLAGS = -Wall `root-config --cflags`
LIBS = `root-config --libs` -Llib -lMinuit -lvectorDict
INCLUDE = -Iinterface -Idict -I../../
OBJS  = $(patsubst src/%.cc,lib/%.o,$(wildcard src/*.cc))
EXECS = $(patsubst bin/%.cc,%,$(wildcard bin/*.cc))
SCRIPTS = $(patsubst scripts/%.cc,%,$(wildcard scripts/*.cc))
EXEOBJS  = $(patsubst bin/%.cc,lib/%.o,$(wildcard bin/*.cc))


all: $(LINKDEFCC) $(LINKDEFOBJ) $(OBJS) $(EXEOBJS) $(SCRIPTS) $(EXECS)

lib/%.o : src/%.cc
	$(CC) -Wall $(CFLAGS) $(INCLUDE) -c $< -o $@

lib/%.o : bin/%.cc
	$(CC) -Wall $(CFLAGS) $(INCLUDE) -c $< -o $@

lib/%.o : scripts/%.cc
	$(CC) -Wall $(CFLAGS) $(INCLUDE) -c $< -o $@

lib/%.o : dict/%.cc
	$(CC) -Wall $(CFLAGS) $(INCLUDE) -c $< -o $@

% : $(OBJS) lib/%.o
	$(LD) $(LIBS) $(OBJS) lib/$@.o -o $@





clean:
	rm -f $(EXECS) $(SCRIPTS) lib/*.o dict/*_dict.cc dict/*_dict.h
