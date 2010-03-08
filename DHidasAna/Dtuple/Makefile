CC = g++
LD = g++
CFLAGS = -Wall `root-config --cflags`
LIBS = `root-config --libs` -Llib -lMinuit
INCLUDE = -Iinterface -Isrc -I../../ -I.
OBJS  = $(patsubst src/%.cc,lib/%.o,$(wildcard src/*.cc))
EXECS = $(patsubst bin/%.cc,%,$(wildcard bin/*.cc))
SCRIPTS = $(patsubst scripts/%.cc,%,$(wildcard scripts/*.cc))
EXEOBJS  = $(patsubst bin/%.cc,lib/%.o,$(wildcard bin/*.cc))
LINKDEFCC  = $(patsubst src/%_LinkDef.h,src/%_dict.cc,$(wildcard src/*_LinkDef.h))
LINKDEFOBJ  = $(patsubst src/%_LinkDef.h,lib/%_dict.o,$(wildcard src/*_LinkDef.h))

all: $(LINKDEFCC) $(LINKDEFOBJ) $(OBJS) $(EXEOBJS) $(SCRIPTS) $(EXECS)

lib/%.o : src/%.cc
	$(CC) -Wall $(CFLAGS) $(INCLUDE) -c $< -o $@

lib/%.o : bin/%.cc
	$(CC) -Wall $(CFLAGS) $(INCLUDE) -c $< -o $@

lib/%.o : scripts/%.cc
	$(CC) -Wall $(CFLAGS) $(INCLUDE) -c $< -o $@

lib/%.o : dict/%.cc
	$(CC) -Wall $(CFLAGS) $(INCLUDE) -c $< -o $@

src/%_dict.cc : interface/%.h src/%_LinkDef.h
	rootcint -f $@ -c -I../../ $<

% : $(OBJS) lib/%.o
	$(LD) $(LIBS) $(OBJS) $(LINKDEFOBJ) lib/$@.o -o $@





clean:
	rm -f $(EXECS) $(SCRIPTS) lib/*.o src/*_dict.cc src/*_dict.h
