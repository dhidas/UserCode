# define compile command

CC = g++


CFLAGS  = -Wall -g -Wno-deprecated 
LDFLAGS = -g
LDFLAGS = -L/usr/local/lib --allow-shlib-undefined
LDFLAGS = -L/usr/local/lib 
SOFLAGS = -shared -g

ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs)

CFLAGS       += $(ROOTCFLAGS)

FOBJECTS=BinaryFileReader.o PHCalibration.o\
	 LangauFitter.o ConfigReader.o  
	 # RocGeometry.o 
	 
GOBJECTS=pulseHeightSpectrum.o LangauFitter.o
	 
HOBJECTS=PHCalibration.o PHCalibrationFit.o	 	 

.cc.o:
	$(CC) $(CFLAGS) -c $<


all: convert_to_tree

noise: noise.cxx
	$(CC) $(CFLAGS) $(ROOTGLIBS) noise.cxx -o noise
	
convert_to_tree: convert_to_tree.cxx $(FOBJECTS) 
	$(CC) $(CFLAGS) -I. $(LDFLAGS) $(ROOTGLIBS) convert_to_tree.cxx -o convert_to_tree $(FOBJECTS)

raw: raw.cxx
	$(CC) raw.cxx -o raw	
	
pulseHeightMain: pulseHeightMain.cxx $(GOBJECTS) 
	$(CC) $(CFLAGS) -I. $(LDFLAGS) $(ROOTGLIBS) pulseHeightMain.cxx -o pulseHeightMain $(GOBJECTS)	
	
fitParDis: fitParDis.cxx	$(HOBJECTS) 
	$(CC) $(CFLAGS) -I. $(LDFLAGS) $(ROOTGLIBS) fitParDis.cxx -o fitParDis $(HOBJECTS)	
	
noiseCompare: noiseCompare.cxx
	$(CC) $(CFLAGS) $(ROOTGLIBS) noiseCompare.cxx -o noiseCompare
	
PHvsVcal: PHvsVcal.cxx
	$(CC) $(CFLAGS) $(ROOTGLIBS) PHvsVcal.cxx -o PHvsVcal
	
charge_vs_bias: charge_vs_bias.cxx
	$(CC) $(CFLAGS)  $(LDFLAGS) $(ROOTGLIBS) charge_vs_bias.cxx -o charge_vs_bias			
	
clean:	
	rm -f *.o charge_vs_bias PHvsVcal noiseCompare fitParDis pulseHeightMain raw convert_to_tree noise 
	#rm ViewerDict.cc ViewerDict.h

