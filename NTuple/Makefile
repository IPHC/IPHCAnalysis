ObjSuf        = o
SrcSuf        = cc
ExeSuf        = 
DllSuf        = so
OutPutOpt     = -o
HeadSuf       = h

PWD           = $(shell pwd)
ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs) -lGenVector -lMinuit -lMinuit2 -lTMVA

# Compiler and linker
CXX           = g++
LD            = g++

# CXX flags
DEFINES       = -DNO_ORCA_CLASSES
CXXFLAGS      = -O -Wall -fPIC $(DEFINES)
CXXFLAGS     += $(ROOTCFLAGS) -I./ -I../../IPHCDataFormat/ -I./LHAPDF/include/ -I/opt/exp_soft/cms/slc5_amd64_gcc462/external/boost/1.50.0-cms/include

# Libs

LDFLAGS       = -g -O -Wall -fPIC
LIBS          = $(ROOTLIBS) -lEG
LIBS         += -L$(PWD)/LHAPDF/lib/ -lLHAPDF -lyaml-cpp 
LIBS         += -L$(PWD)/JetCorrections/ -lJetCorrections
SOFLAGS       = -shared

#------------------------------------------------------------------------------
SOURCES       = $(wildcard Selection/src/*.cc BckgdEstimation/src/*.cc Plots/src/*.cc Tools/src/*.cc tinyxml/*.cc EffEstimation/src/*.cc EventReco/src/*.cc BTagReshaping/src/*.cc)
HEADERS       = $(wildcard Selection/interface/*.h BckgdEstimation/interface/*.h Plots/interface/*.h Tools/interface/*.h tinyxml/*.h EffEstimation/interface/*.h EventReco/interface/*.h BTagReshaping/interface.h)
OBJECTS       = $(SOURCES:.$(SrcSuf)=.$(ObjSuf))
DEPENDS       = $(SOURCES:.$(SrcSuf)=.d)
SOBJECTS      = $(SOURCES:.$(SrcSuf)=.$(DllSuf))

all:  libNTupleAna.so
	@mkdir -p .lib
	cp  libNTupleAna.so .lib/libNTupleAna_`date +"%d-%m-%y_%H-%M-%S"`.so

testDir:
	if [ -d ~/libs ]  then echo "" else mkdir ~/libs  fi

clean:
	@echo "Cleaning..."
	make -C JetCorrections clean
	@rm -f $(OBJECTS) $(DEPENDS) *Dict.* core 

.SUFFIXES: .$(SrcSuf) .C .o .so

libNTupleAna.so: libLHAPDF libJetCorrections $(OBJECTS)
	@echo "Building libNTupleAna..."
	$(LD) $(LIBS) $(SOFLAGS) $(LDFLAGS) $(OBJECTS) -o $@

libJetCorrections: JetCorrections/libJetCorrections.so

JetCorrections/libJetCorrections.so:
	make -C JetCorrections lib

libLHAPDF: LHAPDF/lib/libLHAPDF.so

LHAPDF/lib/libLHAPDF.so:
	cd LHAPDF; python install.py
