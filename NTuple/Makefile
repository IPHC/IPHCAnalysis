ObjSuf        = o
SrcSuf        = cc
ExeSuf        = 
DllSuf        = so
OutPutOpt     = -o
HeadSuf       = h

PWD			  = $(shell pwd)
ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = -lGenVector -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lMinuit -lMinuit2 -lPhysics -lMathCore -lThread -pthread -lm -ldl -rdynamic -lTMVA
ROOTGLIBS     = -lGenVector -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lMinuit -lMinuit2 -lPhysics -lMathCore -lThread -pthread -lm -ldl -rdynamic -lTMVA

# Linux with egcs
DEFINES       = -DNO_ORCA_CLASSES
CXX           = g++
CXXFLAGS	  = -O -Wall -fPIC $(DEFINES)
LD		      = g++
LDFLAGS		  = -g -O -Wall -fPIC
SOFLAGS		  = -shared

CXXFLAGS	 += $(ROOTCFLAGS) -I./ -I../../IPHCDataFormat/ -I./LHAPDF/include/ -I/opt/exp_soft/cms/slc5_amd64_gcc462/external/boost/1.50.0-cms/include
LIBS		  = $(ROOTLIBS) -lEG -L$(PWD)/LHAPDF/lib/ -lLHAPDF -lyaml-cpp -L$(PWD)/JetCorrections/ -lJetCorrections
GLIBS		  = $(ROOTGLIBS)  
#------------------------------------------------------------------------------
SOURCES		  = $(wildcard Selection/src/*.cc BckgdEstimation/src/*.cc Plots/src/*.cc Tools/src/*.cc tinyxml/*.cc EffEstimation/src/*.cc EventReco/src/*.cc BTagReshaping/src/*.cc)
HEADERS		  = $(wildcard Selection/interface/*.h BckgdEstimation/interface/*.h Plots/interface/*.h Tools/interface/*.h tinyxml/*.h EffEstimation/interface/*.h EventReco/interface/*.h BTagReshaping/interface.h)
OBJECTS		  = $(SOURCES:.$(SrcSuf)=.$(ObjSuf))
DEPENDS		  = $(SOURCES:.$(SrcSuf)=.d)
SOBJECTS	  = $(SOURCES:.$(SrcSuf)=.$(DllSuf))


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

libNTupleAna.so: libJetCorrections libLHAPDF $(OBJECTS)
	@echo "Building libNTupleAna..."
	$(LD) -L${ROOTSYS}/lib $(LIBS) $(SOFLAGS) $(LDFLAGS) $(OBJECTS) -o $@

libLHAPDF: LHAPDF/lib/libLHAPDF.so
libJetCorrections: JetCorrections/libJetCorrections.so

JetCorrections/libJetCorrections.so:
	make -C JetCorrections lib

LHAPDF/lib/libLHAPDF.so:
	cd LHAPDF; python install.py
