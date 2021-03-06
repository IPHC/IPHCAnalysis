################################################################################
#
# CondFormats/JetMETObjects Makefile (for standalone use outside CMSSW/SCRAM)
# ---------------------------------------------------------------------------
#
# INSTRUCTIONS:
# =============
# setenv ROOTSYS /path/to/root
# setenv PATH $ROOTSYS/bin:${PATH}
# setenv LD_LIBRARY_PATH $ROOTSYS/lib
#
# mkdir standalone; cd standalone
# setenv STANDALONE_DIR $PWD
# setenv PATH $STANDALONE_DIR/bin:${PATH}
# setenv LD_LIBRARY_PATH $STANDALONE_DIR/lib:${LD_LIBRARY_PATH}
# cvs co -d JetMETObjects CMSSW/CondFormats/JetMETObjects
# cd JetMETObjects
# make
#
# [you might want to stick these into e.g. $STANDALONE_DIR/setup.[c]sh]
#
#             07/11/2009 Philipp Schieferdecker <philipp.schieferdecker@cern.ch>
################################################################################

ifeq ($(STANDALONE_DIR),)
	standalone_dir:=./
	export STANDALONE_DIR:=$(standalone_dir)
endif


TMPDIR       = $(STANDALONE_DIR)/src
LIBDIR       = $(STANDALONE_DIR)/
BINDIR       = $(STANDALONE_DIR)/bin



CXX          = g++


ROOTCXXFLAGS = $(shell $(ROOTSYS)/bin/root-config --cflags)
CXXFLAGS     = -O3 -Wall -fPIC -DSTANDALONE -I. $(ROOTCXXFLAGS)

ROOTLIBS     = $(shell $(ROOTSYS)/bin/root-config --libs)

OBJS         = $(TMPDIR)/JetCorrectorParameters.o \
	       $(TMPDIR)/JetResolution.o \
	       $(TMPDIR)/SimpleJetCorrectionUncertainty.o \
	       $(TMPDIR)/JetCorrectionUncertainty.o \
	       $(TMPDIR)/SimpleJetCorrector.o \
               $(TMPDIR)/FactorizedJetCorrector.o \
               $(TMPDIR)/JetMETObjects_dict.o

LIB          = libJetCorrections.so


all: setup $(OBJS) lib

setup:
	mkdir -p $(TMPDIR)
	mkdir -p $(LIBDIR)
	mkdir -p $(BINDIR)

lib: $(OBJS)
	$(CXX) $(CXXFLAGS) -shared $(OBJS) $(ROOTLIBS) -o $(LIBDIR)/$(LIB)

clean:
	rm -rf $(OBJS) $(LIBDIR)/$(LIB) \
	       $(TMPDIR)/JetMETObjects_dict.h $(TMPDIR)/JetMETObjects_dict.cc


################################################################################
# $(OBJS)
################################################################################

$(TMPDIR)/JetResolution.o: interface/JetResolution.h \
				    src/JetResolution.cc
	$(CXX) $(CXXFLAGS) -c src/JetResolution.cc \
	-o $(TMPDIR)/JetResolution.o 

$(TMPDIR)/JetCorrectorParameters.o: interface/JetCorrectorParameters.h \
				    src/JetCorrectorParameters.cc
	$(CXX) $(CXXFLAGS) -c src/JetCorrectorParameters.cc \
	-o $(TMPDIR)/JetCorrectorParameters.o 


$(TMPDIR)/SimpleJetCorrectionUncertainty.o: interface/SimpleJetCorrectionUncertainty.h \
				    src/SimpleJetCorrectionUncertainty.cc
	$(CXX) $(CXXFLAGS) -c src/SimpleJetCorrectionUncertainty.cc \
	-o $(TMPDIR)/SimpleJetCorrectionUncertainty.o 

$(TMPDIR)/JetCorrectionUncertainty.o: interface/JetCorrectionUncertainty.h \
				    src/JetCorrectionUncertainty.cc
	$(CXX) $(CXXFLAGS) -c src/JetCorrectionUncertainty.cc \
	-o $(TMPDIR)/JetCorrectionUncertainty.o 

$(TMPDIR)/SimpleJetCorrector.o: interface/SimpleJetCorrector.h \
				    src/SimpleJetCorrector.cc
	$(CXX) $(CXXFLAGS) -c src/SimpleJetCorrector.cc \
	-o $(TMPDIR)/SimpleJetCorrector.o 

$(TMPDIR)/FactorizedJetCorrector.o: interface/FactorizedJetCorrector.h \
				    src/FactorizedJetCorrector.cc
	$(CXX) $(CXXFLAGS) -c src/FactorizedJetCorrector.cc \
	-o $(TMPDIR)/FactorizedJetCorrector.o 


$(TMPDIR)/JetMETObjects_dict.o: $(TMPDIR)/JetMETObjects_dict.cc
	$(CXX) $(CXXFLAGS) -I$(TMPDIR) -c $(TMPDIR)/JetMETObjects_dict.cc \
	-o $(TMPDIR)/JetMETObjects_dict.o

$(TMPDIR)/JetMETObjects_dict.cc: interface/JetCorrectorParameters.h \
			         interface/SimpleJetCorrectionUncertainty.h \
			         interface/JetCorrectionUncertainty.h \
				 interface/SimpleJetCorrector.h \
				 interface/FactorizedJetCorrector.h \
				 interface/JetResolution.h \
				 interface/Linkdef.h
	rm -rf $(TMPDIR)/JetMETObjects_dict.h
	rm -rf $(TMPDIR)/JetMETObjects_dict.cc
	$(ROOTSYS)/bin/rootcint -f $(TMPDIR)/JetMETObjects_dict.cc \
	-c -I$(TMPDIR) \
	interface/JetCorrectorParameters.h \
	interface/JetResolution.h \
	interface/SimpleJetCorrector.h \
	interface/FactorizedJetCorrector.h \
	interface/Linkdef.h
