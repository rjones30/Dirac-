#---------------------------------------------------
CXXFLAGS      = -O4 -fPIC $(shell root-config --cflags) -I . \
                          $(shell python-config --includes)
CDBFLAGS      = -g -fPIC $(shell root-config --cflags) -I . \
                         $(shell python-config --includes)
LDFLAGS       = -g -Wl,--export-dynamic
SOFLAGS       = -shared -Wl,--export-dynamic
LD            = g++

ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)
LIBS          = $(ROOTLIBS)
EXTRA_LIBS    = -lsunmath -lcomplex
GLIBS         = $(ROOTGLIBS) -L/usr/X11R6/lib -lXpm -lX11
EXTRA_GLIBS   = -lsunmath -lcomplex

SRCS          = TThreeVectorReal.cxx \
                TThreeVectorComplex.cxx \
                TFourVectorReal.cxx \
                TFourVectorComplex.cxx \
                TLorentzTransform.cxx \
                TLorentzBoost.cxx \
                TThreeRotation.cxx \
                TPauliSpinor.cxx \
                TPauliMatrix.cxx \
                TDiracSpinor.cxx \
                TDiracMatrix.cxx \
                TPhoton.cxx \
                TLepton.cxx \
                TCrossSection.cxx \
                TCrossSection_v1.cxx

OBJS = $(foreach src, $(SRCS), $(subst cxx,o,$(src))) \
       $(foreach src, $(SRCS), $(subst .cxx,Dict.o,$(src)))

.SUFFIXES:	.so .cxx

all: libDirac.so

.cxx.o:
	@g++ -c $(CDBFLAGS) $<

.o.so:
	@echo "Building" $@
	@$(LD) $(SOFLAGS) $< -o $@
	@echo "done"

$(PROGRAM): main.o $(OBJS)
	@echo "Linking $(PROGRAM) ..."
	@$(LD) $(LDFLAGS) main.o $(OBJS) $(GLIBS) -o $(PROGRAM)
	@echo "done"

debug.o: Compton.C

debug: debug.o $(OBJS)
	@echo "Linking debug ..."
	@$(LD) $(LDFLAGS) $< $(OBJS) $(GLIBS) -o $@
	@echo "done"

clean:
	@rm -f $(OBJS) core.* *Dict.* *.o *_rdict.pcm *.so *.d 

libDirac.so: $(OBJS) python_bindings.o
	@echo "Building shared library ..."
	@$(LD) $(SOFLAGS) -Wl,-soname,$@ $^ -o $@ -lboost_python
	@echo "done"

TThreeVectorRealDict.cxx: TThreeVectorReal.h TThreeVectorRealLinkDef.h
	@echo Generating $@
	@rootcling -f $@ -c $^

TThreeVectorComplexDict.cxx: TThreeVectorComplex.h TThreeVectorComplexLinkDef.h
	@echo Generating $@
	@rootcling -f $@ -c $^

TFourVectorRealDict.cxx: TFourVectorReal.h TFourVectorRealLinkDef.h
	@echo Generating $@
	@rootcling -f $@ -c $^

TFourVectorComplexDict.cxx: TFourVectorComplex.h TFourVectorComplexLinkDef.h
	@echo Generating $@
	@rootcling -f $@ -c $^

TLorentzTransformDict.cxx: TLorentzTransform.h TLorentzTransformLinkDef.h
	@echo Generating $@
	@rootcling -f $@ -c $^

TLorentzBoostDict.cxx: TLorentzBoost.h TLorentzBoostLinkDef.h
	@echo Generating $@
	@rootcling -f $@ -c $^

TThreeRotationDict.cxx: TThreeRotation.h TThreeRotationLinkDef.h
	@echo Generating $@
	@rootcling -f $@ -c $^

TPauliSpinorDict.cxx: TPauliSpinor.h TPauliSpinorLinkDef.h
	@echo Generating $@
	@rootcling -f $@ -c $^

TPauliMatrixDict.cxx: TPauliMatrix.h TPauliMatrixLinkDef.h
	@echo Generating $@
	@rootcling -f $@ -c $^

TDiracSpinorDict.cxx: TDiracSpinor.h TDiracSpinorLinkDef.h
	@echo Generating $@
	@rootcling -f $@ -c $^

TDiracMatrixDict.cxx: TDiracMatrix.h TDiracMatrixLinkDef.h
	@echo Generating $@
	@rootcling -f $@ -c $^

TPhotonDict.cxx: TPhoton.h TPhotonLinkDef.h
	@echo Generating $@
	@rootcling -f $@ -c $^

TLeptonDict.cxx: TLepton.h TLeptonLinkDef.h
	@echo Generating $@
	@rootcling -f $@ -c $^

TCrossSectionDict.cxx: TCrossSection.h TCrossSectionLinkDef.h
	@echo Generating $@
	@rootcling -f $@ -c $^

TCrossSection_v1Dict.cxx: TCrossSection_v1.h TCrossSection_v1LinkDef.h
	@echo Generating $@
	@rootcling -f $@ -c $^

TThreeVectorReal.o:	 TThreeVectorReal.h TThreeVectorReal.cxx
TThreeVectorComplex.o:	 TThreeVectorComplex.h TThreeVectorComplex.cxx \
			 TThreeVectorReal.h
TFourVectorReal.o:	 TFourVectorReal.h TFourVectorReal.cxx \
			 TThreeVectorReal.h
TFourVectorComplex.o:	 TFourVectorComplex.h TFourVectorComplex.cxx \
			 TThreeVectorComplex.h TThreeVectorReal.h \
			 TFourVectorReal.h
TLorentzTransform.o:	 TLorentzTransform.h TLorentzTransform.cxx \
			 TLorentzBoost.h TThreeRotation.h \
			 TThreeVectorReal.h TThreeVectorComplex.h \
			 TFourVectorReal.h TFourVectorComplex.h
TLorentzBoost.o:	 TLorentzBoost.h TLorentzBoost.cxx \
			 TThreeRotation.h TLorentzTransform.h \
			 TThreeVectorReal.h TThreeVectorComplex.h \
			 TFourVectorReal.h TFourVectorComplex.h
TThreeRotation.o:	 TThreeRotation.h TThreeRotation.cxx \
			 TLorentzTransform.h \
			 TThreeVectorReal.h TThreeVectorComplex.h \
			 TFourVectorReal.h TFourVectorComplex.h
TPauliSpinor.o:		 TPauliSpinor.h TPauliSpinor.cxx \
			 TPauliMatrix.h \
			 TLorentzTransform.h TLorentzBoost.h TThreeRotation.h \
			 TThreeVectorComplex.h TThreeVectorReal.h \
			 TFourVectorComplex.h TFourVectorReal.h
TPauliMatrix.o:		 TPauliMatrix.h TPauliMatrix.cxx \
			 TPauliSpinor.h \
			 TLorentzTransform.h TLorentzBoost.h TThreeRotation.h \
			 TThreeVectorComplex.h  TThreeVectorReal.h \
			 TFourVectorComplex.h  TFourVectorReal.h
TDiracSpinor.o:		 TDiracSpinor.h TDiracSpinor.cxx \
			 TDiracMatrix.h \
			 TPauliSpinor.h TPauliMatrix.h \
			 TLorentzTransform.h TLorentzBoost.h TThreeRotation.h \
			 TFourVectorComplex.h TFourVectorReal.h \
			 TThreeVectorComplex.h TThreeVectorReal.h
TDiracMatrix.o:		 TDiracMatrix.h TDiracMatrix.cxx \
			 TDiracSpinor.h \
			 TPauliSpinor.h TPauliMatrix.h \
			 TLorentzTransform.h TLorentzBoost.h TThreeRotation.h \
			 TFourVectorComplex.h TFourVectorReal.h \
			 TThreeVectorComplex.h TThreeVectorReal.h
TPhoton.o:		 TPhoton.h TPhoton.cxx \
			 TPauliMatrix.h TPauliSpinor.h \
			 TLorentzTransform.h TLorentzBoost.h TThreeRotation.h \
			 TThreeVectorComplex.h  TThreeVectorReal.h \
			 TFourVectorComplex.h  TFourVectorReal.h
TLepton.o:		 TLepton.h TLepton.cxx \
			 TDiracMatrix.h TDiracSpinor.h \
			 TPauliSpinor.h TPauliMatrix.h \
			 TLorentzTransform.h TLorentzBoost.h TThreeRotation.h \
			 TFourVectorComplex.h TFourVectorReal.h \
			 TThreeVectorComplex.h TThreeVectorReal.h
TCrossSection.o:	 TCrossSection.h TCrossSection.cxx \
			 TLepton.h TPhoton.h \
			 TDiracMatrix.h TDiracSpinor.h \
			 TPauliSpinor.h TPauliMatrix.h \
			 TLorentzTransform.h TLorentzBoost.h TThreeRotation.h \
			 TFourVectorComplex.h TFourVectorReal.h \
			 TThreeVectorComplex.h TThreeVectorReal.h
TCrossSection_v1.o:	 TCrossSection_v1.h TCrossSection_v1.cxx \
			 TLepton.h TPhoton.h \
			 TDiracMatrix.h TDiracSpinor.h \
			 TPauliSpinor.h TPauliMatrix.h \
			 TLorentzTransform.h TLorentzBoost.h TThreeRotation.h \
			 TFourVectorComplex.h TFourVectorReal.h \
			 TThreeVectorComplex.h TThreeVectorReal.h

#---------------------------------------------------
