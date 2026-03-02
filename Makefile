#object_files (testEffi.o)
OBJECTS := $(wildcard *.o)

#root_stuff (root libraries and needed root options)
ROOTLIBS := $(shell root-config --glibs)
ROOTFLAGS := $(shell root-config --cflags --libs)  -lRooFit -lRooFitCore -lGenVector -lMathCore -lTMVA -l TMVAGui 
ROOTCINT := $(shell which rootcint)

#compiling options
DEBUGFLAGS := -O3 -Wall 
CXXFLAGS := $(DEBUGFLAGS) 

#exe_files
EXECUTABLE         := tmva-test-allYears
EXECUTABLE2        := testPlotData
EXECUTABLE3        := signifPlot
EXECUTABLE4        := signifPlotData



#all:  $(EXECUTABLE) $(EXECUTABLE1) $(EXECUTABLE2)
all:  $(EXECUTABLE) $(EXECUTABLE2) $(EXECUTABLE3) $(EXECUTABLE4) 

$(EXECUTABLE): $(EXECUTABLE).cc 
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS) $(ROOTLIBS) $(ROOTFLAGS) -I.

#$(EXECUTABLE1): $(EXECUTABLE1).cc 
#	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS) $(ROOTLIBS) $(ROOTFLAGS) -I.

$(EXECUTABLE2): $(EXECUTABLE2).cc 
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS) $(ROOTLIBS) $(ROOTFLAGS) -I.
	
$(EXECUTABLE3): $(EXECUTABLE3).cc 
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS) $(ROOTLIBS) $(ROOTFLAGS) -I.

$(EXECUTABLE4): $(EXECUTABLE4).cc 
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS) $(ROOTLIBS) $(ROOTFLAGS) -I.
		

#cleaning options
.PHONY: clean cleanall
#clean:
#	rm -f $(OBJECTS) && rm -f $(EXECUTABLE) $(EXECUTABLE1) $(EXECUTABLE2)
clean:
	rm -f $(OBJECTS) && rm -f $(EXECUTABLE) $(EXECUTABLE2) $(EXECUTABLE3) $(EXECUTABLE4) 
