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

all:  $(EXECUTABLE)

$(EXECUTABLE): $(EXECUTABLE).cc 
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS) $(ROOTLIBS) $(ROOTFLAGS) -I.

#cleaning options
.PHONY: clean cleanall
clean:
	rm -f $(OBJECTS) && rm -f $(EXECUTABLE)
