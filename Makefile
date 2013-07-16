TARGET = libSigmaScan.so
HEADERS = CutVarTreeProducer.h SignificanceCalculator.h
OBJECTS = CutVarTreeProducer.cc SignificanceCalculator.cc

INC = -I. -I$(shell root-config --incdir) -I$(CMSSW_BASE)/src
LIBS = $(shell root-config --libs)

all: $(TARGET)

clean:
	rm $(TARGET)

$(TARGET): $(OBJECTS) $(HEADERS)
	g++ -O2 -Wall -fPIC -shared -o $(TARGET) $(INC) $(LIBS) $(OBJECTS)
