include Makefile.arch

#
# stuff to make
#
SOURCES=$(wildcard *.cc)
OBJECTS=$(SOURCES:.cc=.o)
HEADERS=$(SOURCES:.cc=.h)
LIB=rooutil.so
EXTRACXXFLAGS= $(ROOTCFLAGS) -std=c++20 -Wall -Wno-unused-function -g -O2 -fPIC -fno-var-tracking

#
# how to make it
#

%.o: %.cc
	$(CXX) $(CXXFLAGS) $(EXTRACXXFLAGS) $< -c -o $@

$(LIB): $(OBJECTS) $(SOURCES) $(HEADERS)
	$(LD) $(LDFLAGS) -fPIC -ITMultiDrawTreePlayer -Wunused-variable $(SOFLAGS) $(OBJECTS) $(ROOTLIBS) -lTMVA -lEG -lGenVector -lXMLIO -lMLP -lTreePlayer -o $@
	ln -sf rooutil.so librooutil.so

all: $(LIB)
clean:
	rm -f *.o \
	rm -f *.d \
	rm -f *.so \
