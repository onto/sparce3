CXX=g++
CFLAGS=-c -Wall -Wextra -fopenmp -I. -O3
LDFLAGS=-fopenmp
SOURCES1=main.cpp sparce.cpp
SOURCES2=generator.cpp
OBJECTS1=$(SOURCES1:.cpp=.o)
OBJECTS2=$(SOURCES2:.cpp=.o)
EXECUTABLE1=sparce3
EXECUTABLE2=generator

all: $(SOURCES) $(EXECUTABLE1) $(EXECUTABLE2)

$(EXECUTABLE1): $(OBJECTS1)
	$(CXX) $(LDFLAGS) $(OBJECTS1) -o $@

$(EXECUTABLE2): $(OBJECTS2)
	$(CXX) $(LDFLAGS) $(OBJECTS2) -o $@

.cpp.o:
	$(CXX) $(CFLAGS) $< -o $@

clean:
	rm *.o $(EXECUTABLE1) $(EXECUTABLE2)
