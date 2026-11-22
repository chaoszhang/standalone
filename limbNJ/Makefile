CXX      ?= g++
CXXFLAGS ?= -O2 -std=c++20 -Wall -Wextra

BINS = NJ limbNJ1.3 limbNJ1.3.1

all: $(BINS)

NJ: NJ.cpp limbNJ.hpp
	$(CXX) $(CXXFLAGS) -o $@ $<

limbNJ1.3: limbNJ1.3.cpp limbNJ.hpp
	$(CXX) $(CXXFLAGS) -o $@ $<

limbNJ1.3.1: limbNJ1.3.1.cpp limbNJ.hpp
	$(CXX) $(CXXFLAGS) -o $@ $<

clean:
	rm -f $(BINS)

.PHONY: all clean
