OBJ = main.o ConfigurationNetworkGrowth.o ConfigurationNetworkTopology.o
SRC = main.cpp ConfigurationNetworkGrowth.cpp ConfigurationNetworkTopology.cpp
CXX = g++
CFLAGS = -std=c++11 -g -Wall 
CXXFLAGS = -std=c++11 -g -Wall -fopenmp
CPPFLAGS = -I../
LDFLAGS = -lconfig++
VPATH = ../

main: $(OBJ)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) $(OBJ) -o main

%.d: %.cpp
	@set -e; rm -f $@; \
	$(CXX) -MM $(CPPFLAGS) $(CXXFLAGS) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

include $(SRC:.cpp=.d)

clean:
	rm -rf *.o
