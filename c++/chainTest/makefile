OBJ = main.o HH2_buffer.o HH2_test.o HHI_buffer.o poisson_noise.o NetworkGrowthSimulator.o ConfigurationNetworkGrowth.o ConfigurationNetworkTopology.o
SRC = main.cpp HH2_buffer.cpp HH2_test.cpp HHI_buffer.cpp poisson_noise.cpp NetworkGrowthSimulator.cpp ConfigurationNetworkGrowth.cpp ConfigurationNetworkTopology.cpp
CXX = mpic++
CFLAGS = -std=c++11 -g -w 
CXXFLAGS = -std=c++11 -g -w
CPPFLAGS = -I../include
LDFLAGS = -lconfig++
VPATH = ../src/

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
