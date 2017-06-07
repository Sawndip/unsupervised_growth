CC = mpic++
VPATH = ../
CFLAGS = -std=c++0x -g -w
LDFLAGS = -lconfig++
objects = ConfigurationNetworkGenerator.o poisson_noise.o NetworkGenerator.o main.o

all: $(objects)
	$(CC) $(CFLAGS) $^ $(LDFLAGS) -o generateNetwork

main.o: main.cpp
	$(CC) $(CFLAGS) -c $^ -o $@

NetworkGenerator.o: NetworkGenerator.cpp
	$(CC) $(CFLAGS) -c $^ -o $@

ConfigurationNetworkGenerator.o: ConfigurationNetworkGenerator.cpp
	$(CC) $(CFLAGS) -c $^ -o $@

poisson_noise.o: poisson_noise.cpp
	$(CC) $(CFLAGS) -c $^ -o $@

clean:
	$(RM) count *.o *~
