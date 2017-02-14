CC = mpic++
VPATH = ../
CFLAGS = -std=c++0x -g -w
objects = HH2_final.o HHI_final.o poisson_noise.o PoolParallel.o timeIt.o

all: $(objects)
	$(CC) $(CFLAGS) $^ -o out

timeIt.o: timeIt.cpp
	$(CC) $(CFLAGS) -c $^ -o $@

PoolParallel.o: PoolParallel.cpp
	$(CC) $(CFLAGS) -c $^ -o $@

HH2_final.o: HH2_final.cpp
	$(CC) $(CFLAGS) -c $^ -o $@

HHI_final.o: HHI_final.cpp
	$(CC) $(CFLAGS) -c $^ -o $@

poisson_noise.o: poisson_noise.cpp
	$(CC) $(CFLAGS) -c $^ -o $@

clean:
	$(RM) count *.o *~


