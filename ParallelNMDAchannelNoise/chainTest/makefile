CC = mpic++
CFLAGS = -std=c++0x -g -w
VPATH = ../
objects = HH2_final.o HH2_final_pool.o HHI_final.o HHI_final_pool.o poisson_noise.o PoolParallel.o chainTest.o

all: $(objects)
	$(CC) $(CFLAGS) $^ -o chainTest

chainTest.o: chainTest.cpp
	$(CC) $(CFLAGS) -c  $^ -o $@

PoolParallel.o: PoolParallel.cpp
	$(CC) $(CFLAGS) -c  $^ -o $@

HH2_final_pool.o: HH2_final_pool.cpp
	$(CC) $(CFLAGS) -c  $^ -o $@

HHI_final_pool.o: HHI_final_pool.cpp
	$(CC) $(CFLAGS) -c $^ -o $@

HH2_final.o: HH2_final.cpp
	$(CC) $(CFLAGS) -c  $^ -o $@

HHI_final.o: HHI_final.cpp
	$(CC) $(CFLAGS) -c  $^ -o $@

poisson_noise.o: poisson_noise.cpp
	$(CC) $(CFLAGS) -c  $^ -o $@

clean:
	$(RM) count *.o *~


