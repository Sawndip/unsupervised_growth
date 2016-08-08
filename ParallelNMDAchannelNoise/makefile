CC = mpic++
CFLAGS = -std=c++0x -g -w
objects = HH2_final.o HH2_final_pool.o HHI_final.o HHI_final_pool.o poisson_noise.o PoolParallel.o PoolParallel_test.o

all: $(objects)
	$(CC) $(CFLAGS) -o out $(objects)

PoolParallel_test.o: PoolParallel_test.cpp PoolParallel.h poisson_noise.h
	$(CC) $(CFLAGS) -c PoolParallel_test.cpp

PoolParallel.o: PoolParallel.cpp PoolParallel.h HH2_final_pool.h HHI_final_pool.h poisson_noise.h training_current.h
	$(CC) $(CFLAGS) -c PoolParallel.cpp

HH2_final_pool.o: HH2_final.h HH2_final_pool.h HH2_final_pool.cpp
	$(CC) $(CFLAGS) -c HH2_final_pool.cpp

HHI_final_pool.o: HHI_final_pool.cpp HHI_final_pool.h HHI_final.h
	$(CC) $(CFLAGS) -c HHI_final_pool.cpp 

HH2_final.o: HH2_final.cpp HH2_final.h poisson_noise.h HHI_final.h
	$(CC) $(CFLAGS) -c HH2_final.cpp

HHI_final.o: HHI_final.cpp HHI_final.h HH2_final.h poisson_noise.h
	$(CC) $(CFLAGS) -c HHI_final.cpp

poisson_noise.o: poisson_noise.cpp poisson_noise.h
	$(CC) $(CFLAGS) -c poisson_noise.cpp

clean:
	$(RM) count *.o *~


