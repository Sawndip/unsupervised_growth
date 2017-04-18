CC = mpic++
CFLAGS = -std=c++0x -g -Wall
#LDFLAGS = -lconfig++
LDFLAGS = -L /usr/local/lib -l config++ -Bstatic
objects = Configuration.o HH2_final.o HHI_final.o poisson_noise.o PoolParallel.o PoolParallel_test.o

all: $(objects)
	$(CC) $(CFLAGS) $(objects) $(LDFLAGS) -o out 

PoolParallel_test.o: PoolParallel_test.cpp PoolParallel.h poisson_noise.h
	$(CC) $(CFLAGS) -c PoolParallel_test.cpp

PoolParallel.o: PoolParallel.cpp PoolParallel.h HH2_final.h HHI_final.h poisson_noise.h training_current.h
	$(CC) $(CFLAGS) -c PoolParallel.cpp

HH2_final.o: HH2_final.cpp HH2_final.h poisson_noise.h HHI_final.h
	$(CC) $(CFLAGS) -c HH2_final.cpp

HHI_final.o: HHI_final.cpp HHI_final.h HH2_final.h poisson_noise.h
	$(CC) $(CFLAGS) -c HHI_final.cpp

Configuration.o: Configuration.cpp Configuration.h
	$(CC) $(CFLAGS) -c Configuration.cpp

poisson_noise.o: poisson_noise.cpp poisson_noise.h
	$(CC) $(CFLAGS) -c poisson_noise.cpp

clean:
	$(RM) count *.o *~


