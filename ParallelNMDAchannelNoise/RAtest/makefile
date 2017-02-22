CC = mpic++
VPATH = ../
CFLAGS = -std=c++0x -g -w
objects = HH2_final.o HHI_final.o poisson_noise.o main.o

all: $(objects)
	$(CC) $(CFLAGS) $^ $(LDFLAGS) -o out

main.o: main.cpp
	$(CC) $(CFLAGS) -c $^ -o $@

HH2_final.o: HH2_final.cpp
	$(CC) $(CFLAGS) -c $^ -o $@

HHI_final.o: HHI_final.cpp
	$(CC) $(CFLAGS) -c $^ -o $@

poisson_noise.o: poisson_noise.cpp
	$(CC) $(CFLAGS) -c $^ -o $@

clean:
	$(RM) count *.o *~
