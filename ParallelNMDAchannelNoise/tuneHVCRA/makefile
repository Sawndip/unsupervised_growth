CC = mpic++
VPATH = ../
CFLAGS = -std=c++0x -g -w -O2
LDFLAGS = -lconfig++
objects = poisson_noise.o HHI_buffer.o HH2_buffer.o  HH2_test.o main.o

all: $(objects)
	$(CC) $(CFLAGS) $^ $(LDFLAGS) -o exec

main.o: main.cpp
	$(CC) $(CFLAGS) -c $^ -o $@

HH2_test.o: HH2_test.cpp
	$(CC) $(CFLAGS) -c $^ -o $@

HH2_buffer.o: HH2_buffer.cpp
	$(CC) $(CFLAGS) -c $^ -o $@

HHI_buffer.o: HHI_buffer.cpp
	$(CC) $(CFLAGS) -c $^ -o $@

poisson_noise.o: poisson_noise.cpp
	$(CC) $(CFLAGS) -c $^ -o $@

clean:
	$(RM) count *.o *~
