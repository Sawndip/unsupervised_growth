CC = mpic++
VPATH = ../
CFLAGS = -std=c++0x -g -w
LDFLAGS = -lconfig++
objects = Configuration.o main.o

all: $(objects)
	$(CC) $(CFLAGS) $^ $(LDFLAGS) -o out

main.o: main.cpp
	$(CC) $(CFLAGS) -c $^ -o $@

Configuration.o: Configuration.cpp
	$(CC) $(CFLAGS) -c $^ -o $@

clean:
	$(RM) count *.o *~


