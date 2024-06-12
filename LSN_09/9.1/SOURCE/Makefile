CC = g++
CFLAGS = -Wall -O3 -std=c++11
AFLAGS = -larmadillo

main.exe : main.o random.o city.o tsp.o
	$(CC) random.o city.o tsp.o main.o -o main.exe $(AFLAGS)

main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

city.o : city.cpp city.h
	$(CC) -c city.cpp -o city.o $(CFLAGS)

tsp.o : tsp.cpp tsp.h
	$(CC) -c tsp.cpp -o tsp.o $(CFLAGS)

clean :
	-rm *.o main.exe seed.out

run : main.exe
	./main.exe