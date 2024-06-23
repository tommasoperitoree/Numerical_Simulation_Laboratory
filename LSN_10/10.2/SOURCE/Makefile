CC = g++
CCMPI = mpicxx
CFLAGS = -Wall -Wextra -O3 -std=c++11
AFLAGS = -larmadillo

main.exe : main.o random.o city.o tsp.o
	$(CCMPI) random.o city.o tsp.o main.o -o main.exe $(CFLAGS) $(AFLAGS)

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

clean circle :
	-rm ../OUTPUT/CIRCLE/* 
	
clean square :
	-rm ../OUTPUT/SQUARE/*

run : main.exe
	mpiexec -np 8 ./main.exe