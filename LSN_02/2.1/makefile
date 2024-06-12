CC = g++
CFLAGS = -Wall -O3 -std=c++11

main.exe : main.o random.o functions.o lib.o integral.o
	$(CC) random.o functions.o integral.o lib.o main.o -o main.exe 

main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

functions.o : functions.cpp functions.h
	$(CC) -c functions.cpp -o functions.o $(CFLAGS)

lib.o : lib.cpp lib.h
	$(CC) -c lib.cpp -o lib.o $(CFLAGS)

integral.o : integral.cpp integral.h
	$(CC) -c integral.cpp -o integral.o $(CFLAGS)

clean :
	rm *.o main.exe seed.out

run : main.exe
	./main.exe