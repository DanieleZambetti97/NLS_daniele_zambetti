CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o datablocking.o
	$(CC) random.o main.o datablocking.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
datablocking.o : datablocking.cpp datablocking.h
	$(CC) -c datablocking.cpp -o datablocking.o $(CFLAGS)


esegui:
	./main.exe

clean :
	rm *.o main.exe seed.out
