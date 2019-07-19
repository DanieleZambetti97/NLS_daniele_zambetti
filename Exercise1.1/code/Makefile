CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o datablocking.o random.o
	$(CC) main.o datablocking.o random.o -o main.exe $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
datablocking.o : datablocking.cpp datablocking.h
	$(CC) -c datablocking.cpp -o datablocking.o $(CFLAGS)
main.o : main.cpp 
	$(CC) -c main.cpp -o main.o $(CFLAGS)

esegui :
	./main.exe > eliminami.dat

clean :
	rm *.o main.exe seed.out Risultati*.dat
