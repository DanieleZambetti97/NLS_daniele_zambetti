CC = g++
CFLAGS = -Wall -O3 --std=c++11

Monte_Carlo_NVT.exe : Monte_Carlo_NVT.o random.o
	$(CC) random.o Monte_Carlo_NVT.o -o Monte_Carlo_NVT.exe
Monte_Carlo_NVT.o : Monte_Carlo_NVT.cpp
	$(CC) -c Monte_Carlo_NVT.cpp -o Monte_Carlo_NVT.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o Monte_Carlo_NVT.exe risultati/*

esegui:
	cp input.solid input.dat
	./Monte_Carlo_NVT.exe > eliminami.dat
	cp input.liquid input.dat
	./Monte_Carlo_NVT.exe > eliminami.dat
	cp input.gas input.dat
	./Monte_Carlo_NVT.exe > eliminami.dat
