CC = g++
CFLAGS = -Wall -O3 --std=c++11

Monte_Carlo_ISING_1D.exe : Monte_Carlo_ISING_1D.o random.o
	$(CC) random.o Monte_Carlo_ISING_1D.o -o Monte_Carlo_ISING_1D.exe
Monte_Carlo_ISING_1D.o : Monte_Carlo_ISING_1D.cpp
	$(CC) -c Monte_Carlo_ISING_1D.cpp -o Monte_Carlo_ISING_1D.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

clean :
	./clean.sh
	rm *.o Monte_Carlo_ISING_1D.exe eliminami.dat risultati/* risultati/h_0.02/* risultati/h_0.0/* 

esegui_h_0.0 :
	python script_h_0.py

esegui_h_0.02 :
	python script_h_0.02.py

make esegui:
	python script_h_0.py
	python script_h_0.02.py

