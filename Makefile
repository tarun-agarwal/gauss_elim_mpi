program.out: main.o tools.o GEserial.o GEmpi.o
	mpicc -O0 -o program.out main.o tools.o GEserial.o GEmpi.o

main.o: main.c
	gcc -O0 -c main.c

tools.o: tools.c tools.h
	gcc -O0 -c tools.c

GEserial.o: GEserial.c GEserial.h
	gcc -O0 -c GEserial.c

GEmpi.o: GEmpi.c GEmpi.h
	gcc -O0 -c GEmpi.c

clean:
	rm *.o program.out
