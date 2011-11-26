all: main.cpp libpoint.so
	g++ -Wall -o run main.cpp libpoint.so

libpoint.so: point.cpp
	g++ -c -fPIC point.cpp
	g++ -shared -WI, -soname,libpoint.so.1 -o libpoint.so1 point.o
	ln -s -f libpoint.so1 libpoint.so

clean:
	rm -f *.so *.o *.so* run

