++      = g++
VERSION = 0.0.1
SRC     = point.cpp
OBJS    = $(SRC:.cpp=.o)
LINK    = lib_$(VERSION).so
LIB     = lib-$(VERSION).so
SO      = lib.so.$(VERSION)
RUN     = run

all: main.cpp $(LINK)
	$(++) -Wall -o $(RUN) main.cpp $(LINK)

$(LINK): $(SRC)
	$(++) -c -fPIC point.cpp
	$(++) -shared -WI, -soname,$(SO) -o $(LIB) $(OBJS)
	ln -s -f $(LIB) $(LINK)

clean:
	rm -f *.so *.o $(RUN)

