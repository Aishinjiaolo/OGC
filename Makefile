++      = g++
VERSION = 0.0.1
FLOW    = main.cpp
SRC     = point.cpp segment.cpp polygon.cpp
OBJS    = $(SRC:.cpp=.o)
LINK    = lib_$(VERSION).so
LIB     = lib-$(VERSION).so
SO      = lib.so.$(VERSION)
RUN     = run

all: $(FLOW) $(LINK)
	$(++) -Wall -o $(RUN) $(FLOW) $(LINK)

$(LINK): $(SRC)
	$(++) -c -fPIC -Wall $(SRC)
	$(++) -shared -WI, -soname,$(SO) -o $(LIB) $(OBJS)
	ln -s -f $(LIB) $(LINK)

lib: $(LINK)

clean:
	rm -f $(OBJS) $(LINK) $(LIB) $(RUN)

