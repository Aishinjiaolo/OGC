++      = g++
VERSION = 0.0.1
CURD    = .
FLOW    = $(CURD)/main.cpp

SRC_DIR = $(CURD)/src
SRC     = point.cpp segment.cpp polygon.cpp
OBJS    = $(SRC:.cpp=.o)
LINK    = lib_$(VERSION).so
LIB     = lib-$(VERSION).so
SO      = lib.so.$(VERSION)
RUN     = run

TEST_DIR = $(CURD)/test
TEST     = unittest

all: $(FLOW) $(LINK)
	$(++) -Wall -I$(SRC_DIR) -o $(RUN) $(FLOW) $(LINK)

$(LINK): $(SRC_DIR)/*.cpp
	$(++) -c -fPIC -Wall $(SRC_DIR)/*.cpp
	$(++) -shared -WI -I$(SRC_DIR), -soname,$(SO) -o $(LIB) $(OBJS)
	ln -s -f $(LIB) $(LINK)

lib: $(LINK)

clean:
	rm -f $(OBJS) $(LINK) $(LIB) $(RUN)

gtest:
	make -f $(TEST_DIR)/Makefile

gtclean:
	make -f $(TEST_DIR)/Makefile clean
