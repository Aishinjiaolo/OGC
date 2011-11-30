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

UNITTESTS = unittest
GTEST_DIR = /Users/ktlo/work/gtest/gtest-1.6.0
TEST_DIR  = $(CURD)/test
CPPFLAGS += -I$(GTEST_DIR)/include
CXXFLAGS += -g -Wall -Wextra

GTEST_HEADERS = $(GTEST_DIR)/include/gtest/*.h \
                $(GTEST_DIR)/include/gtest/internal/*.h

GTEST_SRCS_ = $(GTEST_DIR)/src/*.cc $(GTEST_DIR)/src/*.h $(GTEST_HEADERS)

all: $(FLOW) $(LINK)
	$(++) -Wall -I$(SRC_DIR) -o $(RUN) $(FLOW) $(LINK)

$(LINK): $(SRC_DIR)/*.cpp
	$(++) -c -fPIC -Wall $(SRC_DIR)/*.cpp
	$(++) -shared -WI -I$(SRC_DIR), -soname,$(SO) -o $(LIB) $(OBJS)
	ln -s -f $(LIB) $(LINK)

lib: $(LINK)

clean:
	rm -f $(OBJS) $(LINK) $(LIB) $(RUN)

gtest: $(UNITTESTS)

gtclean:
	rm -f $(UNITTESTS) gtest.a gtest_main.a *.o

gtest-all.o: $(GTEST_SRCS_)
	$(++) $(CPPFLAGS) -I$(GTEST_DIR) $(CXXFLAGS) -c \
            $(GTEST_DIR)/src/gtest-all.cc

gtest_main.o: $(GTEST_SRCS_)
	$(++) $(CPPFLAGS) -I$(GTEST_DIR) $(CXXFLAGS) -c \
            $(GTEST_DIR)/src/gtest_main.cc

gtest.a: gtest-all.o
	$(AR) $(ARFLAGS) $@ $^

gtest_main.a: gtest-all.o gtest_main.o
	$(AR) $(ARFLAGS) $@ $^

point.o: $(SRC_DIR)/point.cpp $(SRC_DIR)/point.h $(GTEST_HEADERS)
	$(++) $(CPPFLAGS) $(CXXFLAGS) -c $(SRC_DIR)/point.cpp

point_unittest.o: $(TEST_DIR)/point_unittest.cc \
                     $(SRC_DIR)/point.h $(GTEST_HEADERS)
	$(++) $(CPPFLAGS) $(CXXFLAGS) -I$(SRC_DIR) -c $(TEST_DIR)/point_unittest.cc

$(UNITTESTS): point.o point_unittest.o gtest_main.a
	$(++) $(CPPFLAGS) $(CXXFLAGS) -lpthread $^ -o $@
