# setting
++      = g++
VERSION = 0.0.1
CURD    = .
FLOW    = $(CURD)/main.cpp

# my source / library
SRC_DIR = $(CURD)/src
SRC     = point.cpp segment.cpp polygon.cpp
OBJS    = $(SRC:.cpp=.o)
LINK    = lib_$(VERSION).so
LIB     = lib-$(VERSION).so
SO      = lib.so.$(VERSION)
RUN     = run

# my unit test
UNITTESTS = unittest
TEST_DIR  = $(CURD)/test
TEST_SRCS = point_unittest.cc segment_unittest.cc polygon_unittest.cc
TEST_OBJS = $(TEST_SRCS:.cc=.o)

# gtest framework
GTEST_DIR = /Users/ktlo/work/gtest/gtest-1.6.0
CPPFLAGS += -I$(GTEST_DIR)/include
CXXFLAGS += -g -Wall -Wextra

GTEST_HEADERS = $(GTEST_DIR)/include/gtest/*.h \
                $(GTEST_DIR)/include/gtest/internal/*.h
GTEST_SRCS_   = $(GTEST_DIR)/src/*.cc $(GTEST_DIR)/src/*.h $(GTEST_HEADERS)


all: $(FLOW) $(LINK)
	$(++) -Wall -I$(SRC_DIR) -o $(RUN) $(FLOW) $(LINK)

$(LINK): $(SRC_DIR)/*.cpp
	$(++) -c -fPIC -Wall $(SRC_DIR)/*.cpp
	$(++) -shared -WI -I$(SRC_DIR), -soname,$(SO) -o $(LIB) $(OBJS)
	ln -s -f $(LIB) $(LINK)

lib: $(LINK)

clean:
	rm -f $(OBJS) $(LINK) $(LIB) $(RUN)


test: $(UNITTESTS)

tclean:
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

$(OBJS): $(SRC_DIR)/*.cpp $(SRC_DIR)/*.h $(GTEST_HEADERS)
	$(++) $(CPPFLAGS) $(CXXFLAGS) -c $(SRC_DIR)/*.cpp

$(TEST_OBJS): $(TEST_DIR)/*.cc $(SRC_DIR)/*.h $(GTEST_HEADERS)
	$(++) $(CPPFLAGS) $(CXXFLAGS) -I$(SRC_DIR) -c $(TEST_DIR)/*.cc

$(UNITTESTS): $(OBJS) $(TEST_OBJS) gtest_main.a
	$(++) $(CPPFLAGS) $(CXXFLAGS) -lpthread $^ -o $@
