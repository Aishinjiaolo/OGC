# setting
++      = g++
VERSION = 0.0.1
CURD    = .
FLOW    = $(CURD)/main.cpp

# my source / library
SRC_DIR = $(CURD)/src
SRC     = point.cpp segment.cpp polygon.cpp api.cpp
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

# run bundled python flow
PY_FLOW          = run.py
PY_LIB_DIR       = $(CURD)/py_lib
PY_INTERFACE_DIR = $(CURD)/swig
CAPI             = api
PY_BUNDLE_LIB    = _$(CAPI).so
PY_INTERFACE     = api.i
PY_INTERFACE_SRC = $(PY_INTERFACE:.i=_wrap.cxx)
PY_INTERFACE_OBJ = $(PY_INTERFACE_SRC:.cxx=.o)
PY_INCLUDE_DIR   = /System/Library/Frameworks/Python.framework/Versions/2.7/include/python2.7
PY_CONFIG_DIR    = /System/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/config

# build c++ library and main flow
all: $(FLOW) $(LINK)
	$(++) -Wall -I$(SRC_DIR) -o $(RUN) $(FLOW) $(LINK)

$(LINK): $(SRC_DIR)/*.cpp
	$(++) -c -fPIC -Wall $(SRC_DIR)/*.cpp
	$(++) -shared -WI -I$(SRC_DIR), -soname,$(SO) -o $(LIB) $(OBJS)
	ln -s -f $(LIB) $(LINK)

lib: $(LINK)

clean:
	rm -f $(OBJS) $(LINK) $(LIB) $(RUN)


# build gtest
gtest: $(UNITTESTS)

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


# build python flow with capi
pyflow: $(PY_FLOW) $(PY_BUNDLE_LIB) $(LINK)
	$(CURD)/$(PY_FLOW)

$(PY_BUNDLE_LIB): $(PY_INTERFACE_DIR)/$(PY_INTERFACE) $(LINK)
	swig -python -c++ $(PY_INTERFACE_DIR)/$(PY_INTERFACE)
	swig -python -c++ $(PY_INTERFACE_DIR)/polygon.i
	swig -python -c++ $(PY_INTERFACE_DIR)/segment.i
	swig -python -c++ $(PY_INTERFACE_DIR)/point.i
	$(++) -c -fPIC -Wall $(PY_INTERFACE_DIR)/*.cxx -I$(PY_INCLUDE_DIR) -I$(PY_CONFIG_DIR)
	$(++) -bundle -undefined suppress -flat_namespace $(PY_INTERFACE_OBJ)  $(LINK) -o $(PY_INTERFACE_DIR)/$(PY_BUNDLE_LIB)
	$(++) -bundle -undefined suppress -flat_namespace polygon_wrap.o  $(LINK) -o $(PY_INTERFACE_DIR)/_polygon.so
	$(++) -bundle -undefined suppress -flat_namespace segment_wrap.o  $(LINK) -o $(PY_INTERFACE_DIR)/_segment.so
	$(++) -bundle -undefined suppress -flat_namespace point_wrap.o  $(LINK) -o $(PY_INTERFACE_DIR)/_point.so

pyclean:
	rm -f $(PY_LIB_DIR)/*.pyc
	rm -f $(PY_INTERFACE_DIR)/*.py
	rm -f $(PY_INTERFACE_DIR)/*.pyc
	rm -f $(PY_INTERFACE_DIR)/*.cxx
	rm -f $(PY_INTERFACE_DIR)/*.so
	rm -f $(CURD)/*.o
	rm -f $(OBJS) $(LINK) $(LIB) $(PY_INTERFACE_OBJ) $(PY_BUNDLE_LIB)
