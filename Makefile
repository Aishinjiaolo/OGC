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

# bundled python flow
PY_FLOW          = run.py
PY_LIB_DIR       = $(CURD)/py_lib
PY_INTERFACE_DIR = $(CURD)/swig
CAPI             = api polygon segment point
PY_BUNDLE        = py_bundle
PY_BUNDLE_LIB    = $(CAPI:%=_%.so)
PY_INTERFACE     = $(CAPI:=.i)
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


# build capi with swig to run python flow
pyflow: $(PY_FLOW) $(PY_BUNDLE) $(LINK)
	$(CURD)/$(PY_FLOW)

$(PY_BUNDLE): $(PY_INTERFACE_DIR)/*.i $(LINK)
	for interface in $(PY_INTERFACE) ; do \
		swig -python -c++ $(PY_INTERFACE_DIR)/$$interface ; \
	done

	for py_src in $(PY_INTERFACE_SRC) ; do \
		$(++) -c -fPIC -Wall $(PY_INTERFACE_DIR)/$$py_src -I$(PY_INCLUDE_DIR) -I$(PY_CONFIG_DIR) ; \
	done

	object_index=0 ; \
	for py_obj in $(PY_INTERFACE_OBJ) ; do \
		lib_index=0 ; \
		for py_lib in $(PY_BUNDLE_LIB) ; do \
			if [ $$object_index -eq $$lib_index ] ; then \
				$(++) -bundle -undefined suppress -flat_namespace $$py_obj $(LINK) -o $(PY_INTERFACE_DIR)/$$py_lib ; \
			fi ; \
			((lib_index = lib_index + 1)) ; \
		done ; \
		((object_index = object_index + 1)) ; \
	done ; \

pyclean:
	rm -f $(PY_LIB_DIR)/*.pyc
	rm -f $(PY_INTERFACE_DIR)/*.py
	rm -f $(PY_INTERFACE_DIR)/*.pyc
	rm -f $(PY_INTERFACE_DIR)/*.cxx
	rm -f $(PY_INTERFACE_DIR)/*.so
	rm -f $(OBJS) $(LINK) $(LIB) $(PY_INTERFACE_OBJ)
