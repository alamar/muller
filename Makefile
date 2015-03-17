#TARGET=$(shell basename `pwd`)
TARGET = muller
#SOURCES=$(wildcard *.cpp)
# SOURCES= main.cpp organism.cpp world.cpp
SOURCES= organism.cpp world.cpp
OBJECTS=$(SOURCES:%.cpp=%.o)

SWIGFILE = $(TARGET).i
PYWRAPCPP= $(TARGET)_wrap.cpp
PYWRAPO=$(PYWRAPCPP:%.cpp=%.o)
PYTARGET =_$(TARGET).so
PYSOURCES= $(SOURCES) $(PYWRAPCPP)
PYOBJECTS= $(PYSOURCES:%.cpp=%.o)


# CXXFLAGS+= -O3
CXXFLAGS+= -std=c++11 -O3
CXXFLAGS+= -fPIC -I/usr/include/python2.7
#CFLAGS+=$(shell pkg-config --cflags libxslt sqlite3)
#LDFLAGS+=$(shell pkg-config --libs ncurses)
# LDFLAGS+= -lncurses -lblas -g -pg
LDFLAGS+= -g -pg

# all: $(TARGET) $(PYTARGET)
all: $(PYTARGET)

$(OBJECTS): $(SOURCES)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(LDFLAGS) $(OBJECTS) $(LOADLIBES) $(LDLIBS)

clean:
	$(RM) $(PYOBJECTS) $(PYTARGET) $(TARGET) $(PYWRAPCPP)

$(PYWRAPCPP): $(SWIGFILE)
	swig -includeall -c++ -python -o $(PYWRAPCPP) $(SWIGFILE)

$(PYWRAPO): $(PYWRAPCPP)

$(PYTARGET): $(PYOBJECTS) $(PYWRAPO)
	$(CXX) -shared $(CXXFLAGS) -o $(PYTARGET) $(LDFLAGS) $(PYOBJECTS) $(LOADLIBES) $(LDLIBS)

.PHONY: all clean
