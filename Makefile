#TARGET=$(shell basename `pwd`)
TARGET = muller
#SOURCES=$(wildcard *.cpp)
# SOURCES= main.cpp organism.cpp world.cpp
SOURCES= organism.cpp world.cpp main.cpp
OBJECTS= organism.o world.o
MAINO=main.o
HEADERS = $(wildcard *.h)

SWIGFILE = $(TARGET).i
PYWRAPCPP= $(TARGET)_wrap.cpp
PYWRAPO=$(PYWRAPCPP:%.cpp=%.o)
PYTARGET =_$(TARGET).so
PYSOURCES= $(SOURCES) $(PYWRAPCPP)
PYOBJECTS= $(OBJECTS) $(PYRWAPO)


# CXXFLAGS+= -O3
CXXFLAGS+= -g -std=c++11 -O3
CXXFLAGS+= -fPIC -I/usr/include/python2.7
#CFLAGS+=$(shell pkg-config --cflags libxslt sqlite3)
#LDFLAGS+=$(shell pkg-config --libs ncurses)
# LDFLAGS+= -lncurses -lblas -g -pg
LDFLAGS+= -g -pg

# all: $(TARGET) $(PYTARGET)
all: $(PYTARGET)

$(OBJECTS): $(SOURCES) $(HEADERS)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(LDFLAGS) $(OBJECTS) $(LOADLIBES) $(LDLIBS)

clean:
	$(RM) $(PYOBJECTS) $(PYTARGET) $(MAINO) $(TARGET) $(PYWRAPCPP)

$(PYWRAPCPP): $(SWIGFILE)
	swig -includeall -c++ -python -o $(PYWRAPCPP) $(SWIGFILE)

$(PYWRAPO): $(PYWRAPCPP)

$(PYTARGET): $(PYOBJECTS) $(PYWRAPO)
	$(CXX) -shared $(CXXFLAGS) -o $(PYTARGET) $(LDFLAGS) $(PYOBJECTS) $(LOADLIBES) $(LDLIBS)

main: $(OBJECTS) $(MAINO)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(LDFLAGS) $(OBJECTS) $(MAINO) $(LOADLIBES) $(LDLIBS)

.PHONY: all clean
