TARGET=$(shell basename `pwd`)
#SOURCES=$(wildcard *.cpp)
SOURCES= main.cpp
OBJECTS=$(SOURCES:%.cpp=%.o)

# CXXFLAGS+= -O3
CXXFLAGS+= -std=c++11 -O3
#CFLAGS+=$(shell pkg-config --cflags libxslt sqlite3)
#LDFLAGS+=$(shell pkg-config --libs ncurses)
# LDFLAGS+= -lncurses -lblas -g -pg
LDFLAGS+= -g -pg

all: $(TARGET)

$(OBJECTS): $(SOURCES)

$(TARGET): $(OBJECTS)
	$(CXX) -o $(TARGET) $(LDFLAGS) $(OBJECTS) $(LOADLIBES) $(LDLIBS)

clean:
	$(RM) $(OBJECTS) $(TARGET)

.PHONY: all clean