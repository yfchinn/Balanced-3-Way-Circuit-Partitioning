CXX = g++
CXXFLAGS = -O3 -Wall -std=c++14

all: Lab1

Lab1: main.cpp
	$(CXX) $(CXXFLAGS) main.cpp -o Lab1

clean:
	rm -f Lab1