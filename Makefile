CXX = g++
CXXFLAGS = -O3 -std=c++20 -Wall -Wextra -Wpedantic

wave: main.cpp solution.cpp solution.hpp
	$(CXX) $(CXXFLAGS) -o wave main.cpp solution.cpp

clean:
	rm -f wave

.PHONY: clean
