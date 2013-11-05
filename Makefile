test: a
	./a

a: a.cpp array.h input-class.h
	g++ -g -O0 -Wall -D__CHECK -std=c++0x a.cpp -o a
# -pedantic

# Try compiling with http://en.wikipedia.org/wiki/Clang .
