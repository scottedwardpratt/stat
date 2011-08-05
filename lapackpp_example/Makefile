PREFIX = /playpen2/local
CFLAGS = -I $(PREFIX)/include/lapackpp -L $(PREFIX)/lib -llapackpp -g -Wall
CXX = g++
lapackpp_example: lapackpp_example.cxx
	$(CXX) $(CFLAGS) $< -o $@
run: lapackpp_example
	LD_LIBRARY_PATH=$(PREFIX)/lib ./lapackpp_example
