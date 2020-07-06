
CXX := $(if $(CXX),$(CXX),g++)
CXXFLAGS=-Iseqan/include -pthread -fopenmp -O3 -std=c++14 -static-libstdc++

srcs=src/main.cpp src/utils.cpp src/io.cpp src/Settings.cpp src/Jellyfish.cpp src/Read.cpp src/Explorer.cpp src/Trail.cpp src/Trajectory.cpp

talc: $(srcs)
	$(CXX) $(CXXFLAGS) $(srcs) -o talc


