CXX = g++
# FLAGS = -std=c++1z -O0 -g
FLAGS = -std=c++1z -Ofast -march=native
INCLUDES = 
LIBS = 
TARGET = main
SRCS = main.cc delaunay.cc
OBJS = $(patsubst %.cc,%.o,$(SRCS))

.cc.o:
	$(CXX) $(FLAGS) $(INCLUDES) -c $< -o $@

$(TARGET): $(OBJS)
		$(CXX) $(FLAGS) $(FLAG) -o $@ $(OBJS) $(LIBS)

clean:
	rm $(TARGET) $(OBJS)
