# the compiler: gcc for C program, define as g++ for C++
CC = g++

# compiler flags:
#  -g    adds debugging information to the executable file
#  -Wall turns on most, but not all, compiler warnings
CFLAGS  = -std=c++17 -g -Wall

# the build target executable:
TARGET = GCA

all: $(TARGET)

$(TARGET): src/smartPeakDetect.cpp src/FileHandler.cpp src/matrix_arithmetic.cpp src/Levenberg–Marquardt.cpp src/DataSmoothing.cpp src/CSV_iterator.cpp src/main.cpp src/FOKModel.cpp src/File_Manager.cpp src/quick_half_max.cpp src/data_input.cpp
	$(CC) $(CFLAGS) -o $(TARGET) src/smartPeakDetect.cpp src/matrix_arithmetic.cpp src/Levenberg–Marquardt.cpp src/DataSmoothing.cpp src/FileHandler.cpp src/CSV_iterator.cpp src/main.cpp src/FOKModel.cpp src/File_Manager.cpp src/quick_half_max.cpp src/data_input.cpp

clean:
	$(RM) $(TARGET)


