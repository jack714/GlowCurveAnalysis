# the compiler: gcc for C program, define as g++ for C++
CC = g++

# compiler flags:
#  -g    adds debugging information to the executable file
#  -Wall turns on most, but not all, compiler warnings
CFLAGS  = -std=c++17 -g -Wall

# the build target executable:
TARGET = GCA

all: $(TARGET)

$(TARGET): data_input.cpp smartPeakDetect.cpp FileHandler.cpp matrix_arithmetic.cpp Levenberg_Marquardt.cpp DataSmoothing.cpp CSV_iterator.cpp main.cpp FOKModel.cpp File_Manager.cpp quick_half_max.cpp background_subtraction.cpp xls_iterator.cpp
	$(CC) $(CFLAGS) -o $(TARGET) data_input.cpp smartPeakDetect.cpp matrix_arithmetic.cpp Levenberg_Marquardt.cpp DataSmoothing.cpp FileHandler.cpp CSV_iterator.cpp main.cpp FOKModel.cpp File_Manager.cpp quick_half_max.cpp background_subtraction.cpp xls_iterator.cpp

clean:
	$(RM) $(TARGET)


