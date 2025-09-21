cmake_minimum_required(VERSION 4.0)
project(MathPlot LANGUAGES CXX)

set(CMAKE_CXX_COMPILER "clang++")
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_CXX_FLAGS "-stdlib=libc++")
set(CMAKE_EXE_LINKER_FLAGS "-stdlib=libc++ -lc++abi")

set(CMAKE_STATIC_LINKER_FLAGS "-stdlib=libc++")
set(CMAKE_SHARED_LINKER_FLAGS "-stdlib=libc++ -lc++abi")
set(CMAKE_MODULE_LINKER_FLAGS "-stdlib=libc++ -lc++abi")

if(NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    message(WARNING "Not using Clang compiler!")
endif()


include_directories(./include)

set(SOURCE_LIB 
	./source/interpolation.cpp 
	./source/MathAnalys.cpp 
	./source/RungeKutt.cpp
)
add_library(MathPlot SHARED ${SOURCE_LIB})
target_link_libraries(MathPlot)
