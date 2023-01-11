cmake_minimum_required(VERSION 3.19)
project(encryptGCN)

set(CMAKE_CXX_STANDARD 20)

add_executable(encryptGCN main.cpp util.cpp util.h mpcUtil.cpp mpcUtil.h typedef.h test.cpp test.h)

include_directories("/opt/homebrew/Cellar/gmp/6.2.1/lib")
include_directories("/opt/homebrew/Cellar/gmp/6.2.1/include")

target_link_libraries(encryptGCN "/opt/homebrew/Cellar/gmp/6.2.1/lib/libgmp.a")
target_link_libraries(encryptGCN "/opt/homebrew/Cellar/gmp/6.2.1/lib/libgmpxx.a")