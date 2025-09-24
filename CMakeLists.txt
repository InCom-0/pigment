cmake_minimum_required(VERSION 4.0.0)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_SCAN_FOR_MODULES OFF)
project(spectral_cpp)


add_executable(spectral_cpp)
target_sources(spectral_cpp PRIVATE src/main.cpp)
target_link_libraries(spectral_cpp PRIVATE "-lstdc++exp")

if(MINGW AND (CMAKE_BUILD_TYPE STREQUAL "Release"))
    target_link_options(Scratch PUBLIC -static -static-libgcc -static-libstdc++)
endif()

target_sources(spectral_cpp
    PUBLIC
        FILE_SET pub_headers
        TYPE HEADERS
        BASE_DIRS
                include)

target_sources(spectral_cpp
    PRIVATE
        FILE_SET priv_headers
        TYPE HEADERS
        BASE_DIRS
                src/private_inc)