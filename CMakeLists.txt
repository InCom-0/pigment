cmake_minimum_required(VERSION 4.0.0)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_SCAN_FOR_MODULES OFF)
project(inc_pigment)


add_executable(inc_pigment)
target_sources(inc_pigment PRIVATE src/main.cpp)
target_link_libraries(inc_pigment PRIVATE "-lstdc++exp")

if(MINGW AND (CMAKE_BUILD_TYPE STREQUAL "Release"))
    target_link_options(inc_pigment PUBLIC -static -static-libgcc -static-libstdc++)
endif()

target_sources(inc_pigment
    PUBLIC
        FILE_SET pub_headers
        TYPE HEADERS
        BASE_DIRS
                include)

target_sources(inc_pigment
    PRIVATE
        FILE_SET priv_headers
        TYPE HEADERS
        BASE_DIRS
                src/private_inc)