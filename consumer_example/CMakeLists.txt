cmake_minimum_required(VERSION 3.25)
project(corvus_consumer_example LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find the installed corvus package
find_package(corvus REQUIRED)

add_executable(consumer_demo main.cpp)
target_link_libraries(consumer_demo PRIVATE corvus::corvus)
