cmake_minimum_required(VERSION 3.25)
project(corvus_consumer_example LANGUAGES CXX)

# No CMAKE_CXX_STANDARD here on purpose: the C++20 requirement must arrive
# through corvus::corvus's exported cxx_std_20 (that export is what this
# project exists to test).
# Find the installed corvus package
find_package(corvus REQUIRED)

add_executable(consumer_demo main.cpp)
target_link_libraries(consumer_demo PRIVATE corvus::corvus)
