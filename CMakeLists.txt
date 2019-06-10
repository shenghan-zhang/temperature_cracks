cmake_minimum_required(VERSION 3.0.0)

project(thermal-crack)
enable_language(CXX)


find_package(Akantu REQUIRED)

set(CMAKE_CXX_STANDARD 14)

add_akantu_simulation(my_simu my_simu.cc)
