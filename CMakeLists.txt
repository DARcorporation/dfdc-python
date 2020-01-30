cmake_minimum_required(VERSION 3.12)
project(dfdc Fortran)

file(GLOB all_src "src/*.f90")

enable_language(Fortran)
add_library(dfdc SHARED ${all_src})