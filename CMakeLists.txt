cmake_minimum_required(VERSION 3.10)
project(KMCThinFilm)

file(READ "VERSION" KMC_THIN_FILM_VERSION_STRING)
include(KMCThinFilmVersion.cmake)

add_subdirectory(src)
add_subdirectory(doc)
