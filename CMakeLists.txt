cmake_minimum_required(VERSION 3.2)

project(seminar-competition C)

set(CMAKE_CC_FLAGS "${CMAKE_CC_FLAGS} -pg -lm")
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")

file(GLOB ARGORITHMS "src/algorithm/*.c")
set(SRCS
    src/algorithm.c
    src/cpu_time.c
    src/tsp.c
    src/utils.c
    ${ARGORITHMS})

add_executable(${PROJECT_NAME} ${SRCS})
target_link_libraries(${PROJECT_NAME} PRIVATE m)