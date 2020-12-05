cmake_minimum_required(VERSION 3.2)

project(seminar-competition C)

set(CMAKE_CC_FLAGS "${CMAKE_CC_FLAGS} -lm")

set(SRCS
    src/algorithm.c
    src/cpu_time.c
    src/tsp.c
    src/utils.c    )

add_executable(${PROJECT_NAME} ${SRCS})
target_link_libraries(${PROJECT_NAME} PRIVATE m)