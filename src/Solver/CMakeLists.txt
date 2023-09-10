# 这是 Solver 对应的 CMake 文件
cmake_minimum_required(VERSION 3.17)
project(solver)
set(CMAKE_BUILD_TYPE MinSizeRel)

# 判断操作系统类型
if(WIN32)
    # 如果是 Windows
    message("当前操作系统为 Windows")
    set(SHARED_LIBRARY_SUFFIX ".dll")
elseif(UNIX AND NOT APPLE)
    # 如果是 Linux
    message("当前操作系统为 Linux")
    set(SHARED_LIBRARY_SUFFIX ".so")
else()
    message(FATAL_ERROR "不支持的操作系统类型")
endif()

file(GLOB SRC
        ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/dugks/*.cpp)
include_directories(${PROJECT_SOURCE_DIR}/../include)

set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/../../build) # 指定输出文件目录为 build
link_directories(${LIBRARY_OUTPUT_PATH})  # 指定 build 为 dll 文件目录

set(CMAKE_SHARED_LIBRARY_PREFIX "")                 # 去除前缀
add_library(${PROJECT_NAME} SHARED ${SRC})          # 生成动态库
target_link_libraries(${PROJECT_NAME} ${LIBRARY_OUTPUT_PATH}/core${SHARED_LIBRARY_SUFFIX})  # 链接 Core 的 dll/so 文件
target_link_libraries(${PROJECT_NAME} ${LIBRARY_OUTPUT_PATH}/mesh${SHARED_LIBRARY_SUFFIX})  # 链接 Mesh 的 dll/so 文件

# check openMp
FIND_PACKAGE(OpenMP REQUIRED)
if (OPENMP_FOUND)
    message("OpenMP Found.")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()
