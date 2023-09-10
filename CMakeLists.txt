cmake_minimum_required(VERSION 3.17)

project(meso)

add_custom_target(Core)

set(PROJECT_VERSION 1.1.0)
set(CMAKE_CXX_STANDARD 14)
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

# 添加头文件目录
include_directories(./src/include)

# 设置可执行文件输出路径
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/build)

add_executable(${PROJECT_NAME}
        main.cpp
        main.h
        )

# 动态链接库编译时设置
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/build) # 指定输出文件目录为 build
link_directories(${LIBRARY_OUTPUT_PATH})  # 指定 build 为 dll 文件目录
target_link_libraries(${PROJECT_NAME} ${LIBRARY_OUTPUT_PATH}/core${SHARED_LIBRARY_SUFFIX})  # 链接 Core 的 dll/so 文件
target_link_libraries(${PROJECT_NAME} ${LIBRARY_OUTPUT_PATH}/mesh${SHARED_LIBRARY_SUFFIX})  # 链接 Mesh 的 dll/so 文件
target_link_libraries(${PROJECT_NAME} ${LIBRARY_OUTPUT_PATH}/solver${SHARED_LIBRARY_SUFFIX})  # 链接 Solver 的 dll/so 文件

# check openMp
FIND_PACKAGE(OpenMP REQUIRED)
if (OPENMP_FOUND)
    message("OpenMP Found.")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()
