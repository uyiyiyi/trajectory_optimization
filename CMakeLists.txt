#cmake最小版本需求
cmake_minimum_required(VERSION 2.8.3)
set (CMAKE_CXX_STANDARD 17)
#指定编译类型，debug 或者为 release
# debug 会生成相关调试信息，可以使用 GDB 进行
# release 不会生成调试信息。当无法进行调试时查看此处是否设置为 debug.
set(CMAKE_BUILD_TYPE release)
#设置此项目的名称
project(trajectory_optimization) 

find_package(Eigen3 REQUIRED)
find_package(osqp REQUIRED)


# 设置一个名字var_name 的变量，同时给此变量赋值为var_value
# SET(var_name var_value)


# 打印消息
# MESSAGE("MSG") 

#给变量var_name赋值为var_value，comment是此变量的注释，和SET 有类似的功效，用于给某变量设置默认值
# option(var_name "comment" var_value) 

# 添加include路径，也就是头文件路径
include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${EIGEN3_INCLUDE_DIRS}
) 

#生成可执行文件target ，后面填写的是生成此可执行文件所依赖的源文件列表。
add_executable(hello src/hello.cpp)
add_executable(min_snap src/min_snap.cpp)

# Link the OSQP shared library
target_link_libraries(min_snap PRIVATE osqp::osqp)

# or...
# Link the OSQP static library
# target_link_libraries(yourTarget PRIVATE osqp::osqpstatic)



# 调用xxx子目录的CMakeLists.txt执行
# add_subdirectory(xxx) 

# 给编译器添加xxx参数
# add_compile_options(xxx)

# 给编译器添加库目录，
# link_directories(xxx)

# 生成库文件，SHARED代表动态库，STATIC代表静态库， 最后一个参数代表此库的源文件列表
# add_library(lib_name SHARED or STATIC lib_source_code) 

# 给目标添加依赖库
# target_link_libraries(test)