#ifndef INCLUDE_CORE
#include <Core.h>
#endif

#ifndef INCLUDE_MESH
#include <Mesh.h>
#endif

/// 数值方法
#include "utils/Numerical.h"


/**
 * 内置求解器
 * 从此次添加求解器
 **/

#include "solver/dugks_incompressible.h"
#include "solver/dugks_aoki.h"
#include "solver/wbdugks_shakhov.h"

/// 函数
MESH::Mesh GenerateMeshFromConfig(ConfigReader &reader, int mesh_type);
MESH::Mesh GenerateMeshGH(int dimension, double RT);

/// 可执行程序调用入口函数
template <class Solver>
int handle_solver(ConfigReader &reader);

/// 模板函数
#include "tpps/solver.tpp"
