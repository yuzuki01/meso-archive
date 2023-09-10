#include <Core.h>
#include <Mesh.h>
#include <Solver.h>

#if defined(_WIN32)
/// Windows
#define EXPORT_FUNC extern "C"
#define EXPORT_CLASS __declspec(dllexport)

#elif defined(__linux__)
/// Linux
#define EXPORT_FUNC extern "C"
#define EXPORT_CLASS __attribute__((visibility("default")))

#endif

/// MesoKinG 初始化
EXPORT_FUNC void api_init();

/// 可执行文件功能
// create case
EXPORT_FUNC void api_create_case(const char *name);
// set core params
EXPORT_FUNC void api_set_core_params(bool api_debug_mode, int api_save_interval,
                                     int api_residual_interval, int api_max_step,
                                     double api_residual_limit);
// parse mesh
EXPORT_FUNC int api_parse_mesh(const char *path);
// new config-reader handle
EXPORT_FUNC ConfigReader *api_new_config_reader(const char *path);
EXPORT_FUNC const char *api_config_get_value(ConfigReader *reader, const char *key);


/// 模板函数
// 从配置文件初始化求解器
template <class Solver>
Solver *api_solver_init_from_config(ConfigReader *reader);
// get step
template <class Solver>
int api_solver_get_step(Solver *solver);
// get continue_to_run
template <class Solver>
bool api_solver_get_continue_to_run(Solver *solver);
// do step
template <class Solver>
bool api_solver_do_step(Solver *solver);
