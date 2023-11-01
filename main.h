#include <Core.h>
#include <Mesh.h>
#include <Solver.h>


int handle_parse_failed();
int handle_case(const std::string &path);
int handle_mesh(const std::string &path);
int handle_create_case(const std::string &name);
int handle_help();
int handle_debug_func(const std::string &_string);

void handle_debug();
