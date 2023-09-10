#include <API.h>

/// init from config
template<class Solver>
Solver *api_solver_init_from_config(ConfigReader *reader) {
    auto solver = new Solver(*reader);
    try {
        solver->init();
        solver->info();
    } catch (std::exception &ex) {
        pprint::error << "api_solver_init_from_config() caught error: " << ex.what();
        pprint::error(solver->prefix);
        throw std::invalid_argument("api_solver_init_from_config() caught error.");
    }
    return solver;
}
// 显式实例化
template DUGKS_INCOMPRESSIBLE *api_solver_init_from_config(ConfigReader *reader);
template DUGKS_AOKI *api_solver_init_from_config(ConfigReader *reader);

/// get step
template <class Solver>
int api_solver_get_step(Solver *solver) {
    return solver->step;
}
// 显式实例化
template int api_solver_get_step<DUGKS_INCOMPRESSIBLE>(DUGKS_INCOMPRESSIBLE *solver);
template int api_solver_get_step <DUGKS_AOKI>(DUGKS_AOKI *solver);

/// get continue_to_run
template <class Solver>
bool api_solver_get_continue_to_run(Solver *solver) {
    return solver->continue_to_run;
}
// 显式实例化
template bool api_solver_get_continue_to_run<DUGKS_INCOMPRESSIBLE>(DUGKS_INCOMPRESSIBLE *solver);
template bool api_solver_get_continue_to_run<DUGKS_AOKI>(DUGKS_AOKI *solver);

/// do step
template <class Solver>
bool api_solver_do_step(Solver *solver) {
    solver->do_step();
    return solver->continue_to_run;
}
// 显式实例化
template bool api_solver_do_step<DUGKS_INCOMPRESSIBLE>(DUGKS_INCOMPRESSIBLE *solver);
template bool api_solver_do_step<DUGKS_AOKI>(DUGKS_AOKI *solver);
