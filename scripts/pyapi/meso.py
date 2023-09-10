import os
from pyapi import cpp_func, link_library
from pyapi.cpp_types import *

api = link_library("api")


def init():
    @cpp_func(api.api_init)
    def _init(func):
        func()


def set_core_params(debug_mode=False, save_interval=1000,
                    residual_interval=1000, max_step=1000000,
                    residual_limit=1e-6):
    @cpp_func(api.api_set_core_params, [c_bool, c_int, c_int, c_int, c_double])
    def _set(func):
        func(debug_mode, save_interval, residual_interval, max_step, residual_limit)


def get_config_reader(case_name):
    @cpp_func(api.api_new_config_reader, [c_char_p], c_void_p)
    def reader_pt(func):
        return func(case_name.encode())
    return reader_pt


def get_config_value(reader_pt, key):
    @cpp_func(api.api_config_get_value, [c_void_p, c_char_p], c_char_p)
    def _get_value(func):
        return func(reader_pt, key.encode()).decode()
    return _get_value


class Solver:

    def __init__(self, case_name):
        self.config = get_config_reader(case_name)
        solver_type = get_config_value(self.config, "SOLVER")
        if solver_type == "dugks@incompressible":
            init_func = api._Z27api_solver_init_from_configI20DUGKS_INCOMPRESSIBLEEPT_P12ConfigReader
            self.get_step_func = api._Z19api_solver_get_stepI20DUGKS_INCOMPRESSIBLEEiPT_
            self.get_continue_func = api._Z30api_solver_get_continue_to_runI20DUGKS_INCOMPRESSIBLEEbPT_
            self.do_step_func = api._Z18api_solver_do_stepI20DUGKS_INCOMPRESSIBLEEbPT_
        elif solver_type == "dugks@aoki":
            init_func = api._Z27api_solver_init_from_configI10DUGKS_AOKIEPT_P12ConfigReader
            self.get_step_func = api._Z19api_solver_get_stepI10DUGKS_AOKIEiPT_
            self.get_continue_func = api._Z30api_solver_get_continue_to_runI10DUGKS_AOKIEbPT_
            self.do_step_func = api._Z18api_solver_do_stepI10DUGKS_AOKIEbPT_
        else:
            raise RuntimeError("Python script caught unsupported solver type<%s>." % solver_type)

        @cpp_func(init_func, [c_void_p], c_void_p)
        def solver_pt(func):
            return func(self.config)

        self.handle = solver_pt

    def get_step(self):
        @cpp_func(self.get_step_func, [c_void_p], c_int)
        def step(func):
            return func(self.handle)
        return step

    def get_continue(self):
        @cpp_func(self.get_continue_func, [c_void_p], c_bool)
        def continue_to_run(func):
            return func(self.handle)
        return continue_to_run

    def do_step(self):
        @cpp_func(self.do_step_func, [c_void_p], c_bool)
        def solver_do_step(func):
            return func(self.handle)
        return solver_do_step
