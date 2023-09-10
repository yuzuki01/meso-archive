import os
import ctypes
from pyapi.cpp_types import *


def link_library(file_name):
    if os.name == "nt":
        return ctypes.CDLL("%s.dll" % file_name, winmode=0)
    elif os.name == "posix":
        return ctypes.CDLL("./%s.so" % file_name)
    else:
        raise OSError("link_library() caught os<name=%s>" % os.name)


# cpp 函数调用装饰器
def cpp_func(dll_func, arg_types=None, res_type=c_void):
    if arg_types is None:
        arg_types = []
    dll_func.argtypes = arg_types
    dll_func.restype = res_type

    def decorator(func):
        return func(dll_func)

    return decorator
