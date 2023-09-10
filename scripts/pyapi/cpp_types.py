from ctypes import POINTER
from ctypes import c_bool, c_int, c_double, c_void_p
from ctypes import c_char, c_char_p

c_void = c_void_p
c_bool_p = POINTER(c_bool)
c_int_p = POINTER(c_int)
c_double_p = POINTER(c_double)
