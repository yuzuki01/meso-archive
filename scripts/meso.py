from pyapi.meso import *


init()
set_core_params(max_step=2000)

solver = Solver("demo")
flag = True
while flag:
    flag = solver.do_step()
