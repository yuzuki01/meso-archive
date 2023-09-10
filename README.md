# Meso

## 已实现
<hr>
<li>启动参数解析</li>
<li>配置文件读取</li>
<li>支持数学运算的空间向量 Vector 类</li>
<li>2D / 3D 网格读取</li>
<li>2D / 3D 网格界面构建</li>
<li>2D / 3D 网格几何参数计算</li>
<li>2D / 3D 网格法向量构建</li>
<li>支持数组运算的分布函数基类</li>
<li>网格梯度通用型求解函数</li>
<li>不可压缩 DUGKS 求解器 [dugks@incompressible] (D2Q9 & D3Q27)</li>
<li>残差管理</li>
<li>动态链接库接口 (供 Python ctypes 调用)</li>
<li>Thermal Creep 算例复现的求解器 [dugks@aoki] dugks@shakhov 的 2 维特殊情况</li>
<hr>

## 待实现

<hr>
<li>dugks@shakhov 可压缩求解器</li>
<hr>

## Quick Start

### 构建项目

```python build.py```
<p>该脚本能构建出所有动态链接库和可执行文件</p>

```python build_api.py```
<p>该脚本能构建出 api 动态链接库，可以使用 Python ctypes 调用</p>

```python pack.py```
<p>该脚本能够将除构建生成的文件(./build 下所有文件)以外的其他文件，包括自身打包为 .zip 文件。储存在 ./pack 下</p>

### 方法一：可执行程序

<p>进入可执行程序 MesoKinG 所在目录：</p>

```meso <command>```

<p>使用如下指令获取帮助：</p>

```meso -help```

```
Mesoscopic Kinetic Solver
  meso [<command>]
    -h / -help                      Get help info
    --case <NAME>                   Run case <NAME>
    --create <NAME>                 Create a case named <NAME> with 
                                    a default config file
    --max_step <VALUE>              Set max-step to <VALUE>,
                                    default = 1000000
    --save_interval <VALUE>         Set save-interval to <VALUE>,
                                    default = 1000
    --residual_interval <VALUE>     Set residual-interval to <VALUE>,
                                    default = 1000
    --esp <VALUE>                   Set residual-limit to <VALUE>,
                                    default = 1e-6
    --parse_mesh <FILE>             Parse mesh file <FILE>, output
                                    tecplot file as parsed_mesh.dat
```

### 方法二：脚本

   <p>进入动态链接库 api.dll/api.so 所在目录：</p>

```python meso.py```

   <p>这种方法必须要构建出 api.dll/api.so ，通过 Python 调用 cpp 动态链接库内的函数，脚本文件位于 ./scripts 下，可根据需求修改示例脚本文件</p>
   <p>可调用的函数都在 ./src/include/API.h 中给出</p>
