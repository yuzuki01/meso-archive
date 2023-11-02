# Meso

基于气体动理学的 CFD 求解器，主要是 Discrete Unified Gas Kinetic Scheme (DUGKS) 的实现。

## 已实现

| 求解器 | 碰撞模型 | 其他特性 | 备注 |
| :----: | :----: | :----: | :----: |
| dugks@incompressible | BGK |-| Stable |
| dugks@shakhov | BGK-Shakhov | Limiter | UNSTABLE for Ma > 1 |
| dugks@aoki | BGK | 等价于 BGK-Shakhov (K=0, Pr=1) | Ray-effect |
| wbdugks | BGK-Shakhov | Limiter | Sharing with dugks@shakhov |

## 未来计划

### Nov 2, 2023
 - [开发] wbdugks@phase_field 多相流求解器
 - [重要] 未来可能放弃 Windows 平台的兼容


## Quick Start

### 环境配置

|  | version |
|:---:|:---:|
| gcc | \>=13.2.x |
| cmake | \>=3.25.x |
| Python | \>=3.9.x |

### 构建项目

```python build.py```
<p>该脚本能构建出所有动态链接库和可执行文件，在 Windows 10 和 Ubuntu 22.04 下都构建成功过。</p>

```python build_api.py```
<p>该脚本能构建出 api 动态链接库，可以使用 Python ctypes 调用</p>

```python pack.py```
<p>该脚本能够将除构建生成的文件(./build 下所有文件)以外的其他文件，包括自身打包为 .zip 文件。储存在 ./pack 下</p>

### 方法一：可执行程序

<p>进入可执行程序 meso 所在目录：</p>

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
