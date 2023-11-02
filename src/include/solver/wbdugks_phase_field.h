/**************************************************
 *  Well-Balanced DUGKS                           *
 *    with phase-field model                      *
 *    for incompressible multi-phase case         *
 * ---------------------------------------------- *
 *                 Created by MYC on Nov 1, 2023  *
 **************************************************/

#ifndef SOLVER_WBDUGKS_PF
#define SOLVER_WBDUGKS_PF

#endif //SOLVER_WBDUGKS_PF

class WBDUGKS_PF {
public:
    const std::string prefix = "wbdugks@phase-field";
    /// 网格
    MESH::Mesh mesh;            // 物理网格
    MESH::Mesh DVS;             // 速度网格
    int mesh_cell_num;          // 网格单元数
    int mesh_face_num;          // 网格界面数
    double mesh_total_volume;   // 网格总体积
    /// 变量
    std::string case_name;
    bool continue_to_run;       // 是否运行
    bool is_crashed;            // 是否发散
    int step;                   // 迭代步数
    double dt, half_dt;         // 时间步长
    /// 数学、物理参数
    int D;
    double Pe, Cn;
    double RT;
    double sigma;               // 表面张力系数
    double interface_width;     // 界面宽度
    double kappa, beat;
    double density_gas, density_fluid;
    double phi_gas, phi_fluid;  // 序参数
    /// 算法函数
    inline double f_eq(double phi,
                       const MESH::Cell &particle,
                       const Vector &flow_velocity) const;
    inline double g_eq(double density, double pressure,
                       const MESH::Cell &particle,
                       const Vector &flow_velocity) const;

    using Scheme = WBDUGKS_PF;
    /// 算法内部类
    class Cell {
    public:
        Scheme *solver;
        int id;
        Vec(Residual) residual;
        /// 最小二乘法
        LeastSecondParam lsp;
        /// 宏观量
        double order_param, density, pressure;
        Vector velocity;
        /// 分布函数
        VecDouble f_t, f_bp, g_t, g_bp;
        VecDouble source_f, source_g;
        Vec(Vector) slope_f, slope_g;

        /// 函数
        Cell(int cell_id, Scheme *Solver);

        void update_f_bp();
        void update_gradient();
        void update_macro_var();
        void update_f_t();
        void update_source();

        std::string info() const;
    };

    class Interface {
    public:
        Scheme *solver;
        int id;
        /// 宏观量
        double order_param, density, pressure;
        Vector velocity;
        /// 分布函数
        VecDouble f, f_b, g, g_b;

        /// 函数
        Interface(int face_id, Scheme *Solver);
        void update_f_b();
        void update_macro_var();
        void update_f();
        void boundary_condition();
    };
    /// 容器
    Vec(Cell) CELLS, SHADOW_CELLS;
    Vec(Interface) FACES;

    /// 函数
    explicit WBDUGKS_PF(ConfigReader &reader);

    void init();
    void info();
    void do_step();
    void do_save();
    void do_residual();
    void do_crashed(Scheme::Cell &cell);
};
