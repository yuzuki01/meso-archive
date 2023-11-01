#ifndef SOLVER_DUGKS_INCOMPRESSIBLE
#define SOLVER_DUGKS_INCOMPRESSIBLE
#endif
/**
 * 不可压缩 DUGKS
 *      配置文件 Solver = dugks@incompressible
 */

class DUGKS_INCOMPRESSIBLE {
public:
    const std::string prefix = "dugks@incompressible";
    /// 网格
    MESH::Mesh mesh;            // 物理网格
    MESH::Mesh DVS;             // 速度网格
    int mesh_cell_num;          // 网格单元数
    int mesh_face_num;          // 网格界面数
    double mesh_total_volume;   // 网格总体积
    /// 变量
    std::string case_name;
    bool continue_to_run;   // 是否运行
    bool is_crashed;        // 是否发散
    int step;               // 迭代步数
    double dt, half_dt;     // 时间步长
    /// 数学、物理参数
    double Re, Ma;
    double rho, RT;
    double L;               // 特征长度
    double tau;             // 松弛时间
    /// 算法函数
    inline double f_eq(double density,
                       const Vector &particle_velocity,
                       const Vector &flow_velocity) const;

    /// 算法内部类
    class Cell {
    public:
        DUGKS_INCOMPRESSIBLE *solver;
        int id;
        Vec(Residual) residual;
        /// 最小二乘法
        LeastSecondParam lsp;
        /// 宏观量
        double density, temperature;
        Vector velocity, heat_flux;
        /// 分布函数
        VecDouble f_t, f_bp;
        Vec(Vector) slope_f;

        /// 函数
        Cell(int cell_id, DUGKS_INCOMPRESSIBLE *Solver);

        void update_f_bp();

        void update_gradient();

        void update_macro_var();

        void update_f_t();

        std::string info() const;
    };

    class Interface {
    public:
        DUGKS_INCOMPRESSIBLE *solver;
        int id;
        /// 宏观量
        double density, temperature;
        Vector velocity;
        /// 分布函数
        VecDouble f, f_b;

        /// 函数
        Interface(int face_id, DUGKS_INCOMPRESSIBLE *Solver);

        void update_f_b();

        void update_macro_var();

        void update_f();

        void boundary_condition();
    };
    /// 容器
    Vec(Cell) CELLS;
    Vec(Interface) FACES;

    /// 函数
    explicit DUGKS_INCOMPRESSIBLE(ConfigReader &reader);

    void init();

    void info();

    void do_step();

    void do_save();

    void do_residual();

    void do_crashed(DUGKS_INCOMPRESSIBLE::Cell &cell);
};
