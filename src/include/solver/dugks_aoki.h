#ifndef SOLVER_DUGKS_AOKI
#define SOLVER_DUGKS_AOKI
#endif
/**
 * 用于求解参考文献
 * [1] K.Aoki, S.Takata, H.Aikawa, et al. A Rarefied Gas Flow Caused by a Discontinuous Wall Temperature[J],
 *     Physics of Fluids, 2001, 13: 2645–2661.
 *      配置文件 Solver = dugks@aoki
 **/

class DUGKS_AOKI {
public:
    const std::string prefix = "dugks@aoki";
    /// 网格
    MESH::Mesh mesh;
    MESH::Mesh DVS;
    int mesh_cell_num;
    int mesh_face_num;
    double mesh_total_volume;
    /// 变量
    std::string case_name;
    bool continue_to_run;
    bool is_crashed;
    int step, D;
    double dt, half_dt;
    /// 数学、物理参数
    double Kn;
    double rho, R, T;
    double L;
    double Ac;
    /// 算法函数
    inline double tau_f(double density) const;
    inline double g_eq(double density, double temperature,
                       const Vector &particle_velocity,
                       const Vector &flow_velocity) const;
    inline double h_eq(double geq, double temperature) const;
    /// 算法内部类
    class Cell {
    public:
        DUGKS_AOKI *solver;
        int id;
        Vec(Residual) residual;
        /// 最小二乘法
        LeastSecondParam lsp;
        /// 宏观量
        double density, temperature;
        Vector velocity;
        /// 分布函数
        VecDouble g_t, g_bp;
        VecDouble h_t, h_bp;
        Vec(Vector) slope_g, slope_h;

        /// 函数
        Cell(int cell_id, DUGKS_AOKI *Solver);

        void update_f_bp();

        void update_gradient();

        void update_macro_var();

        void update_f_t();

        std::string info() const;
    };

    class Interface {
    public:
        DUGKS_AOKI *solver;
        int id;
        /// 宏观量
        double density, temperature;
        Vector velocity;
        /// 分布函数
        VecDouble g, g_b;
        VecDouble h, h_b;

        /// 函数
        Interface(int face_id, DUGKS_AOKI *Solver);

        void update_f_b();

        void update_macro_var();

        void update_f();

        void boundary_condition();
    };
    /// 容器
    Vec(Cell) CELLS, SHADOW_CELLS;
    Vec(Interface) FACES;

    /// 函数
    explicit DUGKS_AOKI(ConfigReader &reader);

    void init();

    void info();

    void do_step();

    void do_save();

    void do_residual();

    void do_crashed(DUGKS_AOKI::Cell &cell);
};
