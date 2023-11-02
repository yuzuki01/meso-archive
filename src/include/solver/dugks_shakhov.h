#ifndef SOLVER_DUGKS_SHAKHOV
#define SOLVER_DUGKS_SHAKHOV

#endif //SOLVER_DUGKS_SHAKHOV
/**
 * Compressible DUGKS
 *      BGK-Shakhov model
 *      Ventaka limiter
 */

class DUGKS_SHAKHOV {
public:
    const std::string prefix = "dugks@shakhov";
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
    bool limiter;
    int step;
    double dt, half_dt;
    /// 数学物理参数
    int K, D;
    double Kn, Ma, Re, Pr;
    double rho, T;
    double R, L;
    double miu0, vhs_index, gamma, Cv;
    double limiter_k;
    /// 算法函数
    inline double g_maxwell(double density, double temperature, double cc) const;
    inline double g_shakhov(double density, double temperature, double cc, double cq, double gm) const;
    inline double h_maxwell(double temperature, double gm) const;
    inline double h_shakhov(double density, double temperature, double cc, double cq, double gm) const;
    inline double tau_f(double density, double temperature) const;
    inline double venkata_limiter(double wi, double wi_max, double wi_min, double delta_i, double volume) const;
    /// 内部类
    using Scheme = DUGKS_SHAKHOV;
    class Cell {
    public:
        Scheme *solver;
        int id;
        Vec(Residual) residual;
        /// 最小二乘法
        LeastSecondParam lsp;
        /// 宏观量
        double density, temperature;
        Vector velocity, heat_flux;
        /// 分布函数
        // VecDouble g, h;
        VecDouble g_t, g_bp, Wg_max, Wg_min;
        Vec(Vector) slope_g;
        VecDouble h_t, h_bp, Wh_max, Wh_min;
        Vec(Vector) slope_h;
        /// 函数
        Cell(int cell_id, Scheme *Solver);

        void update_f_bp();
        void update_gradient();
        void update_macro_var();
        void update_f_t();

        std::string info() const;
    };

    class Interface {
    public:
        Scheme *solver;
        int id;
        /// 宏观量
        double density, temperature;
        Vector velocity, heat_flux;
        /// 分布函数
        VecDouble g, g_b;
        VecDouble h, h_b;
        /// 函数
        Interface(int face_id, Scheme *Solver);

        void update_f_b();
        void update_macro_var();
        void update_f();
        void boundary_condition();
    };
    /// 容器
    Vec(Cell) CELLS;
    Vec(Interface) FACES;

    /// 函数
    explicit DUGKS_SHAKHOV(ConfigReader &reader);

    void init();
    void info();
    void do_step();
    void do_save();
    void do_residual();
    void do_crashed(Scheme::Cell &cell);
};
