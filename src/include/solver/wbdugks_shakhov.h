#ifndef INCLUDE_WBDUGKS_SHAKHOV
#define INCLUDE_WBDUGKS_SHAKHOV
#endif //INCLUDE_WBDUGKS_SHAKHOV

/**
 * Well-Balance DUGKS
 *  with BGK-Shakhov model
 **/

class WBDUGKS_SHAKHOV {
public:
    const std::string prefix = "wbdugks@shakhov";
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
    int step;
    double dt, half_dt;
    /// 数学、物理参数
    int D, K;
    double Re, Ma, Kn, Fr, Pr;
    double rho, R, L, T, gamma, Cv;
    double miu0, vhs_index;
    Vector gravity;

    /// 算法函数
    inline double g_eq(double density, double temperature,
                       double cc) const;

    inline double h_eq(double temperature, double geq) const;

    inline double g_pr(double density, double temperature,
                       double cc, double cq, double geq) const;

    inline double h_pr(double density, double temperature,
                       double cc, double cq, double geq) const;

    inline double tau_f(double density, double temperature) const;

    /// 算法内部类
    class Cell {
    public:
        WBDUGKS_SHAKHOV *solver;
        int id;
        Vec(Residual) residual;
        /// 最小二乘法
        LeastSecondParam lsp;
        /// 宏观量
        double density, temperature;
        Vector velocity, heat_flux;
        /// 分布函数
        VecDouble g_t, g_bp, h_t, h_bp;
        Vec(Vector) slope_g, slope_h;

        /// 函数
        Cell(int cell_id, WBDUGKS_SHAKHOV *Solver);

        void update_f_bp();

        void update_gradient();

        void update_macro_var();

        void update_f_t();

        std::string info() const;
    };

    class Interface {
    public:
        WBDUGKS_SHAKHOV *solver;
        int id;
        /// 宏观量
        double density, temperature;
        Vector velocity, heat_flux;
        /// 分布函数
        VecDouble g, g_b, h, h_b;

        /// 函数
        Interface(int face_id, WBDUGKS_SHAKHOV *Solver);

        void update_f_b();

        void update_macro_var();

        void update_f();

        void boundary_condition();
    };

    /// 容器
    Vec(Cell) CELLS, SHADOW_CELLS;
    Vec(Interface) FACES;

    /// 函数
    explicit WBDUGKS_SHAKHOV(ConfigReader &reader);

    void init();

    void info();

    void do_step();

    void do_save();

    void do_residual();

    void do_crashed(WBDUGKS_SHAKHOV::Cell &cell);
};
