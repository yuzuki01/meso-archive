/**************************************************
 *  Well-Balanced DUGKS                           *
 *    with phase-field model                      *
 *    for incompressible multi-phase case         *
 * ---------------------------------------------- *
 *                 Created by MYC on Nov 1, 2023  *
 **************************************************/

#include <Solver.h>
#ifdef SOLVER_WBDUGKS_PF


using Scheme = WBDUGKS_PF;
using SCell = Scheme::Cell;
using SFace = Scheme::Interface;


/// 构造函数
SCell::Cell(int cell_id, Scheme *Solver) {
    id = cell_id;
    solver = Solver;
    const int dvs_num = solver->DVS.NELEM;
    order_param = density = pressure = 0.0;
    velocity = {0.0, 0.0, 0.0};
    // 最小二乘法
    lsp = GenerateLeastSecondParam(id, solver->mesh);
    // 分布函数
    f_t.resize(dvs_num, 0.0);
    f_bp.resize(dvs_num, 0.0);
    g_t.resize(dvs_num, 0.0);
    g_bp.resize(dvs_num, 0.0);
    slope_f.resize(dvs_num, {0.0, 0.0, 0.0});
    slope_g.resize(dvs_num, {0.0, 0.0, 0.0});
    source_f.resize(dvs_num, 0.0);
    source_g.resize(dvs_num, 0.0);
    // shrink to fit
    f_t.shrink_to_fit();
    f_bp.shrink_to_fit();
    g_t.shrink_to_fit();
    g_bp.shrink_to_fit();
    slope_f.shrink_to_fit();
    slope_g.shrink_to_fit();
    source_f.shrink_to_fit();
    source_g.shrink_to_fit();
}

SFace::Interface(int face_id, Scheme *Solver) {
    id = face_id;
    solver = Solver;
    const int dvs_num = solver->DVS.NELEM;
    order_param = density = pressure = 0.0;
    f.resize(dvs_num, 0.0);
    f_b.resize(dvs_num, 0.0);
    g.resize(dvs_num, 0.0);
    g_b.resize(dvs_num, 0.0);
    // shrink to fit
    f.shrink_to_fit();
    f_b.shrink_to_fit();
    g.shrink_to_fit();
    g_b.shrink_to_fit();
}

/// 算法初始化
void Scheme::init() {
    int output_count = 0;
    /// 创建对应算法控制体
    for (auto &cell : mesh.CELLS) {
        CELLS.emplace_back(cell.id, this);
        auto &scheme_cell = CELLS.back();
        // 赋值
        for (int i = 0; i < DVS.NELEM; i++) {
            auto &particle = DVS.CELLS[i];

            /// scheme_cell.f_t[i] =
            /// scheme_cell.g_t[i] =
        }
        scheme_cell.update_macro_var();
        if (debug_mode && ++output_count <= 10) {
            pprint::debug << "Generate cell: id=" << scheme_cell.id << " order_param=" << scheme_cell.order_param
                          << " density=" << scheme_cell.density << " pressure=" << scheme_cell.pressure;
            pprint::debug(prefix);
        } else if (debug_mode && ++output_count == 12) {
            pprint::debug << "  ......(10 cells info shown)";
            pprint::debug();
        }
        // 残差
        mesh_total_volume += cell.volume;
        scheme_cell.residual.emplace_back(scheme_cell.density);
        scheme_cell.residual.emplace_back(scheme_cell.order_param);
        scheme_cell.residual.emplace_back(scheme_cell.pressure);
        scheme_cell.residual.emplace_back(scheme_cell.velocity.x);
        scheme_cell.residual.emplace_back(scheme_cell.velocity.y);
        if (D == 3) scheme_cell.residual.emplace_back(scheme_cell.velocity.z);
        scheme_cell.residual.shrink_to_fit();
    }
    pprint::info << "Create scheme cells with mesh cells.";
    pprint::info(prefix);
    /// 创建对应算法界面
    for (auto &face : mesh.FACES) {
        FACES.emplace_back(face.id, this);
    }
    pprint::info << "Create scheme interfaces with mesh interfaces.";
    pprint::info(prefix);
    /// shrink_to_fit
    CELLS.shrink_to_fit();
    FACES.shrink_to_fit();
    /// 变量
    continue_to_run = true;
    pprint::note << "Initialization finished.";
    pprint::note(prefix);
}

/// 平衡态函数
double Scheme::f_eq(double phi, const MESH::Cell &particle, const Vector &flow_velocity) const {
    return particle.volume * phi * (1.0 + (particle.pos * flow_velocity) / RT);
}

double Scheme::g_eq(double density, double pressure, const MESH::Cell &particle, const Vector &flow_velocity) const {
    /// if (particle.)
    return 0.0; ///debug
}

#endif // SOLVER_WBDUGKS_PF
