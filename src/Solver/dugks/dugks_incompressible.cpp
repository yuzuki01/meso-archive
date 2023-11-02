/**
 * 不压缩 DUGKS
 */

#include <Solver.h>
#ifdef SOLVER_DUGKS_INCOMPRESSIBLE

using Scheme = DUGKS_INCOMPRESSIBLE;
using SCell = Scheme::Cell;
using SFace = Scheme::Interface;

/// 构造函数
SCell::Cell(int cell_id, Scheme *Solver) {
    id = cell_id;
    solver = Solver;
    int dvs_num = solver->DVS.NELEM;
    density = 0.0;
    // 最小二乘法
    lsp = GenerateLeastSecondParam(id, solver->mesh);
    // 分布函数
    f_t.resize(dvs_num, 0.0);
    f_bp.resize(dvs_num, 0.0);
    slope_f.resize(dvs_num, Vector(0.0, 0.0, 0.0));
    // shrink to fit
    f_t.shrink_to_fit();
    f_bp.shrink_to_fit();
    slope_f.shrink_to_fit();
}

SFace::Interface(int face_id, DUGKS_INCOMPRESSIBLE *Solver) {
    id = face_id;
    solver = Solver;
    int dvs_num = solver->DVS.NELEM;
    density = 0.0;
    // 分布函数
    f.resize(dvs_num, 0.0);
    f_b.resize(dvs_num, 0.0);
    // shrink to fit
    f.shrink_to_fit();
    f_b.shrink_to_fit();
}

/// 算法初始化
void Scheme::init() {
    /// 创建对应算法控制体
    for (auto &cell : mesh.CELLS) {
        CELLS.emplace_back(cell.id, this);
        auto &scheme_cell = CELLS.back();
        // 赋值
        for (int i = 0; i < DVS.NELEM; i++) {
            auto &particle = DVS.CELLS[i];
            scheme_cell.f_t[i] = f_eq(rho, particle.pos, Vector(0.0, 0.0, 0.0));
        }
        scheme_cell.update_macro_var();
        // 残差
        mesh_total_volume += cell.volume;
        scheme_cell.residual.emplace_back(scheme_cell.density);
        scheme_cell.residual.emplace_back(scheme_cell.velocity.x);
        scheme_cell.residual.emplace_back(scheme_cell.velocity.y);
        if (mesh.dimension() == 3) scheme_cell.residual.emplace_back(scheme_cell.velocity.z);
        scheme_cell.residual.shrink_to_fit();
    }
    /// 创建对应算法界面
    for (auto &face : mesh.FACES) {
        FACES.emplace_back(face.id, this);
    }
    /// 变量
    continue_to_run = true;
    pprint::note << "Initialization finished.";
    pprint::note(prefix);
}

/// 平衡态函数
double Scheme::f_eq(double density, const Vector &particle_velocity, const Vector &flow_velocity) const {
    double uu_RT = flow_velocity * flow_velocity / RT;
    double ku_RT = particle_velocity * flow_velocity / RT;
    return density * (1.0 + ku_RT + ku_RT * ku_RT / 2.0 - uu_RT / 2.0);
}

/// 控制体算法函数
void SCell::update_f_bp() {
    double C = 3.0 * solver->half_dt / (2.0 * solver->tau + solver->dt);
    for (int i = 0; i < solver->DVS.NELEM; i++) {
        auto &particle = solver->DVS.CELLS[i];
        f_bp[i] = (1.0 - C) * f_t[i] + C * solver->f_eq(density, particle.pos, velocity);
    }
}

void SCell::update_gradient() {
    auto &cell = solver->mesh.CELLS[id];
    int near_num = cell.near_cell_id.size();
    for (int i = 0; i < solver->DVS.NELEM; i++) {
        Vector Sfr(0.0, 0.0, 0.0);
        for (int j = 0; j < near_num; j++) {
            auto &near_cell = solver->CELLS[cell.near_cell_id[j]];
            Sfr += lsp.weight[j] * (near_cell.f_bp[i] - f_bp[i]) * lsp.dr[j];
        }
        switch (solver->mesh.dimension()) {
            case 2:
                slope_f[i] = {lsp.Cx * Sfr, lsp.Cy * Sfr, 0.0};
                break;
            case 3:
                slope_f[i] = {lsp.Cx * Sfr, lsp.Cy * Sfr, lsp.Cz * Sfr};
                break;
            default:
                throw std::invalid_argument("update_gradient() dimension error.");
        }
        /// slope_f[i] = {0.0, 0.0, 0.0};   // set zero gradient
    }
}

void SCell::update_macro_var() {
    density = 0.0;
    velocity = {0.0, 0.0, 0.0};
    for (int i = 0; i < solver->DVS.NELEM; i++) {
        auto &particle = solver->DVS.CELLS[i];
        density += particle.volume * f_t[i];
        velocity += particle.volume * f_t[i] * particle.pos;
    }
    velocity /= density;

    if (std::isnan(density)) {
        solver->do_crashed(*this);
    }
}

void SCell::update_f_t() {
    // 计算通量
    VecDouble flux(solver->DVS.NELEM, 0.0);
    // 遍历界面
    auto &mesh_cell = solver->mesh.CELLS[id];
    for (int i = 0; i < solver->DVS.NELEM; i++) {
        auto &particle = solver->DVS.CELLS[i];
        for (int face_id : mesh_cell.face_id) {
            auto &face = solver->FACES[face_id];
            auto &mesh_face = solver->mesh.FACES[face_id];
            auto &nv = (mesh_face.on_cell_id == id) ? mesh_face.on_cell_nv
                                                    : mesh_face.inv_cell_nv;
            flux[i] += (particle.pos * nv) * face.f[i] * (mesh_face.area / mesh_cell.volume * solver->dt);
        }
    }
    // 更新 f_t
    for (int i = 0; i < solver->DVS.NELEM; i++) {
        f_t[i] = (4.0 * f_bp[i] - f_t[i]) / 3.0 - flux[i];
    }
}

std::string SCell::info() const {
    std::stringstream ss;
    ss << solver->mesh.CELLS[id].info_with_pos() << " density=" << density << " velocity=" << velocity.info();
    return ss.str();
}

/// 界面算法函数
void SFace::update_f_b() {
    auto &face = solver->mesh.FACES[id];
    auto &on_cell = solver->CELLS[face.on_cell_id];
    auto &inv_cell = solver->CELLS[face.inv_cell_id];
    auto &nv = face.on_cell_nv;
    for (int i = 0; i < solver->DVS.NELEM; i++) {
        auto &particle = solver->DVS.CELLS[i];
        if (particle.pos * nv >= 0.0) {
            Vector dr = face.pos - solver->mesh.CELLS[face.on_cell_id].pos - solver->half_dt * particle.pos;
            f_b[i] = on_cell.f_bp[i] + on_cell.slope_f[i] * dr;
        } else {
            Vector dr = face.pos - solver->mesh.CELLS[face.inv_cell_id].pos - solver->half_dt * particle.pos;
            f_b[i] = inv_cell.f_bp[i] + inv_cell.slope_f[i] * dr;
        }
    }
}

void SFace::update_macro_var() {
    density = 0.0;
    velocity = {0.0, 0.0, 0.0};
    for (int i = 0; i < solver->DVS.NELEM; i++) {
        auto &particle = solver->DVS.CELLS[i];
        density += particle.volume * f_b[i];
        velocity += particle.volume * f_b[i] * particle.pos;
    }
    velocity /= density;
}

void SFace::update_f() {
    double C = solver->half_dt / (2.0 * solver->tau + solver->half_dt);
    for (int i = 0; i < solver->DVS.NELEM; i++) {
        auto &particle = solver->DVS.CELLS[i];
        f[i] = (1.0 - C) * f_b[i] + C * solver->f_eq(density, particle.pos, velocity);
    }
}

void SFace::boundary_condition() {
    /// default:
    /// "symmetry"          5
    /// "periodic"          6
    auto &face = solver->mesh.FACES[id];
    switch (face.boundary_type) {
        case 0: {
            /// "interface"         0
            return;
        }
        case 1: {
            /// "inlet"             1
            auto &mark = solver->mesh.MARKS[face.mark_id];
            auto &nv = face.inv_cell_nv;
            for (int i = 0; i < solver->DVS.NELEM; i++) {
                auto &particle = solver->DVS.CELLS[i];
                if (particle.pos * nv >= 0.0) {
                    f[i] = solver->f_eq(mark.density, particle.pos, mark.velocity);
                }
            }
        }
            break;
        case 2: {
            /// "outlet"            2
            auto &near_cell = solver->CELLS[face.on_cell_id];
            auto &nv = face.inv_cell_nv;
            for (int i = 0; i < solver->DVS.NELEM; i++) {
                auto &particle = solver->DVS.CELLS[i];
                if (particle.pos * nv >= 0.0) {
                    f[i] = solver->f_eq(near_cell.density, particle.pos, near_cell.velocity);
                }
            }
        }
            break;
        case 3:
            /// "isothermal_wall"   3
        case 4:
            /// "adiabat_wall"      4
        {
            auto &mark = solver->mesh.MARKS[face.mark_id];
            auto &nv = face.inv_cell_nv;
            double rho_w, rho_w0;
            rho_w = rho_w0 = 0.0;
            for (int i = 0; i < solver->DVS.NELEM; i++) {
                auto &particle = solver->DVS.CELLS[i];
                double kn = particle.pos * nv;
                if (kn >= 0.0) {
                    rho_w0 += particle.volume * kn * solver->f_eq(1.0, particle.pos, mark.velocity);
                } else {
                    rho_w -= particle.volume * kn * f[i];
                }
            }
            rho_w /= rho_w0;
            for (int i = 0; i < solver->DVS.NELEM; i++) {
                auto &particle = solver->DVS.CELLS[i];
                if (particle.pos * nv >= 0.0) {
                    f[i] = solver->f_eq(rho_w, particle.pos, mark.velocity);
                }
            }
        }
            break;
        default:
            pprint::error << "Caught unsupported boundary condition.";
            pprint::error(solver->prefix);
            throw std::invalid_argument("Caught unsupported boundary condition.");
    }
}

/// 算法函数
void Scheme::do_step() {
#pragma omp parallel for
    for (int j = 0; j < mesh_cell_num; j++) {
        auto &cell = CELLS[j];
        cell.update_f_bp();
    }
    //pprint::debug << "fbp";
    //pprint::debug("DEBUG");
#pragma omp parallel for
    for (int j = 0; j < mesh_cell_num; j++) {
        auto &cell = CELLS[j];
        cell.update_gradient();
    }

    //pprint::debug << "gradient";
    //pprint::debug("DEBUG");
#pragma omp parallel for
    for (int j = 0; j < mesh_face_num; j++) {
        auto &face = FACES[j];
        face.update_f_b();
    }

    //pprint::debug << "fb";
    //pprint::debug("DEBUG");
#pragma omp parallel for
    for (int j = 0; j < mesh_face_num; j++) {
        auto &face = FACES[j];
        face.update_macro_var();
    }

    //pprint::debug << "cell macro var";
    //pprint::debug("DEBUG");
#pragma omp parallel for
    for (int j = 0; j < mesh_face_num; j++) {
        auto &face = FACES[j];
        face.update_f();
    }

    //pprint::debug << "f";
    //pprint::debug("DEBUG");
#pragma omp parallel for
    for (int j = 0; j < mesh_face_num; j++) {
        auto &face = FACES[j];
        face.boundary_condition();
    }

    //pprint::debug << "bc";
    //pprint::debug("DEBUG");
#pragma omp parallel for
    for (int j = 0; j < mesh_cell_num; j++) {
        auto &cell = CELLS[j];
        cell.update_f_t();
    }

    //pprint::debug << "ft";
    //pprint::debug("DEBUG");
#pragma omp parallel for
    for (int j = 0; j < mesh_cell_num; j++) {
        auto &cell = CELLS[j];
        cell.update_macro_var();
    }

    //pprint::debug << "face macro var";
    //pprint::debug("DEBUG");
    /// OpenMP loop end
    ++step;
    if (debug_mode) {
        pprint::debug << "Run step " << step;
        pprint::debug(prefix);
    }
    /// step over from here
    // break
    if (step >= max_step) continue_to_run = false;
    if (is_crashed) continue_to_run = false;
    pprint::process_bar(step, residual_interval);
    // residual
    if (step % residual_interval == 0) {
        do_residual();
    }
    // save
    if (step % save_interval == 0) {
        do_save();
    }
}


/// 功能函数
void Scheme::do_save() {
    std::stringstream ss;
    ss << "./case/" << case_name << "/result/step_" << step << ".dat";
    std::ofstream fp;
    fp.open(ss.str(), std::ios::out | std::ios::trunc);
    // check
    if (!fp.is_open()) {
        pprint::error << "Cannot write to file: " << ss.str();
        pprint::error(prefix);
        throw std::invalid_argument("Cannot write to file.");
    }
    switch (mesh.dimension()) {
        case 2:
            fp << GEOM::tecplot_file_header(mesh, {"density", "velocity-x", "velocity-y"});
            break;
        case 3:
            fp << GEOM::tecplot_file_header(mesh, {"density", "velocity-x", "velocity-y", "velocity-z"});
            break;
        default:
            fp << "Error dimension.";
            fp.close();
            return;
    }
    // write node
    GEOM::tecplot_node_write(fp, mesh);
    // write data
    int count;
    count = 0;
    fp << std::endl << "## density" << std::endl;
    for (auto &cell : CELLS) {
        fp << "\t" << std::setprecision(DATA_PRECISION) << cell.density;
        if (count++ >= LINE_DATA_NUM) {
            fp << std::endl;
            count = 0;
        }
    }
    count = 0;
    fp << std::endl << "## velocity-x" << std::endl;
    for (auto &cell : CELLS) {
        fp << "\t" << std::setprecision(DATA_PRECISION) << cell.velocity.x;
        if (count++ >= LINE_DATA_NUM) {
            fp << std::endl;
            count = 0;
        }
    }
    count = 0;
    fp << std::endl << "## velocity-y" << std::endl;
    for (auto &cell : CELLS) {
        fp << "\t" << std::setprecision(DATA_PRECISION) << cell.velocity.y;
        if (count++ >= LINE_DATA_NUM) {
            fp << std::endl;
            count = 0;
        }
    }
    if (mesh.dimension() == 3) {
        count = 0;
        fp << std::endl << "## velocity-z" << std::endl;
        for (auto &cell : CELLS) {
            fp << "\t" << std::setprecision(DATA_PRECISION) << cell.velocity.z;
            if (count++ >= LINE_DATA_NUM) {
                fp << std::endl;
                count = 0;
            }
        }
    }
    // write geom
    fp << std::endl << "## Geom" << std::endl;
    for (auto &cell : mesh.CELLS) {
        fp << GEOM::tecplot_cell_format(cell);
    }
    // write end
    fp.close();
    pprint::highlight << "Step = " << step << "  Save file: " << ss.str() << "\n";
    pprint::highlight(prefix);
}

void Scheme::do_residual() {
    VecDouble result(CELLS[0].residual.size(), 0.0);
    // 遍历控制体
    for (auto &cell : CELLS) {
        double weight = mesh.CELLS[cell.id].volume / mesh_total_volume;
        // density
        cell.residual[0].update(cell.density);
        result[0] += weight * cell.residual[0].compute();
        // velocity-x
        cell.residual[1].update(cell.velocity.x);
        result[1] += weight * cell.residual[1].compute();
        // velocity-y
        cell.residual[2].update(cell.velocity.y);
        result[2] += weight * cell.residual[2].compute();
        if (mesh.dimension() == 3) {
            // velocity-z
            cell.residual[3].update(cell.velocity.z);
            result[3] += weight * cell.residual[3].compute();
        }
    }
    pprint::info << "Residual:   step = " << step;
    pprint::info(prefix);
    switch (mesh.dimension()) {
        case 2:
            pprint::info << output_data_to_console({"density", "velocity-x", "velocity-y"},
                                                   result);
            pprint::info();
            break;
        case 3:
            pprint::info << output_data_to_console({"density", "velocity-x", "velocity-y", "velocity-z"},
                                                   result);
            pprint::info();
            break;
        default:
            pprint::warn << "do_residual() error.";
            pprint::warn();
            break;
    }
}

void Scheme::do_crashed(Scheme::Cell &cell) {
    is_crashed = true;
    pprint::warn << "Caught crashed value: " << cell.info();
    pprint::warn(prefix);
}

#endif // SOLVER_DUGKS_INCOMPRESSIBLE
