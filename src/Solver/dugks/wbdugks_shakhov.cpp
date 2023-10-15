/**
 * Well-Balance DUGKS
 *  with Shakhov-BGK model
 **/

#include <Solver.h>


using Scheme = WBDUGKS_SHAKHOV;
using SCell = Scheme::Cell;
using SFace = Scheme::Interface;


/// 构造函数
SCell::Cell(int cell_id, WBDUGKS_SHAKHOV *Solver) {
    id = cell_id;
    solver = Solver;
    int dvs_num = solver->DVS.NELEM;
    density = temperature = 0.0;
    // 最小二乘法
    lsp = GenerateLeastSecondParam(id, solver->mesh);
    // 分布函数
    g_t.resize(dvs_num, 0.0);
    g_bp.resize(dvs_num, 0.0);
    slope_g.resize(dvs_num, Vector(0.0, 0.0, 0.0));
    h_t.resize(dvs_num, 0.0);
    h_bp.resize(dvs_num, 0.0);
    slope_h.resize(dvs_num, Vector(0.0, 0.0, 0.0));
    // shrink to fit
    g_t.shrink_to_fit();
    g_bp.shrink_to_fit();
    slope_g.shrink_to_fit();
    h_t.shrink_to_fit();
    h_bp.shrink_to_fit();
    slope_h.shrink_to_fit();
}

SFace::Interface(int face_id, WBDUGKS_SHAKHOV *Solver) {
    id = face_id;
    solver = Solver;
    int dvs_num = solver->DVS.NELEM;
    density = temperature = 0.0;
    // 分布函数
    g.resize(dvs_num, 0.0);
    g_b.resize(dvs_num, 0.0);
    h.resize(dvs_num, 0.0);
    h_b.resize(dvs_num, 0.0);
    // shrink to fit
    g.shrink_to_fit();
    g_b.shrink_to_fit();
    h.shrink_to_fit();
    h_b.shrink_to_fit();
}

/// 算法初始化
void Scheme::init() {
    int output_count = 0;
    /// 创建对应算法控制体
    for (auto &cell : mesh.CELLS) {
        CELLS.emplace_back(cell.id, this);
        // shadow cell
        if (cell.have_shadow_cell) {
            SHADOW_CELLS.emplace_back(cell.shadow_cell_id, this);
        }
        auto &scheme_cell = CELLS.back();
        // 赋值
        for (int i = 0; i < DVS.NELEM; i++) {
            auto &particle = DVS.CELLS[i];
            double g = g_eq(rho, T, particle.pos * particle.pos);
            double h = h_eq(g, T);
            scheme_cell.g_t[i] = g;
            scheme_cell.h_t[i] = h;
            if (cell.have_shadow_cell) {
                SHADOW_CELLS.back().g_t[i] = g;
                SHADOW_CELLS.back().h_t[i] = h;
            }
        }
        scheme_cell.update_macro_var();
        if (cell.have_shadow_cell) SHADOW_CELLS.back().update_macro_var();
        if (debug_mode && ++output_count <= 10) {
            pprint::debug << "Generate cell: id=" << scheme_cell.id << " density=" << scheme_cell.density
                          << " temperature=" << scheme_cell.temperature;
            pprint::debug(prefix);
        } else if (debug_mode && ++output_count == 12) {
            pprint::debug << "  ......(10 cells info shown)";
            pprint::debug();
        }
        // 残差
        mesh_total_volume += cell.volume;
        scheme_cell.residual.emplace_back(scheme_cell.density);
        scheme_cell.residual.emplace_back(scheme_cell.temperature);
        scheme_cell.residual.emplace_back(scheme_cell.velocity.x);
        scheme_cell.residual.emplace_back(scheme_cell.velocity.y);
        if (mesh.dimension() == 3) scheme_cell.residual.emplace_back(scheme_cell.velocity.z);
        scheme_cell.residual.emplace_back(scheme_cell.heat_flux.x);
        scheme_cell.residual.emplace_back(scheme_cell.heat_flux.y);
        if (mesh.dimension() == 3) scheme_cell.residual.emplace_back(scheme_cell.heat_flux.z);
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
    SHADOW_CELLS.shrink_to_fit();
    /// 变量
    continue_to_run = true;
    pprint::note << "Initialization finished.";
    pprint::note(prefix);
}

/// 平衡态函数
double
Scheme::g_eq(double density, double temperature, double cc) const {
    double RT2 = 2.0 * R * temperature;
    double A = Pi * RT2;
    switch (D) {
        case 2:
            break;
        case 3:
            A = pow(A, 3.0 / 2.0);
            break;
        default:
            throw std::invalid_argument("Scheme::g_eq unsupported dimension.");
    }
    return density * exp(-cc / RT2) / A;
}

double
Scheme::h_eq(double temperature, double geq) const {
    return (K - mesh.dimension() + 3.0) * R * temperature * geq;
}

double
Scheme::g_pr(double density, double temperature, double cc, double cq, double geq) const {
    double RT = R * temperature;
    return geq * (1.0 - Pr) * (cc / RT - D - 2.0) * (cq / RT / RT / density) / 5.0;
}

double
Scheme::h_pr(double density, double temperature, double cc, double cq, double geq) const {
    double RT = R * temperature;
    return geq * (1.0 - Pr) * ((cc / RT - D) * (K - D + 3.0) - 2.0 * K) * (cq / RT / density) / 5.0;
}

/// 松弛时间函数
double Scheme::tau_f(double density, double temperature) const {
    return miu0 * pow(temperature / T, vhs_index) / (density * R * temperature);
}

/// 控制体算法函数
void SCell::update_f_bp() {
    double tau = solver->tau_f(density, temperature);
    double C = 3.0 * solver->half_dt / (2.0 * tau + solver->dt);
    double cc, cq, geq;
    for (int i = 0; i < solver->DVS.NELEM; i++) {
        auto &particle = solver->DVS.CELLS[i];
        cc = (particle.pos - velocity) * (particle.pos - velocity);
        cq = (particle.pos - velocity) * heat_flux;
        geq = solver->g_eq(density, temperature, cc);
        g_bp[i] = (1.0 - C) * g_t[i] +
                  C * (geq + solver->g_pr(density, temperature, cc, cq, geq));
        h_bp[i] = (1.0 - C) * h_t[i] +
                  C * (solver->h_eq(temperature, geq) + solver->h_pr(density, temperature, cc, cq, geq));
    }
    // shadow mesh
    auto &cell = solver->mesh.CELLS[id];
    if (cell.have_shadow_cell) {
        solver->SHADOW_CELLS[cell.shadow_cell_id].g_bp = g_bp;
        solver->SHADOW_CELLS[cell.shadow_cell_id].h_bp = h_bp;
    }
}

void SCell::update_gradient() {
    auto &cell = solver->mesh.CELLS[id];
    int near_num = cell.near_cell_id.size();
    for (int i = 0; i < solver->DVS.NELEM; i++) {
        Vector Sgr(0.0, 0.0, 0.0), Shr(0.0, 0.0, 0.0);
        for (int j = 0; j < near_num; j++) {
            auto &near_cell = solver->CELLS[cell.near_cell_id[j]];
            Sgr += lsp.weight[j] * (near_cell.g_bp[i] - g_bp[i]) * lsp.dr[j];
            Shr += lsp.weight[j] * (near_cell.h_bp[i] - h_bp[i]) * lsp.dr[j];
        }
        // shadow mesh
        if (cell.have_shadow_cell) {
            auto &near_cell = solver->SHADOW_CELLS[cell.shadow_cell_id];
            Sgr += lsp.weight[near_num] * (near_cell.g_bp[i] - g_bp[i]) * lsp.dr[near_num];
            Shr += lsp.weight[near_num] * (near_cell.h_bp[i] - h_bp[i]) * lsp.dr[near_num];
        }
        switch (solver->mesh.dimension()) {
            case 2:
                slope_g[i] = {lsp.Cx * Sgr, lsp.Cy * Sgr, 0.0};
                slope_h[i] = {lsp.Cx * Shr, lsp.Cy * Shr, 0.0};
                break;
            case 3:
                throw std::invalid_argument("Caught dimension error when computing gradient.");
            default:
                throw std::invalid_argument("update_gradient() dimension error.");
        }
        /// slope_g[i] = slope_h[i] = {0.0, 0.0, 0.0};   // set zero gradient
    }
}

void SCell::update_macro_var() {
    density = temperature = 0.0;
    velocity = heat_flux = {0.0, 0.0, 0.0};
    for (int i = 0; i < solver->DVS.NELEM; i++) {
        auto &particle = solver->DVS.CELLS[i];
        density += particle.volume * g_t[i];
        velocity += particle.volume * g_t[i] * particle.pos;
    }
    velocity /= density;
    for (int i = 0; i < solver->DVS.NELEM; i++) {
        auto &particle = solver->DVS.CELLS[i];
        Vector c = particle.pos - velocity;
        double ccg_h = (c * c) * g_t[i] + h_t[i];
        temperature += particle.volume * ccg_h;
        heat_flux += particle.volume * ccg_h * c;
    }
    temperature = (temperature - density * (velocity * velocity)) / (2.0 * solver->Cv);
    double tau2 = 2.0 * solver->tau_f(density, temperature);
    heat_flux *= tau2 / (tau2 + solver->dt * solver->Pr);

    if (std::isnan(density) || std::isnan(temperature) || temperature <= 0.0) {
        solver->do_crashed(*this);
    }
}

void SCell::update_f_t() {
    // 计算通量
    VecDouble flux_g(solver->DVS.NELEM, 0.0), flux_h(solver->DVS.NELEM, 0.0);
    // 遍历界面
    auto &mesh_cell = solver->mesh.CELLS[id];
    for (int i = 0; i < solver->DVS.NELEM; i++) {
        auto &particle = solver->DVS.CELLS[i];
        for (int face_id : mesh_cell.face_id) {
            auto &face = solver->FACES[face_id];
            auto &mesh_face = solver->mesh.FACES[face_id];
            auto &nv = (mesh_face.on_cell_id == id) ? mesh_face.on_cell_nv
                                                    : mesh_face.inv_cell_nv;
            double C = (particle.pos * nv) * (mesh_face.area / mesh_cell.volume * solver->dt);
            flux_g[i] += C * face.g[i];
            flux_h[i] += C * face.h[i];
        }
    }
    // 更新 f_t
    for (int i = 0; i < solver->DVS.NELEM; i++) {
        g_t[i] = (4.0 * g_bp[i] - g_t[i]) / 3.0 - flux_g[i];
        h_t[i] = (4.0 * h_bp[i] - h_t[i]) / 3.0 - flux_h[i];
    }
}

std::string SCell::info() const {
    std::stringstream ss;
    ss << solver->mesh.CELLS[id].info_with_pos() << " density=" << density << " temperature=" << temperature
       << " velocity=" << velocity.info();
    return ss.str();
}

/// 界面算法函数
void SFace::update_f_b() {
    auto &face = solver->mesh.FACES[id];
    auto &on_cell = solver->CELLS[face.on_cell_id];
    auto &mesh_on_cell = solver->mesh.CELLS[face.on_cell_id];
    auto &inv_cell = solver->CELLS[face.inv_cell_id];
    auto &nv = face.on_cell_nv;
    for (int i = 0; i < solver->DVS.NELEM; i++) {
        auto &particle = solver->DVS.CELLS[i];
        if (particle.pos * nv >= 0.0) {
            Vector dr = face.pos - mesh_on_cell.pos - solver->half_dt * particle.pos;
            g_b[i] = on_cell.g_bp[i] + on_cell.slope_g[i] * dr;
            h_b[i] = on_cell.h_bp[i] + on_cell.slope_h[i] * dr;
        } else if (mesh_on_cell.have_shadow_cell) {
            // shadow mesh
            Vector dr = face.pos - solver->mesh.SHADOW_CELLS[mesh_on_cell.shadow_cell_id].pos -
                        solver->half_dt * particle.pos;
            g_b[i] = on_cell.g_bp[i] + solver->SHADOW_CELLS[mesh_on_cell.shadow_cell_id].slope_g[i] * dr;
            h_b[i] = on_cell.h_bp[i] + solver->SHADOW_CELLS[mesh_on_cell.shadow_cell_id].slope_h[i] * dr;
        } else {
            Vector dr = face.pos - solver->mesh.CELLS[face.inv_cell_id].pos - solver->half_dt * particle.pos;
            g_b[i] = inv_cell.g_bp[i] + inv_cell.slope_g[i] * dr;
            h_b[i] = inv_cell.h_bp[i] + inv_cell.slope_h[i] * dr;
        }
    }
}

void SFace::update_macro_var() {
    density = temperature = 0.0;
    velocity = heat_flux = {0.0, 0.0, 0.0};
    for (int i = 0; i < solver->DVS.NELEM; i++) {
        auto &particle = solver->DVS.CELLS[i];
        density += particle.volume * g_b[i];
        velocity += particle.volume * g_b[i] * particle.pos;
    }
    velocity /= density;
    for (int i = 0; i < solver->DVS.NELEM; i++) {
        auto &particle = solver->DVS.CELLS[i];
        Vector c = particle.pos - velocity;
        double ccg_h = (c * c) * g_b[i] + h_b[i];
        temperature += particle.volume * ccg_h;
        heat_flux += particle.volume * ccg_h * c;
    }
    temperature = (temperature - density * (velocity * velocity)) / (2.0 * solver->Cv);
    double tau2 = 2.0 * solver->tau_f(density, temperature);
    heat_flux *= tau2 / (tau2 + solver->half_dt * solver->Pr);
}

void SFace::update_f() {
    double tau = solver->tau_f(density, temperature);
    double C = solver->half_dt / (2.0 * tau + solver->half_dt);
    double cc, cq, geq;
    for (int i = 0; i < solver->DVS.NELEM; i++) {
        auto &particle = solver->DVS.CELLS[i];
        cc = (particle.pos - velocity) * (particle.pos - velocity);
        cq = (particle.pos - velocity) * heat_flux;
        geq = solver->g_eq(density, temperature, cc);
        g[i] = (1.0 - C) * g_b[i] +
               C * (geq + solver->g_pr(density, temperature, cc, cq, geq));
        h[i] = (1.0 - C) * h_b[i] +
               C * (solver->h_eq(temperature, geq) + solver->h_pr(density, temperature, cc, cq, geq));
    }
}

void SFace::boundary_condition() {
    /// default:
    /// "symmetry"          5
    /// "periodic"          6
    auto &face = solver->mesh.FACES[id];
    switch (face.boundary_type) {
        case 0: /// "interface"         0
            break;
        case 3:
            /// "isothermal_wall"   3
        {
            auto &mark = solver->mesh.MARKS[face.mark_id];
            auto &nv = face.inv_cell_nv;
            double cc, rho_w, rho_w_reflect;
            rho_w = rho_w_reflect = 0.0;
            for (int i = 0; i < solver->DVS.NELEM; i++) {
                auto &particle = solver->DVS.CELLS[i];
                double kn = particle.pos * nv;
                if (kn < 0.0) {
                    rho_w -= particle.volume * kn * g[i];
                } else {
                    cc = (particle.pos - mark.velocity) * (particle.pos - mark.velocity);
                    rho_w_reflect +=
                            particle.volume * kn * solver->g_eq(1.0, mark.temperature, cc);
                }
            }
            rho_w /= rho_w_reflect;
            for (int i = 0; i < solver->DVS.NELEM; i++) {
                auto &particle = solver->DVS.CELLS[i];
                if (particle.pos * nv >= 0.0) {
                    cc = (particle.pos - mark.velocity) * (particle.pos - mark.velocity);
                    g[i] = solver->g_eq(rho_w, mark.temperature, cc);
                    h[i] = solver->h_eq(g[i], mark.temperature);
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
            fp << GEOM::tecplot_file_header(mesh,
                                            {"density", "temperature", "velocity-x", "velocity-y", "heat_flux-x",
                                             "heat_flux-y"});
            break;
        case 3:
            fp << GEOM::tecplot_file_header(mesh,
                                            {"density", "temperature", "velocity-x", "velocity-y", "velocity-z",
                                             "heat_flux-x", "heat_flux-y", "heat_flux-z"});
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
    fp << std::endl << "## temperature" << std::endl;
    for (auto &cell : CELLS) {
        fp << "\t" << std::setprecision(DATA_PRECISION) << cell.temperature;
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
    count = 0;
    fp << std::endl << "## heat_flux-x" << std::endl;
    for (auto &cell : CELLS) {
        fp << "\t" << std::setprecision(DATA_PRECISION) << cell.heat_flux.x;
        if (count++ >= LINE_DATA_NUM) {
            fp << std::endl;
            count = 0;
        }
    }
    count = 0;
    fp << std::endl << "## heat_flux-y" << std::endl;
    for (auto &cell : CELLS) {
        fp << "\t" << std::setprecision(DATA_PRECISION) << cell.heat_flux.y;
        if (count++ >= LINE_DATA_NUM) {
            fp << std::endl;
            count = 0;
        }
    }
    if (mesh.dimension() == 3) {
        count = 0;
        fp << std::endl << "## heat_flux-z" << std::endl;
        for (auto &cell : CELLS) {
            fp << "\t" << std::setprecision(DATA_PRECISION) << cell.heat_flux.z;
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
        // temperature
        cell.residual[1].update(cell.temperature);
        result[1] += weight * cell.residual[1].compute();
        // velocity-x
        cell.residual[2].update(cell.velocity.x);
        result[2] += weight * cell.residual[2].compute();
        // velocity-y
        cell.residual[3].update(cell.velocity.y);
        result[3] += weight * cell.residual[3].compute();
        if (mesh.dimension() == 3) {
            // velocity-z
            cell.residual[4].update(cell.velocity.z);
            result[4] += weight * cell.residual[4].compute();
        }
    }
    pprint::info << "Residual:   step = " << step;
    pprint::info(prefix);
    switch (mesh.dimension()) {
        case 2:
            pprint::info << output_data_to_console(
                    {"density", "temperature", "velocity-x", "velocity-y", "heat_flux-x", "heat_flux-y"},
                    result);
            pprint::info();
            break;
        case 3:
            pprint::info << output_data_to_console(
                    {"density", "temperature", "velocity-x", "velocity-y", "velocity-z", "heat_flux-x", "heat_flux-y",
                     "heat_flux-z"},
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

