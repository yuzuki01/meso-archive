/**
 * DUGKS
 *  with Shakhov-BGK model
 **/

#include <Solver.h>
#ifdef SOLVER_DUGKS_SHAKHOV

using Scheme = DUGKS_SHAKHOV;
using SCell = Scheme::Cell;
using SFace = Scheme::Interface;


/// 构造函数
SCell::Cell(int cell_id, Scheme *Solver) {
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
    // limiter
    if (solver->limiter) {
        Wg_min.resize(dvs_num, 0.0);
        Wg_max.resize(dvs_num, 0.0);
        Wh_min.resize(dvs_num, 0.0);
        Wh_max.resize(dvs_num, 0.0);
        // shrink_to_fit
        Wg_min.shrink_to_fit();
        Wg_max.shrink_to_fit();
        Wh_min.shrink_to_fit();
        Wh_max.shrink_to_fit();
    } else {
        Wg_min.resize(0);
        Wg_max.resize(0);
        Wh_min.resize(0);
        Wh_max.resize(0);
    }
    // shrink to fit
    g_t.shrink_to_fit();
    g_bp.shrink_to_fit();
    slope_g.shrink_to_fit();
    h_t.shrink_to_fit();
    h_bp.shrink_to_fit();
    slope_h.shrink_to_fit();
}

SFace::Interface(int face_id, Scheme *Solver) {
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
        auto &scheme_cell = CELLS.back();
        // 赋值
        for (int i = 0; i < DVS.NELEM; i++) {
            auto &particle = DVS.CELLS[i];
            double gm = g_maxwell(rho, T, particle.pos * particle.pos);
            double hm = h_maxwell(T, gm);
            scheme_cell.g_t[i] = gm;
            scheme_cell.h_t[i] = hm;
        }
        scheme_cell.update_macro_var();
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
    /// 变量
    continue_to_run = true;
    pprint::note << "Initialization finished.";
    pprint::note(prefix);
}

/// 平衡态函数
double Scheme::g_maxwell(double density, double temperature, double cc) const {
    double RT2 = 2.0 * R * temperature;
    double A = Pi * RT2;
    if (D == 3) A = pow(RT2, 3.0 / 2.0);
    return density * exp(-cc / RT2) / A;
}

double Scheme::h_maxwell(double temperature, double gm) const {
    return (K + 3.0 - D) * (R * temperature) * gm;
}

double Scheme::g_shakhov(double density, double temperature, double cc, double cq, double gm) const {
    double RT = R * temperature;
    double p = density * RT;
    return gm * (1.0 + (1.0 - Pr) * (cq / 5.0 / p / RT) * (cc / RT - D - 2.0));
}

double Scheme::h_shakhov(double density, double temperature, double cc, double cq, double gm) const {
    double RT = R * temperature;
    double p = density * RT;
    return h_maxwell(temperature, gm) + (1.0 - Pr) * gm * (cq / 5.0 / p) * ((cc / RT - D) * (K + 3.0 - D) - 2.0 * K);
}

/// 松弛时间函数
double Scheme::tau_f(double density, double temperature) const {
    return miu0 * pow(temperature / T, vhs_index) / (density * R * temperature);
}

/// Ventaka限制器
double Scheme::venkata_limiter(double wi, double wi_max, double wi_min, double delta_i, double volume) const {
    double a, c;
    if (delta_i == 0.0) return 1.0;
    if (delta_i > 0.0) {
        a = wi_max - wi;
    } else {
        a = wi - wi_min;
    }
    c = (limiter_k * limiter_k * limiter_k) * pow(volume, 3.0 / 2.0);
    return (a * a + 2.0 * a * delta_i + c) / (a * a + 2.0 * delta_i * delta_i + a * delta_i + c);
}

/// 控制体算法函数
void SCell::update_f_bp() {
    double C = (solver->half_dt + solver->dt) / (2.0 * solver->tau_f(density, temperature) + solver->dt);
    double cc, cq, gm;
    for (int i = 0; i < solver->DVS.NELEM; i++) {
        auto &particle = solver->DVS.CELLS[i];
        cc = (particle.pos - velocity) * (particle.pos - velocity);
        cq = (particle.pos - velocity) * heat_flux;
        gm = solver->g_maxwell(density, temperature, cc);
        g_bp[i] = (1.0 - C) * g_t[i] +
                  C * solver->g_shakhov(density, temperature, cc, cq, gm);
        h_bp[i] = (1.0 - C) * h_t[i] +
                  C * solver->h_shakhov(density, temperature, cc, cq, gm);
    }
}

void SCell::update_gradient() {
    auto &cell = solver->mesh.CELLS[id];
    int near_num = cell.near_cell_id.size();
    for (int i = 0; i < solver->DVS.NELEM; i++) {
        Vector Sgr(0.0, 0.0, 0.0), Shr(0.0, 0.0, 0.0);
        if (solver->limiter) {
            Wg_max[i] = Wg_min[i] = g_bp[i];
            Wh_max[i] = Wh_min[i] = h_bp[i];
        }
        for (int j = 0; j < near_num; j++) {
            auto &near_cell = solver->CELLS[cell.near_cell_id[j]];
            if (solver->limiter) {
                Wg_max[i] = Wg_max[i] > near_cell.g_bp[i] ? Wg_max[i] : near_cell.g_bp[i];
                Wg_min[i] = Wg_min[i] < near_cell.g_bp[i] ? Wg_min[i] : near_cell.g_bp[i];
                Wh_max[i] = Wh_max[i] > near_cell.h_bp[i] ? Wh_max[i] : near_cell.h_bp[i];
                Wh_min[i] = Wh_min[i] < near_cell.h_bp[i] ? Wh_min[i] : near_cell.h_bp[i];
            }
            Sgr += lsp.weight[j] * (near_cell.g_bp[i] - g_bp[i]) * lsp.dr[j];
            Shr += lsp.weight[j] * (near_cell.h_bp[i] - h_bp[i]) * lsp.dr[j];
        }
        slope_g[i] = {lsp.Cx * Sgr, lsp.Cy * Sgr, lsp.Cz * Sgr};
        slope_h[i] = {lsp.Cx * Shr, lsp.Cy * Shr, lsp.Cz * Shr};
        /// slope_g[i] = slope_h[i] = {0.0, 0.0, 0.0};   // set zero gradient
    }
}

void SCell::update_macro_var() {
    Vector c;
    density = temperature = 0.0;
    velocity = heat_flux = {0.0, 0.0, 0.0};
    for (int i = 0; i < solver->DVS.NELEM; i++) {
        auto &particle = solver->DVS.CELLS[i];
        density += particle.volume * g_t[i];
        velocity += particle.volume * g_t[i] * particle.pos;
        temperature += particle.volume * (particle.pos_square * g_t[i] + h_t[i]);
    }
    velocity /= density;
    temperature = ((temperature / density) - (velocity * velocity)) / (2.0 * solver->Cv);
    for (int i = 0; i < solver->DVS.NELEM; i++) {
        auto &particle = solver->DVS.CELLS[i];
        c = particle.pos - velocity;
        heat_flux += particle.volume * ((c * c) * g_t[i] + h_t[i]) * c;
    }
    double tau = solver->tau_f(density, temperature);
    heat_flux *= tau / (2.0 * tau + solver->dt * solver->Pr);

    if (std::isnan(density)) {
        solver->do_crashed(*this);
    }
    if (temperature < 0.0) {
        pprint::warn << "Caught negative temperature" << info();
        pprint::warn();
    }
}

void SCell::update_f_t() {
    const int dvs_num = solver->DVS.NELEM;
    // 计算通量
    VecDouble flux_g(dvs_num, 0.0), flux_h(dvs_num, 0.0);
    // 遍历界面
    double Adt_V, knAdt_V;
    auto &mesh_cell = solver->mesh.CELLS[id];
    for (int face_id : mesh_cell.face_id) {
        auto &face = solver->FACES[face_id];
        auto &mesh_face = solver->mesh.FACES[face_id];
        auto &nv = (mesh_face.on_cell_id == id) ? mesh_face.on_cell_nv
                                                : mesh_face.inv_cell_nv;
        Adt_V = mesh_face.area / mesh_cell.volume * solver->dt;
        for (int i = 0; i < dvs_num; i++) {
            auto &particle = solver->DVS.CELLS[i];
            knAdt_V = (particle.pos * nv) * Adt_V;
            flux_g[i] += knAdt_V * face.g[i];
            flux_h[i] += knAdt_V * face.h[i];
        }
    }
    // 更新 f_t
    for (int i = 0; i < dvs_num; i++) {
        g_t[i] = (4.0 * g_bp[i] - g_t[i]) / 3.0 - flux_g[i];
        h_t[i] = (4.0 * h_bp[i] - h_t[i]) / 3.0 - flux_h[i];
    }
}

std::string SCell::info() const {
    std::stringstream ss;
    ss << solver->mesh.CELLS[id].info_with_pos() << " density=" << density << " temperature=" << temperature
       << " velocity=" << velocity.info() << " heat_flux=" << heat_flux.info();
    return ss.str();
}

/// 界面算法函数
void SFace::update_f_b() {
    Vector dr;
    double g0, h0, delta_g, delta_h;
    double phi_g, phi_h;
    auto &face = solver->mesh.FACES[id];
    auto &on_cell = solver->CELLS[face.on_cell_id];
    auto &mesh_on_cell = solver->mesh.CELLS[face.on_cell_id];
    auto &inv_cell = solver->CELLS[face.inv_cell_id];
    auto &mesh_inv_cell = solver->mesh.CELLS[face.inv_cell_id];
    auto &nv = face.on_cell_nv;
    for (int i = 0; i < solver->DVS.NELEM; i++) {
        auto &particle = solver->DVS.CELLS[i];
        if (particle.pos * nv >= 0.0) {
            dr = face.pos - mesh_on_cell.pos - solver->half_dt * particle.pos;
            g0 = on_cell.g_bp[i];
            h0 = on_cell.h_bp[i];
            delta_g = on_cell.slope_g[i] * dr;
            delta_h = on_cell.slope_h[i] * dr;
            if (solver->limiter) {
                phi_g = solver->venkata_limiter(on_cell.g_bp[i], on_cell.Wg_max[i], on_cell.Wg_min[i],
                                                delta_g, mesh_on_cell.volume);
                phi_h = solver->venkata_limiter(on_cell.h_bp[i], on_cell.Wh_max[i], on_cell.Wh_min[i],
                                                delta_h, mesh_on_cell.volume);
            } else {
                phi_g = phi_h = 1.0;
            }
        } else {
            dr = face.pos - mesh_inv_cell.pos - solver->half_dt * particle.pos;
            g0 = inv_cell.g_bp[i];
            h0 = inv_cell.h_bp[i];
            delta_g = inv_cell.slope_g[i] * dr;
            delta_h = inv_cell.slope_h[i] * dr;
            if (solver->limiter) {
                phi_g = solver->venkata_limiter(inv_cell.g_bp[i], inv_cell.Wg_max[i], inv_cell.Wg_min[i],
                                                delta_g, mesh_inv_cell.volume);
                phi_h = solver->venkata_limiter(inv_cell.h_bp[i], inv_cell.Wh_max[i], inv_cell.Wh_min[i],
                                                delta_h, mesh_inv_cell.volume);
            } else {
                phi_g = phi_h = 1.0;
            }
        }
        g_b[i] = g0 + delta_g * phi_g;
        h_b[i] = h0 + delta_h * phi_h;
    }
}

void SFace::update_macro_var() {
    Vector c;
    density = temperature = 0.0;
    velocity = heat_flux = {0.0, 0.0, 0.0};
    for (int i = 0; i < solver->DVS.NELEM; i++) {
        auto &particle = solver->DVS.CELLS[i];
        density += particle.volume * g_b[i];
        velocity += particle.volume * g_b[i] * particle.pos;
        temperature += particle.volume * (particle.pos_square * g_b[i] + h_b[i]);
    }
    velocity /= density;
    temperature = ((temperature / density) - (velocity * velocity)) / (2.0 * solver->Cv);
    for (int i = 0; i < solver->DVS.NELEM; i++) {
        auto &particle = solver->DVS.CELLS[i];
        c = particle.pos - velocity;
        heat_flux += particle.volume * ((c * c) * g_b[i] + h_b[i]) * c;
    }
    double tau = solver->tau_f(density, temperature);
    heat_flux *= tau / (2.0 * tau + solver->half_dt * solver->Pr);
}

void SFace::update_f() {
    double C = solver->half_dt / (2.0 * solver->tau_f(density, temperature) + solver->half_dt);
    double cc, cq, geq;
    for (int i = 0; i < solver->DVS.NELEM; i++) {
        auto &particle = solver->DVS.CELLS[i];
        cc = (particle.pos - velocity) * (particle.pos - velocity);
        cq = (particle.pos - velocity) * heat_flux;
        geq = solver->g_maxwell(density, temperature, cc);
        g[i] = (1.0 - C) * g_b[i] +
               C * solver->g_shakhov(density, temperature, cc, cq, geq);
        h[i] = (1.0 - C) * h_b[i] +
               C * solver->h_shakhov(density, temperature, cc, cq, geq);
    }
}

void SFace::boundary_condition() {
    auto &face = solver->mesh.FACES[id];
    switch (face.boundary_type) {
        case BC_INTERFACE:
            return;
        case BC_INLET: {
            double gm, cc;
            Vector c;
            auto &near_cell = solver->CELLS[solver->mesh.FACES[id].on_cell_id];
            auto &mark = solver->mesh.MARKS[face.mark_id];
            auto &nv = face.inv_cell_nv;
            for (int i = 0; i < solver->DVS.NELEM; i++) {
                auto &particle = solver->DVS.CELLS[i];
                c = particle.pos - mark.velocity;
                cc = c * c;
                gm = solver->g_maxwell(mark.density, mark.temperature, cc);
                if (particle.pos * nv >= 0.0) {
                    g[i] = gm;
                    h[i] = solver->h_maxwell(mark.temperature, gm);
                } else {
                    /// f_out = f_eq(W_inlet) + f_neq_near(W_near)
                    c = particle.pos - near_cell.velocity;
                    cc = c * c;
                    double gm_near = solver->g_maxwell(near_cell.density, near_cell.temperature, cc);
                    g[i] = gm + near_cell.g_bp[i] - gm_near;
                    h[i] = solver->h_maxwell(mark.temperature, gm)
                            + near_cell.h_bp[i] - solver->h_maxwell(near_cell.temperature, gm_near);
                }
            }
            return;
        }
        case BC_OUTLET: {
            double gm;
            Vector c;
            auto &near_cell = solver->CELLS[solver->mesh.FACES[id].on_cell_id];
            auto &mark = solver->mesh.MARKS[face.mark_id];
            auto &nv = face.inv_cell_nv;
            for (int i = 0; i < solver->DVS.NELEM; i++) {
                auto &particle = solver->DVS.CELLS[i];
                if (particle.pos * nv >= 0.0) {
                    c = particle.pos - near_cell.velocity;
                    gm = solver->g_maxwell(near_cell.density, near_cell.temperature, c * c);
                    g[i] = gm;
                    h[i] = solver->h_maxwell(near_cell.temperature, gm);
                }
            }
            return;
        }
        case BC_ISOTHERMAL_WALL:
        {
            double gm, cc;
            Vector c;
            auto &near_cell = solver->CELLS[solver->mesh.FACES[id].on_cell_id];
            auto &mark = solver->mesh.MARKS[face.mark_id];
            auto &nv = face.inv_cell_nv;
            double rho_w, rho_w_reflect;
            rho_w = rho_w_reflect = 0.0;
            for (int i = 0; i < solver->DVS.NELEM; i++) {
                auto &particle = solver->DVS.CELLS[i];
                double kn = particle.pos * nv;
                if (kn < 0.0) {
                    rho_w -= particle.volume * kn * g[i];
                } else {
                    cc = (particle.pos - mark.velocity) * (particle.pos - mark.velocity);
                    rho_w_reflect +=
                            particle.volume * kn * solver->g_maxwell(1.0, mark.temperature, cc);
                }
            }
            rho_w /= rho_w_reflect;
            for (int i = 0; i < solver->DVS.NELEM; i++) {
                auto &particle = solver->DVS.CELLS[i];
                if (particle.pos * nv >= 0.0) {
                    c = particle.pos - mark.velocity;
                    cc = c * c;
                    gm = solver->g_maxwell(rho_w, mark.temperature, cc);
                    g[i] = gm;
                    h[i] = solver->h_maxwell(mark.temperature, gm);
                }
            }
            return;
        }
        case BC_ADIABAT_WALL:
        {
            double gm, cc;
            Vector c;
            auto &near_cell = solver->CELLS[solver->mesh.FACES[id].on_cell_id];
            auto &mark = solver->mesh.MARKS[face.mark_id];
            auto &nv = face.inv_cell_nv;
            double rho_w, rho_w_reflect;
            rho_w = rho_w_reflect = 0.0;
            for (int i = 0; i < solver->DVS.NELEM; i++) {
                auto &particle = solver->DVS.CELLS[i];
                double kn = particle.pos * nv;
                if (kn < 0.0) {
                    rho_w -= particle.volume * kn * g[i];
                } else {
                    cc = (particle.pos - mark.velocity) * (particle.pos - mark.velocity);
                    rho_w_reflect +=
                            particle.volume * kn * solver->g_maxwell(1.0, near_cell.temperature, cc);
                }
            }
            rho_w /= rho_w_reflect;
            for (int i = 0; i < solver->DVS.NELEM; i++) {
                auto &particle = solver->DVS.CELLS[i];
                if (particle.pos * nv >= 0.0) {
                    c = particle.pos - mark.velocity;
                    cc = c * c;
                    gm = solver->g_maxwell(rho_w, near_cell.temperature, cc);
                    g[i] = gm;
                    h[i] = solver->h_maxwell(near_cell.temperature, gm);
                }
            }
            return;
        }
        default:
            pprint::error << "Caught unsupported boundary condition.";
            pprint::error(solver->prefix);
            throw std::invalid_argument("Caught unsupported boundary condition.");
    }
}

/// 算法函数
void Scheme::do_step() {
    /// non-forcing DUGKS
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
                                            {"density", "temperature",
                                             "velocity-x", "velocity-y",
                                             "heat_flux-x", "heat_flux-y"});
            break;
        case 3:
            fp << GEOM::tecplot_file_header(mesh,
                                            {"density", "temperature",
                                             "velocity-x", "velocity-y", "velocity-z",
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
        int i = 0;
        // density
        cell.residual[i].update(cell.density);
        result[i] += weight * cell.residual[i].compute();
        i++;
        // temperature
        cell.residual[i].update(cell.temperature);
        result[i] += weight * cell.residual[i].compute();
        i++;
        // velocity-x
        cell.residual[i].update(cell.velocity.x);
        result[i] += weight * cell.residual[i].compute();
        i++;
        // velocity-y
        cell.residual[i].update(cell.velocity.y);
        result[i] += weight * cell.residual[i].compute();
        i++;
        if (D == 3) {
            // velocity-z
            cell.residual[i].update(cell.velocity.z);
            result[i] += weight * cell.residual[i].compute();
            i++;
        }
        // heat_flux-x
        cell.residual[i].update(cell.heat_flux.x);
        result[i] += weight * cell.residual[i].compute();
        i++;
        // heat_flux-y
        cell.residual[i].update(cell.heat_flux.y);
        result[i] += weight * cell.residual[i].compute();
        i++;
        // heat_flux-z
        if (D == 3) {
            cell.residual[i].update(cell.heat_flux.x);
            result[i] += weight * cell.residual[i].compute();
            i++;
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

#endif // SOLVER_DUGKS_SHAKHOV
