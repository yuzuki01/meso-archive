/**
 * 求解器类构造函数
 */

#include <Solver.h>


/// 网格生成
MESH::Mesh GenerateMeshFromConfig(ConfigReader &reader, int mesh_type) {
    switch (mesh_type) {
        case MeshTypePHY: {
            MESH::Mesh mesh(mesh_type, reader["PHY_MESH"]);
            std::stringstream ss;
            ss << "./case/" << reader["CASE_NAME"] << "/" << reader["PHY_MESH"];
            MESH::NEUReader neuReader(ss.str());
            neuReader.parse(mesh);
            mesh.BuildMesh();
            return mesh;
        }
        case MeshTypeDVS: {
            MESH::Mesh mesh(mesh_type, reader["DVS_MESH"]);
            if (reader["DVS_MESH"] == "NULL") {
                pprint::error << "Caught NULL when reading mesh from config.";
                pprint::error("Mesh");
                throw std::invalid_argument("Caught NULL when reading mesh from config.");
            }
            std::stringstream ss;
            ss << "./case/" << reader["CASE_NAME"] << "/" << reader["DVS_MESH"];
            MESH::NEUReader neuReader(ss.str());
            neuReader.parse(mesh);
            mesh.BuildMesh();
            return mesh;
        }
        default:
            pprint::error << "GenerateMeshFromConfig() caught invalid mesh type.";
            pprint::error("Solver");
            throw std::invalid_argument("GenerateMeshFromConfig() caught invalid mesh type.");
    }
}


MESH::Mesh GenerateMeshGH(int dimension, double RT) {
    // D1Q3
    const double pv[3] = {-sqrt(3.0 * RT), 0.0, sqrt(3.0 * RT)};
    const double pw[3] = {1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0};
    MESH::Mesh mesh(MeshTypeDVS, "DVS-GH");
    switch (dimension) {
        case 2: {
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    Vector velocity(pv[i], pv[j], 0.0);
                    double wi = pw[i] * pw[j];
                    mesh.CELLS.emplace_back(velocity, wi, mesh.CELLS.size());
                    mesh.update_max_discrete_velocity(velocity.magnitude());
                }
            }
            break;
        }
        case 3: {
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++)
                    for (int k = 0; k < 3; k++) {
                        Vector velocity(pv[i], pv[j], pv[k]);
                        double wi = pw[i] * pw[j] * pw[k];
                        mesh.CELLS.emplace_back(velocity, wi, mesh.CELLS.size());
                        mesh.update_max_discrete_velocity(velocity.magnitude());
                    }
            }
            break;
        }
        default:
            throw std::invalid_argument("GenerateMeshGH() caught invalid dimension.");
    }
    int Q = 1, D = dimension;
    for (int i = 0; i < D; i++) Q *= 3;
    std::stringstream ss;
    ss << "D" << D << "Q" << Q;
    mesh.name = ss.str();
    mesh.status = 3;
    mesh.NELEM = mesh.CELLS.size();
    mesh.NDFVL = mesh.NDFCD = dimension;
    mesh.shrink_to_fit();
    pprint::note << "Generate " << ss.str() << " DVS mesh.";
    pprint::note("Mesh");
    return mesh;
}


/// 算法类构造函数
DUGKS_INCOMPRESSIBLE::DUGKS_INCOMPRESSIBLE(ConfigReader &reader) {
    /// 参数赋值
    case_name = reader["CASE_NAME"];
    continue_to_run = false;
    is_crashed = false;
    step = 0;
    Re = stod(reader["Re"]);
    Ma = stod(reader["Ma"]);
    RT = stod(reader["R"]) * stod(reader["TEMPERATURE"]);
    rho = stod(reader["DENSITY"]);
    L = stod(reader["LENGTH"]);
    /// 物理网格
    mesh = GenerateMeshFromConfig(reader, MeshTypePHY);
    DVS = GenerateMeshGH(mesh.dimension(), RT);
    for (auto &mark : reader.marks) {
        mesh.set_mark_params(mark);
    }
    mesh_total_volume = 0.0;
    /// 计算参数
    mesh_cell_num = mesh.CELLS.size();
    mesh_face_num = mesh.FACES.size();
    tau = rho * Ma / Re * L / sqrt(RT);
    dt = stod(reader["CFL"]) * mesh.min_mesh_size / DVS.max_discrete_velocity;
    half_dt = dt / 2.0;
    /// file
    bool result = create_dir("./case/" + case_name + "/result") && create_dir("./case/" + case_name + "/cache");
    if (!result) {
        pprint::error << "Cannot create dir.";
        pprint::error(prefix);
        throw std::invalid_argument("Cannot create dir.");
    }
    /// OpenMp
    int thread_num = stoi(reader["THREAD_NUM"]);
    omp_set_num_threads(thread_num);
    pprint::note << "Set OpenMP thread num = " << thread_num
                 << " (available:" << omp_get_num_procs() << ").";
    pprint::note(prefix);
}

void DUGKS_INCOMPRESSIBLE::info() {
    pprint::note << "Solver info:";
    pprint::note(prefix);
    pprint::info << output_data_to_console({"Re", "Ma", "RT", "tau"},
                                           {Re, Ma, RT, tau});
    pprint::info();
    pprint::info << output_data_to_console({"Density", "Length", "MinMS", "MaxDV"},
                                           {rho, L, mesh.min_mesh_size, DVS.max_discrete_velocity});
    pprint::info();
    mesh.info();
    DVS.info();
}

DUGKS_AOKI::DUGKS_AOKI(ConfigReader &reader) {
    /// 参数赋值
    case_name = reader["CASE_NAME"];
    continue_to_run = false;
    is_crashed = false;
    step = 0;
    Kn = stod(reader["Kn"]);
    R = stod(reader["R"]);
    T = stod(reader["TEMPERATURE"]);
    rho = stod(reader["DENSITY"]);
    L = stod(reader["LENGTH"]);
    /// 网格
    mesh = GenerateMeshFromConfig(reader, MeshTypePHY);
    DVS = GenerateMeshFromConfig(reader, MeshTypeDVS);
    for (auto &mark : reader.marks) {
        mesh.set_mark_params(mark);
    }
    /// Generate shadow cell
    for (auto &mark : mesh.MARKS) {
        if (mark.type == 5 || mark.type == 6) {
            // mark.type is "symmetry" or "periodic"
            for (auto &mark_elem : mark.MARK_ELEM) {
                MESH::construct_shadow_cell(mesh.FACES[mark_elem.face_id], mesh);
            }
        }
    }
    mesh_total_volume = 0.0;
    /// 计算参数
    mesh_cell_num = mesh.CELLS.size();
    mesh_face_num = mesh.FACES.size();
    D = mesh.dimension();
    if (D != 2) throw std::invalid_argument("Solver dugks@aoki caught unsupported dimension.");
    Ac = (2.0 * sqrt(2.0 * R * T / Pi) / (Kn * rho * L));
    dt = stod(reader["CFL"]) * mesh.min_mesh_size / DVS.max_discrete_velocity;
    half_dt = dt / 2.0;
    /// file
    bool result = create_dir("./case/" + case_name + "/result") && create_dir("./case/" + case_name + "/cache");
    if (!result) {
        pprint::error << "Cannot create dir.";
        pprint::error(prefix);
        throw std::invalid_argument("Cannot create dir.");
    }
    /// OpenMp
    int thread_num = stoi(reader["THREAD_NUM"]);
    omp_set_num_threads(thread_num);
    pprint::note << "Set OpenMP thread num = " << thread_num
                 << " (available:" << omp_get_num_procs() << ").";
    pprint::note(prefix);
}

void DUGKS_AOKI::info() {
    pprint::note << "Solver info:";
    pprint::note(prefix);
    pprint::info << output_data_to_console({"Kn", "R", "T", "Ac"},
                                           {Kn, R, T, Ac});
    pprint::info();
    pprint::info << output_data_to_console({"Density", "Length", "MinMS", "MaxDV"},
                                           {rho, L, mesh.min_mesh_size, DVS.max_discrete_velocity});
    pprint::info();
    mesh.info();
    DVS.info();
}
