/******************************************
 *        Gauss Integral Functions        *
 ******************************************/

#include <Solver.h>


/*************************************************
 * GHV - GH 型积分的高斯点，乘于 sqrt(RT) 可以得到速度 *
 * GHW - GH 型积分的权重                           *
 *************************************************/
Map(int, VecDouble) GHV = {
        {3, {-sqrt(3.0), 0.0, sqrt(3.0)}},
};
Map(int, VecDouble) GHW = {
        {3, {1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0}},
};


MESH::Mesh GenerateMesh_GaussHermit(ConfigReader &reader, int mesh_type) {
    double RT = stod(reader["R"]) * stod(reader["TEMPERATURE"]);
    int D = stoi(reader["DIMENSION"]);
    int gauss_point_num = stoi(reader["GAUSS_HERMITE_POINT"]);
    // D1Q3
    auto &pv = GHV[gauss_point_num];
    auto &pw = GHW[gauss_point_num];
    MESH::Mesh mesh(mesh_type, "DVS-GaussHermit");
    switch (D) {
        case 2: {
            for (int i = 0; i < gauss_point_num; i++) {
                for (int j = 0; j < gauss_point_num; j++) {
                    Vector velocity(pv[i], pv[j], 0.0);
                    velocity *= sqrt(RT);
                    double wi = pw[i] * pw[j];
                    int inv_cell_id = gauss_point_num * (gauss_point_num - 1 - i) +
                                      (gauss_point_num - 1 - j);
                    mesh.CELLS.emplace_back(velocity, wi, mesh.CELLS.size(), inv_cell_id);
                    mesh.update_max_discrete_velocity(velocity.magnitude());
                }
            }
            break;
        }
        case 3: {
            for (int i = 0; i < gauss_point_num; i++) {
                for (int j = 0; j < gauss_point_num; j++)
                    for (int k = 0; k < gauss_point_num; k++) {
                        Vector velocity(pv[i], pv[j], pv[k]);
                        velocity *= sqrt(RT);
                        double wi = pw[i] * pw[j] * pw[k];
                        int inv_cell_id = gauss_point_num * gauss_point_num * (gauss_point_num - 1 - i) +
                                          gauss_point_num * (gauss_point_num - 1 - j) +
                                          (gauss_point_num - 1 - k);
                        mesh.CELLS.emplace_back(velocity, wi, mesh.CELLS.size(), inv_cell_id);
                        mesh.update_max_discrete_velocity(velocity.magnitude());
                    }
            }
            break;
        }
        default:
            throw std::invalid_argument("GenerateMesh_Gauss() caught invalid dimension.");
    }
    int Q = 1;
    for (int i = 0; i < D; i++) Q *= gauss_point_num;
    std::stringstream ss;
    ss << "D" << D << "Q" << Q;
    mesh.name = ss.str();
    mesh.status = 3;
    mesh.NELEM = mesh.CELLS.size();
    mesh.NDFVL = mesh.NDFCD = D;
    mesh.shrink_to_fit();
    pprint::note << "Generate " << ss.str() << " DVS mesh.";
    pprint::note("Mesh");
    return mesh;
}
