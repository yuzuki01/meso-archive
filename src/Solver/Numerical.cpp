#include <Solver.h>


/// 最小二乘法
LeastSecondParam GenerateLeastSecondParam(int cell_id, MESH::Mesh &mesh) {
    LeastSecondParam result;
    const int dimension = mesh.dimension();
    double Sxx, Sxy, Sxz, Syy, Syz, Szz;
    Sxx = Sxy = Sxz = Syy = Syz = Szz = 0.0;
    auto &cell = mesh.CELLS[cell_id];
    int near_cell_num = cell.near_cell_id.size();
    result.dr.resize(near_cell_num + (cell.have_shadow_cell?1:0));
    result.weight.resize(near_cell_num + (cell.have_shadow_cell?1:0));
    for (int i = 0; i < near_cell_num; i++) {
        auto &near_cell = mesh.CELLS[cell.near_cell_id[i]];
        Vector dr = near_cell.pos - cell.pos;
        double wi = 1.0 / (dr * dr);
        result.dr[i] = dr;
        result.weight[i] = wi;
        // sum
        Sxx += wi * dr.x * dr.x;
        Sxy += wi * dr.x * dr.y;
        Sxz += wi * dr.x * dr.z;
        Syy += wi * dr.y * dr.y;
        Syz += wi * dr.y * dr.z;
        Szz += wi * dr.z * dr.z;
    }
    // shadow cell
    if (cell.have_shadow_cell) {
        auto &near_cell = mesh.SHADOW_CELLS[cell.shadow_cell_id];
        Vector dr = near_cell.pos - cell.pos;
        double wi = 1.0 / (dr * dr);
        result.dr[near_cell_num] = dr;
        result.weight[near_cell_num] = wi;
        // sum
        Sxx += wi * dr.x * dr.x;
        Sxy += wi * dr.x * dr.y;
        Sxz += wi * dr.x * dr.z;
        Syy += wi * dr.y * dr.y;
        Syz += wi * dr.y * dr.z;
        Szz += wi * dr.z * dr.z;
    }
    result.dr.shrink_to_fit();
    result.weight.shrink_to_fit();
    switch (dimension) {
        case 2: {
            double FM = Sxx * Syy - Sxy * Sxy;
            result.Cx = {Syy / FM, -Sxy / FM, 0.0};
            result.Cy = {-Sxy / FM, Sxx / FM, 0.0};
            result.Cz = {0.0, 0.0, 0.0};
            return result;
        }
        case 3: {
            double FM = -Sxz * Sxz * Syy + 2.0 * Sxy * Sxz * Syz - Sxx * Syz * Syz - Sxy * Sxy * Szz +
                        Sxx * Syy * Szz;
            double SxyyzNxzyy = Sxy * Syz - Sxz * Syy;
            double SxzyzNxyzz = Sxz * Syz - Sxy * Szz;
            double SyyzzNyzyz = Syy * Szz - Syz * Syz;
            double SxxzzNxzxz = Sxx * Szz - Sxz * Sxz;
            double SxyxzNxxyz = Sxy * Sxz - Sxx * Syz;
            double SxxyyNxyxy = Sxx * Syy - Sxy * Sxy;
            result.Cx = {
                    SyyzzNyzyz / FM, SxzyzNxyzz / FM, SxyyzNxzyy / FM
            };
            result.Cy = {
                    SxzyzNxyzz / FM, SxxzzNxzxz / FM, SxyxzNxxyz / FM
            };
            result.Cz = {
                    SxyyzNxzyy / FM, SxyxzNxxyz / FM, SxxyyNxyxy / FM
            };
            return result;
        }
        default:
            pprint::error << "GenerateLeastSecondParam() caught invalid dimension.";
            pprint::error("Mesh");
            throw std::invalid_argument("GenerateLeastSecondParam() caught invalid dimension.");
    }
}


/// 残差
void Residual::update(double value) {
    old = now;
    now = value;
}

double Residual::compute() {
    return fabs((now - old) / (old == 0.0 ? 1.0 : old));
}
