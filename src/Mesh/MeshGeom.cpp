#include <Mesh.h>


/**
 * Vector 类
 **/

double &Vector::operator[](int i) {
    switch (i) {
        case 0:
            return x;
        case 1:
            return y;
        case 2:
            return z;
        default:
            throw std::out_of_range("Index out of range for Vector.");
    }
}

Vector &Vector::operator=(const Vector &other) {
    x = other.x;
    y = other.y;
    z = other.z;
    return *this;
}

Vector &Vector::operator=(const std::initializer_list<double> &init_list) {
    switch (init_list.size()) {
        case 3:
            z = *(init_list.begin() + 2);
        case 2:
            y = *(init_list.begin() + 1);
        case 1:
            x = *init_list.begin();
            break;
        default:
            throw std::invalid_argument("Vector: initializer_list has a invalid size.");
    }
    return *this;
}

Vector Vector::operator+(const Vector &other) const {
    // 相加
    return Vector(x + other.x, y + other.y, z + other.z);
}

Vector Vector::operator-(const Vector &other) const {
    // 相减
    return Vector(x - other.x, y - other.y, z - other.z);
}

Vector Vector::operator-() {
    // 取反
    return Vector(-x, -y, -z);
}

Vector operator*(double k, const Vector &v) {
    // 乘于系数
    return Vector(v.x * k, v.y * k, v.z * k);
}

double operator*(const Vector &v1, const Vector &v2) {
    // 点乘
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

Vector Vector::operator*(double k) const {
    // 乘于系数
    return Vector(x * k, y * k, z * k);
}

Vector Vector::operator/(double k) const {
    // 除于系数
    if (k != 0.0) {
        return Vector(x / k, y / k, z / k);
    } else {
        return Vector(INFINITY, INFINITY, INFINITY);
    }
}

Vector &Vector::operator+=(Vector &other) {
    x += other.x;
    y += other.y;
    z += other.z;
    return *this;
}

Vector &Vector::operator-=(Vector &other) {
    x -= other.x;
    y -= other.y;
    z -= other.z;
    return *this;
}

Vector &Vector::operator+=(const Vector &other) {
    x += other.x;
    y += other.y;
    z += other.z;
    return *this;
}

Vector &Vector::operator-=(const Vector &other) {
    x -= other.x;
    y -= other.y;
    z -= other.z;
    return *this;
}

Vector &Vector::operator*=(double k) {
    x *= k;
    y *= k;
    z *= k;
    return *this;
}

Vector &Vector::operator/=(double k) {
    x /= k;
    y /= k;
    z /= k;
    return *this;
}

double Vector::operator*(Vector &other) {
    // 点乘
    return x * other.x + y * other.y + z * other.z;
}

Vector Vector::operator^(Vector &other) {
    // 实现叉乘运算
    double result_x = y * other.z - z * other.y;
    double result_y = z * other.x - x * other.z;
    double result_z = x * other.y - y * other.x;
    return Vector(result_x, result_y, result_z);
}

double Vector::magnitude() {
    // 取模
    return sqrt(x * x + y * y + z * z);
}

void Vector::norm() {
    // 归一
    double len = magnitude();
    if (len != 0.0) {
        x /= len;
        y /= len;
        z /= len;
    } else {
        x = INFINITY;
        y = INFINITY;
        z = INFINITY;
    }
}

std::string Vector::info() const {
    std::stringstream ss;
    ss << "<" << std::setprecision(OUT_PRECISION) << x << ", " << y << ", " << z << ">";
    return ss.str();
}

/**
 * GEOM 几何类函数
 **/
int GEOM::node_num(int elem_type) {
    switch (elem_type) {
        case EDGE:
            return 2;
        case QUAD:
        case TETRA:
            return 4;
        case TRIA:
            return 3;
        case BRICK:
            return 8;
        case WEDGE:
            return 6;
        case PYRAM:
            return 5;
        default:
            throw std::invalid_argument("GEOM::node_num() caught unsupported geom type.");
    }
}


int GEOM::face_num(int elem_type) {
    switch (elem_type) {
        case QUAD:
        case TETRA:
            return 4;
        case TRIA:
            return 3;
        case BRICK:
            return 6;
        case WEDGE:
        case PYRAM:
            return 5;
        default:
            throw std::invalid_argument("GEOM::face_num() caught unsupported geom type.");
    }
}

double GEOM::angle_vectors(Vector &v1, Vector &v2) {
    double m1 = v1.magnitude(), m2 = v2.magnitude();
    return acos((v1 * v2) / (m1 * m2));
}

double GEOM::angle_with_X_axis(Vector &vector) {
    Vector unit_x = Vector(1.0, 0.0, 0.0);
    return angle_vectors(vector, unit_x);
}

double GEOM::angle_with_Y_axis(Vector &vector) {
    Vector unit_y = Vector(0.0, 1.0, 0.0);
    return angle_vectors(vector, unit_y);
}

double GEOM::angle_with_Z_axis(Vector &vector) {
    Vector unit_z = Vector(0.0, 0.0, 1.0);
    return angle_vectors(vector, unit_z);
}

std::string GEOM::tecplot_file_header(MESH::Mesh &mesh, const std::initializer_list<std::string> &var_names) {
    int var_num = var_names.size();
    std::stringstream ss;
    ss << "VARIABLES= \"X\", \"Y\"";
    if (mesh.dimension() == 3) ss << ", \"Z\"";
    for (auto &var : var_names) ss << ", \"" << var << "\"";
    ss << std::endl;
    ss << "ZONE N=" << mesh.NUMNP << ", E=" << mesh.NELEM << ", VARLOCATION=([1-" << mesh.dimension()
       << "]=NODAL, ["
       << (var_num == 1 ? std::to_string(mesh.dimension() + 1) : std::to_string(mesh.dimension() + 1) + "-" +
                                                                 std::to_string(mesh.dimension() + var_num))
       << "]=CELLCENTERED), DATAPACKING=BLOCK, ZONETYPE="
       << (mesh.dimension() == 2 ? "FEQUADRILATERAL" : "FEBRICK") << std::endl << std::endl;
    return ss.str();
}

std::string GEOM::tecplot_cell_format(MESH::Cell &cell) {
    VecInt node(cell.node_id.size());
    for (int i = 0; i < cell.node_id.size(); i++) {
        node[i] = cell.node_id[i] + 1;
    }
    std::stringstream ss;
    switch (cell.type) {
        case QUAD:
            ss << " " << node[0] << " " << node[1] << " " << node[2] << " " << node[3];
            break;
        case TRIA:
            ss << " " << node[0] << " " << node[1] << " " << node[2] << " " << node[0];
            break;
        case BRICK:
            ss << " " << node[0] << " " << node[1] << " " << node[3] << " " << node[2]
               << " " << node[4] << " " << node[5] << " " << node[7] << " " << node[6];
            break;
        case WEDGE:
            ss << " " << node[0] << " " << node[1] << " " << node[2] << " " << node[2]
               << " " << node[3] << " " << node[4] << " " << node[5] << " " << node[5];
            break;
        case TETRA:
            ss << " " << node[0] << " " << node[1] << " " << node[2] << " " << node[2]
               << " " << node[3] << " " << node[3] << " " << node[3] << " " << node[3];
            break;
        case PYRAM:
            ss << " " << node[0] << " " << node[1] << " " << node[3] << " " << node[2]
               << " " << node[4] << " " << node[4] << " " << node[4] << " " << node[4];
            break;
        default:
            ss << "Caught unsupported mesh type when handling with " << cell.info_with_pos();
            pprint::error << ss.str();
            pprint::error("tecplot_cell_format");
            throw std::invalid_argument(ss.str());
    }
    ss << std::endl;
    return ss.str();
}

void GEOM::tecplot_node_write(std::ofstream &fp, MESH::Mesh &mesh) {
    int count;
    fp << std::endl;
    // write NODE
    for (int i = 0; i < mesh.dimension(); i++) {
        count = 0;
        fp << std::endl << "## Coordinate " << i + 1 << "/" << mesh.dimension() << std::endl;
        for (auto &node : mesh.NODES) {
            fp << "\t" << std::setprecision(DATA_PRECISION) << node.pos[i];
            if (count++ >= LINE_DATA_NUM) {
                fp << std::endl;
                count = 0;
            }
        }
    }
    fp << std::endl;
}


/// local func
Vector position(Vec(MESH::Node) &node_set) {
    // 声明坐标
    Vector vec(0.0, 0.0, 0.0);
    for (MESH::Node node : node_set) {
        vec += node.pos;
    }
    return vec / double(node_set.size());
}

double area(int elem_type, Vec(MESH::Node) &node_set) {
    // 判断几何体类型
    switch (elem_type) {
        // 1D case
        case GEOM::EDGE: {
            Vector v01 = node_set[1].pos - node_set[0].pos;     // <0->1>
            return v01.magnitude();
        }
            // 2D case
        case GEOM::QUAD: {
            Vector v01 = node_set[1].pos - node_set[0].pos,     // <0->1>
            v03 = node_set[3].pos - node_set[0].pos,     // <0->3>
            v21 = node_set[1].pos - node_set[2].pos,     // <2->1>
            v23 = node_set[3].pos - node_set[2].pos;     // <2->3>
            return (v01 ^ v03).magnitude() / 2. + (v21 ^ v23).magnitude() / 2.;
        }
        case GEOM::TRIA: {
            Vector v01 = node_set[1].pos - node_set[0].pos,     // <0->1>
            v02 = node_set[2].pos - node_set[0].pos;     // <0->2>
            return (v01 ^ v02).magnitude() / 2.;
        }
        default:
            throw std::invalid_argument("area() caught unsupported geom type.");
    }
}

double volume(int elem_type, Vec(MESH::Node) &node_set) {
    // 判断几何体类型
    switch (elem_type) {
        // 2D case
        case GEOM::QUAD: {
            Vector v01 = node_set[1].pos - node_set[0].pos,     // <0->1>
            v03 = node_set[3].pos - node_set[0].pos,     // <0->3>
            v21 = node_set[1].pos - node_set[2].pos,     // <2->1>
            v23 = node_set[3].pos - node_set[2].pos;     // <2->3>
            return (v01 ^ v03).magnitude() / 2. + (v21 ^ v23).magnitude() / 2.;
        }
        case GEOM::TRIA: {
            Vector v01 = node_set[1].pos - node_set[0].pos,     // <0->1>
            v02 = node_set[2].pos - node_set[0].pos;     // <0->2>
            return (v01 ^ v02).magnitude() / 2.;
        }
            // 3D case
        case GEOM::BRICK: {
            // 三棱柱 02-13-46 部分
            // 三棱柱 再拆分出 三棱锥 0-124 和 四棱锥 2-1364
            Vector v01 = node_set[1].pos - node_set[0].pos,
                    v02 = node_set[2].pos - node_set[0].pos,
                    v04 = node_set[4].pos - node_set[0].pos;
            double tetra_0_124 = fabs(v01 * (v02 ^ v04)) / 6.;
            Vector v12 = node_set[2].pos - node_set[1].pos,
                    v13 = node_set[3].pos - node_set[1].pos,
                    v14 = node_set[4].pos - node_set[1].pos,
                    v63 = node_set[3].pos - node_set[6].pos,
                    v64 = node_set[4].pos - node_set[6].pos;
            double pyram_2_1364 = fabs(v12 * (v13 ^ v14)) / 6. + fabs(v12 * (v63 ^ v64)) / 6.;
            // 三棱柱 57-13-46
            Vector v51 = node_set[1].pos - node_set[5].pos,
                    v54 = node_set[4].pos - node_set[5].pos,
                    v57 = node_set[7].pos - node_set[5].pos;
            double tetra_5_147 = fabs(v51 * (v54 ^ v57)) / 6.;
            Vector v17 = node_set[7].pos - node_set[1].pos;
            double pyram_7_1364 = fabs(v17 * (v13 ^ v14)) / 6. + fabs(v17 * (v63 ^ v64)) / 6.;
            return tetra_0_124 + tetra_5_147 + pyram_2_1364 + pyram_7_1364;
        }
        case GEOM::WEDGE: {
            Vector v01 = node_set[1].pos - node_set[0].pos,
                    v02 = node_set[2].pos - node_set[0].pos,
                    v03 = node_set[3].pos - node_set[0].pos;
            double tetra_0_123 = fabs(v01 * (v02 ^ v03)) / 6.;
            Vector v14 = node_set[4].pos - node_set[1].pos,
                    v12 = node_set[2].pos - node_set[1].pos,
                    v54 = node_set[4].pos - node_set[5].pos,
                    v52 = node_set[2].pos - node_set[5].pos,
                    v13 = node_set[3].pos - node_set[1].pos;
            double pyram_3_1254 = fabs(v13 * (v12 ^ v14)) / 6. + fabs(v13 * (v52 ^ v54));
            return tetra_0_123 + pyram_3_1254;
        }
        case GEOM::TETRA: {
            Vector v01 = node_set[1].pos - node_set[0].pos,
                    v02 = node_set[2].pos - node_set[0].pos,
                    v03 = node_set[3].pos - node_set[0].pos;
            return fabs(v01 * (v02 ^ v03)) / 6.;
        }
        case GEOM::PYRAM: {
            Vector v01 = node_set[1].pos - node_set[0].pos,
                    v02 = node_set[2].pos - node_set[0].pos,
                    v04 = node_set[4].pos - node_set[0].pos,
                    v31 = node_set[1].pos - node_set[3].pos,
                    v32 = node_set[2].pos - node_set[3].pos;
            return fabs(v04 * (v01 ^ v02)) / 6. + fabs(v04 * (v31 ^ v32)) / 6.;
        }
        default:
            throw std::invalid_argument("volume() caught unsupported geom type.");
    }
}

/// Mesh.h 全局函数
Vector MESH::compute_position(Cell &cell, Vec(Node) &NODES) {
    Vec(MESH::Node) node_set;
    node_set.reserve(cell.node_id.size());
    for (int id : cell.node_id) {
        node_set.push_back(NODES[id]);
    }
    return position(node_set);
}

Vector MESH::compute_position(Interface &face, Vec(Node) &NODES) {
    Vec(MESH::Node) node_set;
    node_set.reserve(face.node_id.size());
    for (int id : face.node_id) {
        node_set.push_back(NODES[id]);
    }
    return position(node_set);
}

double MESH::compute_volume(Cell &cell, Vec(Node) &NODES) {
    Vec(MESH::Node) node_set;
    node_set.reserve(cell.node_id.size());
    for (int id : cell.node_id) {
        node_set.push_back(NODES[id]);
    }
    return volume(cell.type, node_set);
}

double MESH::compute_area(Interface &face, Vec(Node) &NODES) {
    Vec(MESH::Node) node_set;
    node_set.reserve(face.node_id.size());
    for (int id : face.node_id) {
        node_set.push_back(NODES[id]);
    }
    return area(face.type, node_set);
}

void MESH::construct_shadow_cell(Interface &face, Mesh &mesh) {
    auto &cell = mesh.CELLS[face.on_cell_id];
    int shadow_cell_id = mesh.SHADOW_CELLS.size();
    Vector shadow_cell_pos = 2.0 * face.pos - cell.pos;
    mesh.SHADOW_CELLS.emplace_back(shadow_cell_pos, cell.volume, shadow_cell_id, -1);
    cell.have_shadow_cell = true;
    cell.shadow_cell_id = shadow_cell_id;
}
