/**
 * MeshBuild
 *  实现界面构建、几何参数计算的函数
 **/
#include <Mesh.h>

/// 局部变量
const Vec(VecInt) node_set_quad = {
        /// QUAD
        {0, 1}, // face 1
        {1, 2}, // face 2
        {2, 3}, // face 3
        {3, 0}  // face 4
};
const Vec(VecInt) node_set_tria = {
        /// TRIA
        {0, 1}, // face 1
        {1, 2}, // face 2
        {2, 0}  // face 3
};
const Vec(VecInt) node_set_brick = {
        /// BRICK
        {0, 1, 5, 4}, // face 1
        {1, 3, 7, 5}, // face 2
        {3, 2, 6, 7}, // face 3
        {2, 0, 4, 6}, // face 4
        {1, 0, 2, 3}, // face 5
        {4, 5, 7, 6}  // face 6
};
const Vec(VecInt) node_set_wedge = {
        /// WEDGE
        {0, 1, 4, 3}, // face 1
        {1, 2, 5, 4}, // face 2
        {2, 0, 3, 5}, // face 3
        {0, 2, 1},    // face 4
        {3, 4, 5}     // face 5
};
const Vec(VecInt) node_set_tetra = {
        /// TETRA
        {1, 0, 2}, // face 1
        {0, 1, 3}, // face 2
        {1, 2, 3}, // face 3
        {2, 0, 3}  // face 4
};
const Vec(VecInt) node_set_pyram = {
        /// PYRAM
        {0, 2, 3, 1},  // face 1
        {0, 1, 4},    // face 2
        {1, 3, 4},    // face 3
        {3, 2, 4},    // face 4
        {2, 0, 4}     // face 5
};

/// 局部函数
void GenerateFace(int elem_type, MESH::Cell &cell, int on_cell_face_id, MapInt &map, Vec(MESH::Interface) &FACES);

void BuildInterface2D(MESH::Mesh &mesh);

void BuildInterface3D(MESH::Mesh &mesh);

void LinkFaceToMark(MESH::Mesh &mesh);

void LinkNearCell(MESH::Mesh &mesh);


/// Mesh.h 功能实现
using namespace MESH;

void Mesh::BuildInterface() {
    int D = dimension();
    // 生成界面
    switch (D) {
        case 2:
            BuildInterface2D(*this);
            break;
        case 3:
            BuildInterface3D(*this);
            break;
        default:
            pprint::error << "MESH::BuildInterface() caught unsupported dimension.";
            pprint::error("Mesh");
            throw std::invalid_argument("MESH::BuildInterface() caught unsupported dimension.");
    }
    // 绑定界面与边界
    LinkFaceToMark(*this);
    // 绑定公用界面的单元体
    LinkNearCell(*this);
    status = 2;

    pprint::info << "Built " << D << "D-interfaces successfully.";
    pprint::info("Mesh");
}

void Mesh::BuildCellGeom() {
    // 遍历单元体
    for (auto &cell : CELLS) {
        // 计算格心坐标
        cell.pos = compute_position(cell, NODES);
        // 计算体积
        cell.volume = compute_volume(cell, NODES);
        double pos_mod = cell.pos.magnitude();
        /// 为 PHY 网格时
        double ssc = pow((Pi_over_2 / NDFCD) * cell.volume, 1.0 / NDFCD);
        update_min_mesh_size(ssc);
        /// 为 DVS 网格时
        update_max_discrete_velocity(pos_mod);
        cell.pos_square = cell.pos * cell.pos;
        }
    status = 3;
    pprint::info << "Built " << dimension() << "D-cell-geom successfully.";
    pprint::info("Mesh");
}

void Mesh::BuildFaceGeom() {
    // 遍历界面
    for (auto &face : FACES) {
        // 计算格心坐标
        face.pos = compute_position(face, NODES);
        // 计算面积
        face.area = compute_area(face, NODES);
    }
    status = 4;
    pprint::info << "Built " << dimension() << "D-face-geom successfully.";
    pprint::info("Mesh");
}

void Mesh::BuildNormalVector() {
    int D = dimension();
    // 遍历界面
    for (auto &face : FACES) {
        // 引用界面直接绑定的单元体
        if (face.on_cell_id < 0) throw std::invalid_argument("Interface did not link to a cell.");
        auto &cell = CELLS[face.on_cell_id];
        // 获取单元指向界面单位向量
        Vector vcf = face.pos - cell.pos;
        // 根据维度分类
        if (D == 2) {
            /// 2D case - z = 0
            // 获取界面共线单位向量
            Vector vfn = NODES[face.node_id[0]].pos - NODES[face.node_id[1]].pos;
            vfn.norm();
            // 计算 vfn 与 x, y 轴夹角
            double vfn_x_axis = GEOM::angle_with_X_axis(vfn),
                   vfn_y_axis = GEOM::angle_with_Y_axis(vfn),
                   vn_x_axis;
            if (vfn_y_axis <= Pi_over_2) {
                vn_x_axis = Pi_over_2 + vfn_x_axis;
            } else {
                vn_x_axis = Pi_over_2 - vfn_x_axis;
            }
            // 生成法向量
            Vector nv = {cos(vn_x_axis), sin(vn_x_axis), 0.0};
            // 判断正负
            if (nv * vcf >= 0.0) {
                face.on_cell_nv = nv;
                face.inv_cell_nv = -nv;
            } else {
                face.on_cell_nv = -nv;
                face.inv_cell_nv = nv;
            }
        } else if (D == 3) {
            /// 3D case - z != 0
            // 获取界面相接两边向量
            Vector v01 = NODES[face.node_id[1]].pos - NODES[face.node_id[0]].pos,
                   v12 = NODES[face.node_id[2]].pos - NODES[face.node_id[1]].pos;
            // 获取界面单位法向量
            Vector nv = v01 ^ v12;
            nv.norm();
            // 判断正负
            if (nv * vcf >= 0.0) {
                face.on_cell_nv = nv;
                face.inv_cell_nv = -nv;
            } else {
                face.on_cell_nv = -nv;
                face.inv_cell_nv = nv;
            }
        } else throw std::invalid_argument("Mesh::BuildNormalVector() caught invalid dimension.");
    }
    status = 5;
    pprint::info << "Built " << D << "D-normal-vectors successfully.";
    pprint::info("Mesh");
}

void Mesh::BuildMesh() {
    double cost = clock();
    switch (type) {
        case MeshTypeDVS_ParseAsPHY:
        case MeshTypePHY:     // 物理空间网格
            BuildCellGeom();
            BuildInterface();
            BuildFaceGeom();
            BuildNormalVector();
            break;
        case MeshTypeDVS:     // 速度空间网格
            BuildCellGeom();
            break;
        default:
            throw std::invalid_argument("BuildMesh() caught unsupported mesh type.");
    }
    // release unused vec
    for (auto &cell : CELLS) cell.shrink_to_fit();
    shrink_to_fit();
    cost = (clock() - cost) / 1000.0;
    pprint::note << "Build mesh(" << name << ") successfully, cost: " << cost << " sec";
    pprint::note("Mesh");
}

/// 局部函数功能实现
void GenerateFace(int elem_type, Cell &cell, int on_cell_face_id, MapInt &map, Vec(Interface) &FACES) {
    // 生成 node_set
    VecInt node_set;
    node_set.reserve(GEOM::node_num(elem_type));
    const Vec(VecInt) *node_set_elem = nullptr;
    switch (cell.type) {
        case GEOM::QUAD:
            node_set_elem = &node_set_quad;
            break;
        case GEOM::TRIA:
            node_set_elem = &node_set_tria;
            break;
        case GEOM::BRICK:
            node_set_elem = &node_set_brick;
            break;
        case GEOM::WEDGE:
            node_set_elem = &node_set_wedge;
            break;
        case GEOM::TETRA:
            node_set_elem = &node_set_tetra;
            break;
        case GEOM::PYRAM:
            node_set_elem = &node_set_pyram;
            break;
        default:
            throw std::invalid_argument("GenerateFace() caught unsupported geom type.");
    }
    for (int it : (*node_set_elem)[on_cell_face_id]) node_set.push_back(cell.node_id[it]);

    // 生成界面 key
    std::string key = generate_key_from_node_set(node_set);
    /// pprint::debug << "key = " << key;
    /// pprint::debug("GenerateFace");
    // 检测是否已存在
    auto it = map.find(key);
    if (it != map.end()) {
        // 存在，将 face 与 cell 建立关系
        auto &face = FACES[map[key]];
        // 绑定界面所在控制体
        face.inv_cell_id = cell.id;
        face.inv_cell_face = on_cell_face_id;
        cell.face_id[on_cell_face_id] = face.id;
        /// pprint::debug << "face " << map[key] << " existed.";
        /// pprint::debug("GenerateFace");
        return;
    }
    // 不存在，创建界面
    FACES.emplace_back(FACES.size(), elem_type, node_set);
    auto &face = FACES.back();
    // 绑定界面所在控制体
    face.on_cell_id = cell.id;
    face.on_cell_face = on_cell_face_id;
    face.inv_cell_id = cell.id;
    face.inv_cell_face = on_cell_face_id;
    cell.face_id[on_cell_face_id] = face.id;
    // 哈希表登记界面信息
    map[key] = face.id;
}

void BuildInterface2D(Mesh &mesh) {
    // 定义界面哈希表   2 维情况只存在 edge 界面
    MapInt face_map_edge;
    // 遍历所有控制体
    for (auto &cell : mesh.CELLS) {
        switch (cell.type) {
            // 控制体可接受 QUAD TRIA
            case GEOM::QUAD:
                // Quad - face 1
                GenerateFace(GEOM::EDGE, cell, 0,
                             face_map_edge, mesh.FACES);
                // Quad - face 2
                GenerateFace(GEOM::EDGE, cell, 1,
                             face_map_edge, mesh.FACES);
                // Quad - face 3
                GenerateFace(GEOM::EDGE, cell, 2,
                             face_map_edge, mesh.FACES);
                // Quad - face 4
                GenerateFace(GEOM::EDGE, cell, 3,
                             face_map_edge, mesh.FACES);
                break;
            case GEOM::TRIA:
                // Tria - face 1
                GenerateFace(GEOM::EDGE, cell, 0,
                             face_map_edge, mesh.FACES);
                // Tria - face 2
                GenerateFace(GEOM::EDGE, cell, 1,
                             face_map_edge, mesh.FACES);
                // Tria - face 3
                GenerateFace(GEOM::EDGE, cell, 2,
                             face_map_edge, mesh.FACES);
                break;
            default:
                throw std::invalid_argument("BuildInterface2D() caught invalid geom type.");
        }
    }
}

void BuildInterface3D(Mesh &mesh) {
    // 定义界面哈希表   3 维情况存在 quad, tria 界面
    MapInt face_map_quad;
    MapInt face_map_tria;
    // 遍历所有控制体
    for (auto &cell : mesh.CELLS) {
        switch (cell.type) {
            case GEOM::BRICK:
                // Brick - face 1
                GenerateFace(GEOM::QUAD, cell, 0, face_map_quad, mesh.FACES);
                // Brick - face 2
                GenerateFace(GEOM::QUAD, cell, 1, face_map_quad, mesh.FACES);
                // Brick - face 3
                GenerateFace(GEOM::QUAD, cell, 2, face_map_quad, mesh.FACES);
                // Brick - face 4
                GenerateFace(GEOM::QUAD, cell, 3, face_map_quad, mesh.FACES);
                // Brick - face 5
                GenerateFace(GEOM::QUAD, cell, 4, face_map_quad, mesh.FACES);
                // Brick - face 6
                GenerateFace(GEOM::QUAD, cell, 5, face_map_quad, mesh.FACES);
                break;
            case GEOM::WEDGE:
                // Wedge - face 1
                GenerateFace(GEOM::QUAD, cell, 0, face_map_quad, mesh.FACES);
                // Wedge - face 2
                GenerateFace(GEOM::QUAD, cell, 1, face_map_quad, mesh.FACES);
                // Wedge - face 3
                GenerateFace(GEOM::QUAD, cell, 2, face_map_quad, mesh.FACES);
                // Wedge - face 4
                GenerateFace(GEOM::TRIA, cell, 3, face_map_tria, mesh.FACES);
                // Wedge - face 5
                GenerateFace(GEOM::TRIA, cell, 4, face_map_tria, mesh.FACES);
                break;
            case GEOM::TETRA:
                // Tetra - face 1
                GenerateFace(GEOM::TRIA, cell, 0, face_map_tria, mesh.FACES);
                // Tetra - face 2
                GenerateFace(GEOM::TRIA, cell, 1, face_map_tria, mesh.FACES);
                // Tetra - face 3
                GenerateFace(GEOM::TRIA, cell, 2, face_map_tria, mesh.FACES);
                // Tetra - face 4
                GenerateFace(GEOM::TRIA, cell, 3, face_map_tria, mesh.FACES);
                break;
            case GEOM::PYRAM:
                // Pyram - face 1
                GenerateFace(GEOM::QUAD, cell, 0, face_map_quad, mesh.FACES);
                // Pyram - face 2
                GenerateFace(GEOM::TRIA, cell, 1, face_map_tria, mesh.FACES);
                // Pyram - face 3
                GenerateFace(GEOM::TRIA, cell, 2, face_map_tria, mesh.FACES);
                // Pyram - face 4
                GenerateFace(GEOM::TRIA, cell, 3, face_map_tria, mesh.FACES);
                // Pyram - face 5
                GenerateFace(GEOM::TRIA, cell, 4, face_map_tria, mesh.FACES);
                break;
            default:
                throw std::invalid_argument("BuildInterface3D() caught invalid geom type.");
        }
    }
}

void LinkFaceToMark(MESH::Mesh &mesh) {
    if (mesh.MARKS.empty()) return;
    // 遍历 MARKS
    for (auto &mark : mesh.MARKS) {
        // 遍历 MARK_ELEM
        for (auto &mark_elem : mark.MARK_ELEM) {
            // 取出 mark_elem 对应的 cell 和 face
            auto& cell = mesh.CELLS[mark_elem.cell_id];
            auto& face = mesh.FACES[cell.face_id[mark_elem.on_cell_face_id]];
            // 绑定
            face.mark_id = mark.id;
            mark_elem.face_id = face.id;
        }
    }
}

void LinkNearCell (MESH::Mesh &mesh) {
    ///* by face
    Vec(VecInt) face_map(mesh.FACES.size(), VecInt(0));
    for (auto &cell : mesh.CELLS) {
        for (int face_id : cell.face_id) {
            face_map[face_id].push_back(cell.id);
        }
    }
    for (auto &it : face_map) {
        for (int cell_id : it) {
            auto &cell = mesh.CELLS[cell_id];
            for (int near_id : it) {
                if (cell_id != near_id) cell.near_cell_id.push_back(near_id);
            }
        }
    }

    /*
    Vec(VecInt) node_map(mesh.NUMNP, VecInt(0));
    // 遍历控制体
    for (auto &cell : mesh.CELLS) {
        for (int node_id : cell.node_id) {
            node_map[node_id].push_back(cell.id);
        }
    }
    // 绑定控制体
    for (auto &it : node_map) {
        for (int cell_id : it) {
            auto  &cell = mesh.CELLS[cell_id];
            for (int near_id : it) {
                if (near_id != cell_id) cell.near_cell_id.push_back(near_id);
            }
        }
    }
     */
}
