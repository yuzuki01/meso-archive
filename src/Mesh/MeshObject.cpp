/**
 * MeshObjects
 *  实现 MESH 中
 *      Node, Cell, Interface, MarkElem, Mark
 **/


#include <Mesh.h>

using namespace MESH;

/// 全局变量
MapInt MESH::MarkType = {
        {"interface",       0},
        {"inlet",           1},
        {"outlet",          2},
        {"isothermal_wall", 3},
        {"adiabat_wall",    4},
        {"symmetry",        5},
        {"periodic",        6}
};

/// 结点 Node
Node::Node(int node_id, const Vector &position) {
    id = node_id;
    pos = position;
}

Node::Node(const std::vector<std::string> &list) {
    int len = list.size();
    if (len >= 3) {
        id = stoi(list[0]) - 1;
        double x = stod(list[1]);
        double y = stod(list[2]);
        double z = 0.0;
        if (len >= 4) z = stod(list[3]);
        pos = {x, y, z};
    } else {
        throw std::invalid_argument("neu_node_formatter: invalid list.size().");
    }
}

/// 单元 Cell
Cell::Cell(const VecStr &list) {
    id = stoi(list[0]) - 1;
    type = stoi(list[1]);
    int node_num = stoi(list[2]);
    for (int i = 0; i < node_num; i++) {
        node_id.push_back(stoi(list[i + 3]) - 1);
    }
    if (int(node_id.size()) != GEOM::node_num(type)) {
        pprint::error << "Cell<id=" << id + 1 << "> caught invalid node id set.";
        pprint::error("MESH::Cell");
        throw std::invalid_argument("Cell caught invalid node id set.");
    }
    volume = -1.;
    have_shadow_cell = false;   // 是否拥有影格子
    shadow_cell_id = -1;        // 对应影格子编号
    // 确认界面个数
    face_num = GEOM::face_num(type);
    face_id.resize(face_num, -1);
}

Cell::Cell(const Vector &particle_velocity, double weight, int cell_id) {
    id = cell_id;
    type = -1;
    volume = weight;
    have_shadow_cell = false;   // 是否拥有影格子
    shadow_cell_id = -1;        // 对应影格子编号
    pos = particle_velocity;
}

void Cell::shrink_to_fit() {
    node_id.shrink_to_fit();
    face_id.shrink_to_fit();
    near_cell_id.shrink_to_fit();
}

std::string Cell::info() const {
    std::stringstream ss;
    ss << "Cell<tec_id=" << id + 1 << ">";
    return ss.str();
}

std::string Cell::info_with_pos() const {
    std::stringstream ss;
    ss << info() << " , position:" << pos.info();
    return ss.str();
}

/// 界面 Interface
Interface::Interface(int face_id, int elem_type, const std::vector<int> &node_id_vec) {
    id = face_id;
    type = elem_type;
    node_id = node_id_vec;
    boundary_type = 0;
    mark_id = -1;
    on_cell_id = -1;
    on_cell_face = -1;
    inv_cell_id = -1;
    inv_cell_face = -1;
    have_shadow_cell = false;   // 是否拥有影格子
    shadow_cell_id = -1;        // 对应影格子编号
    area = -1.;
    // shrink_to_fit
    node_id.shrink_to_fit();
}

Interface::Interface() {
    have_shadow_cell = false;
    id = type = on_cell_id = shadow_cell_id = -1;
    area = -1.;
}

std::string Interface::info() const {
    std::stringstream ss;
    ss << "Interface<tec_id=" << id + 1 << ">";
    return ss.str();
}

std::string Interface::info_all() const {
    std::stringstream ss;
    ss << info() << " on Cell<tec_id=" << on_cell_id + 1 << ",face=" << on_cell_face + 1 << ",nv="
       << on_cell_nv.info() << ">";
    return ss.str();
}

/// 边界元素 MarkElem
MarkElem::MarkElem(const VecStr &list) {
    cell_id = stoi(list[0]) - 1;
    cell_type = stoi(list[1]);
    on_cell_face_id = stoi(list[2]) - 1;
}

/// 边界 Mark
Mark::Mark(const std::vector<std::string> &list) {
    int len = list.size();
    if (len > 0) name = list[0];
    if (len > 1) type = stoi(list[1]);
    if (len > 2) elem_num = stoi(list[2]);
}


void Mark::info() const {
    pprint::info << "      Mark<" << name << ">" << "\n"
                 << "       - elements: " << elem_num << "\n"
                 << "       - type: " << type_name << "<" << type << ">" << "\n"
                 << "       - density: " << density << " temperature: " << temperature << "\n"
                 << "       - velocity: " << velocity.info();
    pprint::info();
}

/// 网格对象 Mesh
void Mesh::load(const std::string &file_path) {
    NEUReader reader(file_path);
    reader.parse(*this);
}

void Mesh::parse_meshInfo(const VecStr &vector) {
    NUMNP = stoi(vector[0]);
    NELEM = stoi(vector[1]);
    NGRPS = stoi(vector[2]);
    NBSETS = stoi(vector[3]);
    NDFCD = stoi(vector[4]);
    NDFVL = stoi(vector[5]);

    // 预分配空间
    NODES.reserve(NUMNP);
    CELLS.reserve(NELEM);
}

int Mesh::dimension() const {
    return NDFCD;
}

void Mesh::info() {
    if (name == "NULL") return;
    std::string type_name, state = "unknown";
    switch (type) {
        case MeshTypePHY:
            type_name = "PHY";
            switch (status) {
                case 1:
                    state = "waiting for this->BuildInterface()";
                    break;
                case 2:
                    state = "waiting for this->BuildCellGeom()";
                    break;
                case 3:
                    state = "waiting for this->BuildFaceGeom()";
                    break;
                case 4:
                    state = "waiting for this->BuildNormalVector()";
                    break;
                case 5:
                    state = "all done";
                    break;
                default:
                    state = "waiting for NEUReader.parse()";
                    break;
            }
            break;
        case MeshTypeDVS:
            type_name = "DVS";
            switch (status) {
                case 1:
                    state = "waiting for this->BuildCellGeom()";
                    break;
                case 3:
                    state = "all done";
                    break;
                default:
                    state = "waiting for NEUReader.parse()";
                    break;
            }
            break;
        default:
            type_name = "ERROR";
            break;
    }
    pprint::note << "Mesh info:";
    pprint::note("Mesh");
    pprint::info << "   STATE = " << state << "\n"
                 << "   FILE = " << name << "\n"
                 << "   TYPE = " << type_name << "\n"
                 << "   NODE = " << NUMNP << "\n"
                 << "   CELL = " << NELEM << "\n"
                 << "   FACE = " << int(FACES.size());
    pprint::info();
    if (type == 0) {
        pprint::info << "   MARK = " << NBSETS;
        pprint::info();
        for (auto &mark : MARKS) mark.info();
    }
}

void Mesh::update_max_discrete_velocity(double value) {
    if (max_discrete_velocity <= 0.0) {
        max_discrete_velocity = value;
    } else {
        max_discrete_velocity = (max_discrete_velocity >= value) ? max_discrete_velocity : value;
    }
}

void Mesh::update_min_mesh_size(double value) {
    if (min_mesh_size <= 0.0) {
        min_mesh_size = value;
    } else {
        min_mesh_size = (min_mesh_size <= value) ? min_mesh_size : value;
    }
}

void Mesh::shrink_to_fit() {
    NODES.shrink_to_fit();
    NUMNP = NODES.size();
    CELLS.shrink_to_fit();
    SHADOW_CELLS.shrink_to_fit();
    NELEM = CELLS.size();
    FACES.shrink_to_fit();
    pprint::highlight << "Done shrink_to_fit() to Mesh object<" << name << ">, contained object addresses changed.";
    pprint::highlight("Mesh");
}


void Mesh::set_mark_params(const BoundaryParam &bc_param) {
    for (auto& mark : MARKS) {
        if (mark.name == bc_param.name) {
            /// match
            pprint::info << "Set mark params for Mark<" << bc_param.name << "> from Mesh<" << name << ">.";
            pprint::info("Mesh");
            mark.type_name = bc_param.type_name;
            mark.density = bc_param.density;
            mark.temperature = bc_param.temperature;
            mark.velocity = bc_param.velocity;
            mark.type = bc_param.type;
            for (auto &mark_elem : mark.MARK_ELEM) {
                FACES[mark_elem.face_id].boundary_type = mark.type;
            }
            return;
        }
    }
    /// not match
    pprint::warn << "Caught mark<" << bc_param.name << "> which is not read from mesh file.";
    pprint::warn("Mesh");
}
