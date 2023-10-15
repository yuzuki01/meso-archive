#ifndef INCLUDE_CORE
#include <Core.h>
#endif

#define INCLUDE_MESH

#define MeshNode Vec(MESH::Node)
#define MeshCell Vec(MESH::Cell)
#define MeshFace Vec(MESH::Interface)
#define MeshMark Vec(MESH::Mark)
#define MeshMlEm Vec(MESH::MarkElem)
#define MeshPoin MESH::Mesh*
#define MeshTypePHY 0
#define MeshTypeDVS 1
#define MeshTypeDVS_ParseAsPHY 2


/// 3 维向量
class Vector {
public:
    double x, y, z;

    Vector(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {}

    Vector() : x(0.0), y(0.0), z(0.0) {}

    ~Vector() = default;

    double &operator[](int i);                      // 赋值
    Vector &operator=(const Vector &other);
    Vector &operator=(const std::initializer_list<double> &init_list);
    Vector operator+(const Vector &other) const;    // 相加
    Vector operator-(const Vector &other) const;    // 相减
    Vector operator-();                             // 取反
    double operator*(Vector &other);                // 点乘
    Vector operator^(Vector &other);                // 叉乘
    Vector operator*(double k) const;               // 乘于系数
    Vector operator/(double k) const;               // 除于系数
    Vector &operator+=(Vector &other);              // 自加
    Vector &operator-=(Vector &other);              // 自减
    Vector &operator+=(const Vector &other);        // 自加
    Vector &operator-=(const Vector &other);        // 自减
    Vector &operator*=(double k);                   // 自乘
    Vector &operator/=(double k);                   // 自除
    double magnitude();                             // 取模
    void norm();                                    // 归一

    std::string info() const;
};

Vector operator*(double k, const Vector &v);            // 乘于系数
double operator*(const Vector &v1, const Vector &v2);   // 点乘


/// 网格
namespace MESH {
    // .neu 文件读取器
    class NEUReader;

    // 网格元素对象
    class Node;

    class Cell;

    class Interface;

    class MarkElem;

    class Mark;

    // 网格容器对象
    class Mesh;

    // Mark.type map
    extern MapInt MarkType;

    // MeshGeom
    Vector compute_position(Cell &cell, Vec(Node) &NODES);
    Vector compute_position(Interface &face, Vec(Node) &NODES);
    double compute_area(Interface &face, Vec(Node) &NODES);
    double compute_volume(Cell &cell, Vec(Node) &NODES);
    void construct_shadow_cell(Interface &face, Mesh &mesh);
}

/// 算例类
struct BoundaryParam {
    std::string name;
    std::string type_name;
    int type;
    double density = 0.0;
    double temperature = 0.0;
    Vector velocity = {0.0, 0.0, 0.0};
};


/// 几何
namespace GEOM {
    /**
     * 1 = Edge             线段
     * 2 = Quadrilateral    四边形
     * 3 = Triangle         三角形
     * 4 = Brick            立方体
     * 5 = Wedge (Prism)    三棱柱
     * 6 = Tetrahedron      三棱锥
     * 7 = Pyramid          四棱锥
     */
    const int EDGE = 1;
    const int QUAD = 2;
    const int TRIA = 3;
    const int BRICK = 4;
    const int WEDGE = 5;
    const int TETRA = 6;
    const int PYRAM = 7;

    int node_num(int elem_type);    // 获取几何体结点个数
    int face_num(int elem_type);    // 获取几何体界面个数
    /// 运算
    double angle_vectors(Vector &v1, Vector &v2);       // 计算两个向量夹角
    double angle_with_X_axis(Vector &vector);     // 计算向量与 x 轴夹角
    double angle_with_Y_axis(Vector &vector);     // 计算向量与 y 轴夹角
    double angle_with_Z_axis(Vector &vector);     // 计算向量与 z 轴夹角
    /// 输出
    std::string tecplot_file_header(MESH::Mesh &mesh, const std::initializer_list<std::string> &var_names);
    void tecplot_node_write(std::ofstream &fp, MESH::Mesh &mesh);
    std::string tecplot_cell_format(MESH::Cell &cell);
}

/// neu 读取器
class MESH::NEUReader : public Reader {
protected:
    const std::string end_mark = "ENDOFSECTION";
public:
    explicit NEUReader(const std::string &file_path) : Reader("NEU Reader", file_path) {};

    void parse(MESH::Mesh &mesh);
};

/// 结点
class MESH::Node {
public:
    int id;                         // 编号
    Vector pos = {0., 0., 0.};  // 坐标
    explicit Node(int node_id, const Vector &position);
    explicit Node(const VecStr &list);
};

/// 单元体
class MESH::Cell {
public:
    int id;                         // 编号
    int type;                       // 类型
    int face_num;                   // 界面个数
    bool have_shadow_cell;          // 是否拥有影格子
    int shadow_cell_id;             // 对应影格子编号
    Vector pos = {0., 0., 0.};      // 坐标
    double pos_square = 0.;         // 适用于 DVS 网格，表示离散速度的二次方
    VecInt node_id;                 // 构成结点编号
    VecInt face_id;                 // 构成界面编号
    VecInt near_cell_id;            // 相邻控制体编号
    double volume;                  // 体积

    // 按照 Gambit 格式     控制体编号|控制体类型|结点个数|结点编号
    explicit Cell(const VecStr &list);
    // 自定义生成离散速度空间构造函数
    Cell(const Vector &particle_velocity, double weight, int cell_id);

    void shrink_to_fit();

    std::string info() const;
    std::string info_with_pos() const;
};

/// 界面
class MESH::Interface {
public:
    int boundary_type;          // 边界类型
    int mark_id;                // 界面所在边界编号
    int id;                     // 编号
    int type;                   // 类型
    int on_cell_id;             // 界面所在控制体编号
    int inv_cell_id ;           // 另一控制体编号
    int on_cell_face;           // 界面在界面所在控制体上的编号
    int inv_cell_face;          // 界面在另一控制体上的编号
    bool have_shadow_cell;      // 是否拥有影格子
    int shadow_cell_id;         // 影格子编号
    Vector on_cell_nv;          // 界面对应 on_cell_face_id 的法向量
    Vector inv_cell_nv;         // 界面对应另一单元体或壁面指向流体域的法向量
    double area;                // 面积
    Vector pos = {0., 0., 0.};  // 坐标
    VecInt node_id;              // 构成结点编号

    // 格式 界面编号|界面类型|结点编号
    Interface(int face_id, int elem_type, const VecInt &node_id_vec);
    Interface();

    std::string info() const;
    std::string info_all() const;
};

/// 边界单元
class MESH::MarkElem {
public:
    int cell_id;
    int cell_type;
    int on_cell_face_id;
    int face_id = -1;

    // 符合 Gambit 格式 控制体编号|控制体类型|控制体上界面编号
    explicit MarkElem(const VecStr &list);
};

/// 边界
class MESH::Mark {
public:
    int id;                     // 边界编号
    int type;                   // 边界类型
    int elem_num;               // 边界元素个数
    std::string name;           // 边界名称
    std::string type_name;      // 边界类型名称
    MeshMlEm MARK_ELEM;         // 边界元素容器

    // 边界物理参数
    double density = 0., temperature = 0.;
    Vector velocity = {0., 0., 0.};

    // 构造函数
    explicit Mark(const VecStr &list);

    // 边界信息
    void info() const;
};

/// 网格对象
class MESH::Mesh {
public:
    int status = 0;
    /**
     * 网格类型
     *  0 - 物理空间网格
     *  1 - 速度空间网格
     **/
    int type = -1;
    int NUMNP = 0;  // number of node points
    int NELEM = 0;  // number of elements
    int NGRPS = 0;  // number of element groups
    int NBSETS = 0; // number of B.C. sets
    int NDFCD = 0;  // number of coordinate directions
    int NDFVL = 0;  // number of velocity components

    std::string name;
    MeshNode NODES;
    MeshFace FACES;
    MeshCell CELLS, SHADOW_CELLS;
    MeshMark MARKS;

    double max_discrete_velocity = -1.0;
    double min_mesh_size = -1.0;

    Mesh() = default;
    Mesh(int mesh_type, std::string mesh_name) : name(std::move(mesh_name)), type(mesh_type) {};
    void parse_meshInfo(const VecStr &vector);
    void load(const std::string &file_path);
    int dimension() const;
    void info();
    void update_max_discrete_velocity(double value);
    void update_min_mesh_size(double value);
    void shrink_to_fit();
    void set_mark_params(const BoundaryParam &bc_param);

    /// 网格几何构建  MeshBuild.cpp
    void BuildInterface();        // 构建界面
    void BuildCellGeom();         // 构建控制体格心坐标、体积
    void BuildFaceGeom();         // 构建界面格心坐标、面积
    void BuildNormalVector();     // 构建法向量
    void BuildMesh();             // 完全构建网格
};


/**
 * 算例管理
 **/

class ConfigReader : Reader {
protected:
    const std::string end_mark = "END";
    const std::string def_mark = "DEF";
    const std::string def_end_mark = "END_OF_DEF";
public:
    using ParamMap = std::unordered_map<std::string, std::string>;
    ParamMap case_info = {
            {"CASE_NAME", "NULL"},
            {"PHY_MESH", "NULL"},
            {"DVS_MESH", "NULL"},
            {"SOLVER", "NULL"},
            {"THREAD_NUM", "1"}
    };
    ParamMap param;
    Vec(BoundaryParam) marks;
    explicit ConfigReader(const std::string &case_name);
    void info();
    std::string operator[](const std::string &key);
};


/// 全局函数
void OutputMesh(const std::string &path, MESH::Mesh &mesh);
