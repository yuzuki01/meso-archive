/// 最小二乘法
struct LeastSecondParam {
    Vec(Vector) dr;
    VecDouble weight;
    Vector Cx, Cy, Cz;
};

LeastSecondParam GenerateLeastSecondParam(int cell_id, MESH::Mesh &mesh);

/// 残差
class Residual {
protected:
    double now, old;
public:
    explicit Residual(double value) : now(value), old(value) {};
    void update(double value);
    double compute() const;
};

/// 高斯型积分
MESH::Mesh GenerateMesh_GaussHermit(ConfigReader &reader, int mesh_type);
