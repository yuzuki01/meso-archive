#include <Core.h>

bool is_core_init = false;
bool debug_mode = false;
int save_interval = 1000;
int residual_interval = 1000;
int max_step = 1000000;
double residual_limit = 1e-6;


void meso_init() {
    if (is_core_init) return;
    // init pprint
    pprint::init();
    // create ./case
    if (!create_dir("./case")) {
        pprint::error << "Cannot create dir: ./case";
        pprint::error("Init");
        exit(-1);
    }
    is_core_init = true;
}

bool create_dir(const std::string &path) {
#if OS == 1   // linux
    if (access(path.c_str(), 0) == -1) mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    return (access(path.c_str(), 6) == 0);
#elif OS == 0 // windows
    if (_access(path.c_str(), 0) == -1) _mkdir(path.c_str());
    return (_access(path.c_str(), 6) == 0);
#else
    return false;
#endif
}

VecStr split(const std::string &_str) {
    std::istringstream iss(_str);
    std::vector<std::string> data;
    std::string token;

    while (iss >> token) {
        data.push_back(token);
    }

    return data;
}

bool str_vec_cmp(const VecStr &data, const VecStr &reference) {
    int i = data.size(), j = reference.size();
    if (i < j) return false;
    for (i = 0; i < j; i++) {
        if (data[i] != reference[i]) return false;
    }
    return true;
}

bool is_confirmed(const std::string &context) {
    std::cout << std::endl << context << " [y/n]" << std::endl;
    std::string input;
    std::cin >> input;
    if (input == "y") return true;
    return false;
}

std::string arg_parser_param(int argc, char **argv, const std::string &target, const std::string &_default) {
    for (int i = 0; i < argc; i++) {
        if (argv[i] == target) {
            if (i + 1 < argc) {
                return std::string(argv[i + 1]);
            } else return _default;
        }
    }
    return _default;
}

bool arg_parser_switch(int argc, char **argv, const std::string &target) {
    for (int i = 0; i < argc; i++) {
        if (argv[i] == target) return true;
    }
    return false;
}

std::string generate_key_from_node_set(const VecInt &vector) {
    if (vector.empty()) return "null";
    // 找出最小值
    int min_value = *std::min_element(vector.begin(), vector.end());
    // 找出最小值索引
    auto min_value_pos = std::find(vector.begin(), vector.end(), min_value);
    // 生成新容器
    VecInt result;
    // 预分配空间
    result.reserve(vector.size());
    // 将 vector 中最小值及其之后的元素依次添加到 result 中
    for (auto it = min_value_pos; it != vector.end(); ++it) {
        result.push_back(*it);
    }
    // 将 vector 中最小值之前的元素依次添加到 result 中
    for (auto it = vector.begin(); it != min_value_pos; ++it) {
        result.push_back(*it);
    }
    // 成环方向
    int len = result.size();
    if (len > 3)
        if (result[1] > result.back()) {
            VecInt tmp;
            tmp.reserve(result.size());
            tmp.push_back(result.front());
            for (int i = len - 1; i > 0; i--) {
                tmp.push_back(result[i]);
            }
            result = tmp;
        }
    // 输出为-隔开的字符串
    std::stringstream ss;
    for (size_t i = 0; i < result.size(); i++) {
        if (i != 0) ss << "-";
        ss << result[i];
    }
    return ss.str();
}

std::string output_data_to_console(const VecStr &name, const VecDouble &value) {
    if (name.size() != value.size()) {
        pprint::warn << "Caught lists with different sizes to output.";
        pprint::warn("Core");
        return "";
    }
    const int width = 15;
    std::stringstream ssn, ssv;
    ssn << std::right;
    for (auto &var : name) {
        ssn << std::setw(width) << var;
    }
    ssn << std::endl;
    ssv << std::right;
    for (auto &var : value) {
        ssv << std::scientific << std::setprecision(OUT_PRECISION) << std::setw(width) << var;
    }
    return ssn.str() + ssv.str();
}

void output_default_config(const std::string &name) {
    std::stringstream ss;
    ss << "./case/" << name << "/config.txt";
    std::ofstream fp;
    fp.open(ss.str(), std::ios::out | std::ios::trunc);
    // check
    if (!fp.is_open()) {
        pprint::error << "Cannot write to file: " << ss.str();
        pprint::error("Core");
        throw std::invalid_argument("Cannot write to file.");
    } else {
        fp << "% ####################\n";
        fp << "% # 这是一个示例文件 #\n";
        fp << "% ####################\n";
        fp << "\n";
        fp << "CASE-INFO\n";
        fp << "\n";
        fp << "% 算例名称 <必须>\n";
        fp << " CASE_NAME       demo\n";
        fp << "% 物理空间网格 <必须>\n";
        fp << " PHY_MESH        demo.neu\n";
        fp << "% 速度空间网格 <必须>\n";
        fp << " DVS_MESH        NULL\n";
        fp << "% 求解器类型 <必须>\n";
        fp << " SOLVER          solver\n";
        fp << "% 并行线程数 <必须>\n";
        fp << " THREAD_NUM      10\n";
        fp << "\n";
        fp << "END\n";
        fp << "\n";
        fp << "PARAM-INFO\n";
        fp << "% 流体力学参数 <根据求解器选择变量>\n";
        fp << "% dugks@incompressible    requires    Re, Ma\n";
        fp << "\n";
        fp << " Re             -1\n";
        fp << " Ma             -1\n";
        fp << " Kn             -1\n";
        fp << "\n";
        fp << "% CFL 数 <必须>\n";
        fp << " CFL            0.9\n";
        fp << "% 无量纲气体常数 <必须>\n";
        fp << " R              0.5\n";
        fp << "% 无量纲特征长度 <必须>\n";
        fp << " LENGTH         1.0\n";
        fp << "% 无量纲参考密度 <必须>\n";
        fp << " DENSITY        1.0\n";
        fp << "% 无量纲参考温度 <必须>\n";
        fp << " TEMPERATURE    1.0\n";
        fp << "END\n";
        fp << "\n";
        fp << "BOUNDARY-INFO\n";
        fp << "% 定义边界信息文本块\n";
        fp << "\n";
        fp << "DEF\n";
        fp << "    % 边界名称 <必须>\n";
        fp << "    NAME        MarkName\n";
        fp << "    % 边界类型 <必须>\n";
        fp << "    TYPE        isothermal_wall\n";
        fp << "    % 边界密度 <非必须>\n";
        fp << "    % DENSITY   1.0\n";
        fp << "    % 边界温度 <非必须>\n";
        fp << "    TEMPERATURE 1.0\n";
        fp << "    % 边界速度分量 <非必须, 默认 velocity=<0,0,0>>\n";
        fp << "    VELOCITY-X      0.0\n";
        fp << "    % VELOCITY-Y    0.0\n";
        fp << "    % VELOCITY-Z    0.0\n";
        fp << "END_OF_DEF\n";
        fp << "\n";
        fp << "END";
    }
    fp.close();
}

void clear_log(const std::string &case_name) {
    std::stringstream ss;
    ss << "./case/" << case_name << "/cache/log.txt";
    std::ofstream fp;
    fp.open(ss.str(), std::ios::out | std::ios::trunc);
    // check
    if (!fp.is_open()) {
        pprint::error << "Cannot write to file: " << ss.str();
        pprint::error("Core");
        throw std::invalid_argument("Cannot write to file.");
    }
    fp.close();
}

void output_log(const std::string &case_name, const std::string &log_context, const std::string &prefix) {
    std::stringstream ss;
    ss << "./case/" << case_name << "/cache/log.txt";
    std::ofstream fp;
    fp.open(ss.str(), std::ios::out);
    // check
    if (!fp.is_open()) {
        pprint::error << "Cannot write to file: " << ss.str();
        pprint::error("Core");
        throw std::invalid_argument("Cannot write to file.");
    } else {
        fp << "[" << prefix << "] " << log_context << std::endl;
    }
    fp.close();
}
