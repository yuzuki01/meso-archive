#ifndef INCLUDE_CORE
#define INCLUDE_CORE

/**
 * 核心库
 **/
#include <ctime>
#include <cmath>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <exception>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <utility>
#include <omp.h>

#if defined(_WIN32)
/// Windows
#define OS 0

#include <direct.h>

#elif defined(__linux__)
/// Linux
#define OS 1
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#else
/// Unsupported platform
#define OS -1

#endif


/**
 * 宏定义
 */

#define Pi                3.14159265358979323846
#define Pi_over_2         1.57079632679489661923
#define LINE_DATA_NUM 25
#define DATA_PRECISION 17
#define OUT_PRECISION 6
#define Vec(x) std::vector<x>
#define VecStr Vec(std::string)
#define VecInt Vec(int)
#define VecDouble Vec(double)
#define MapInt std::unordered_map<std::string, int>


/**
 * 类
 */


class Reader {
protected:
    std::string reader_name;
    VecStr lines;
public:
    explicit Reader(const std::string &name, const std::string &file_path);

    ~Reader() = default;
};

namespace pprint {
    void init();

    // output func
    class PPrinter {
    private:
        std::string color_mark;
        std::stringstream ss;
    public:
        explicit PPrinter(int output_lv);

        PPrinter &operator<<(const std::string &message);

        PPrinter &operator<<(int value);

        PPrinter &operator<<(double value);

        void operator()();

        void operator()(const std::string &prefix);
    };

    extern PPrinter info;
    extern PPrinter note;
    extern PPrinter warn;
    extern PPrinter error;
    extern PPrinter debug;
    extern PPrinter highlight;

    void process_bar(int progress, int total);
}


/**
 * 全局函数和变量
 */

extern bool is_core_init;
extern bool debug_mode;
extern int save_interval;
extern int residual_interval;
extern int max_step;
extern double residual_limit;

void meso_init();

bool create_dir(const std::string &path);

VecStr split(const std::string &_str);

bool str_vec_cmp(const VecStr &data, const VecStr &reference);

bool is_confirmed(const std::string &context);

std::string generate_key_from_node_set(const VecInt &vector);

std::string output_data_to_console(const VecStr &name, const VecDouble &value);

void output_default_config(const std::string &case_name);

void clear_log(const std::string &case_name);

void output_log(const std::string &case_name, const std::string &log_context, const std::string &prefix);

/**
 * 参数解析
 */
std::string arg_parser_param(int argc, char **argv, const std::string &target, const std::string &_default);

bool arg_parser_switch(int argc, char **argv, const std::string &target);


/**
 *  自定义头文件
 */

/// NULL

#endif