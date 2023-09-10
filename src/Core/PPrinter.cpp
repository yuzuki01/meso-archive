#include <Core.h>

std::string c_none = "\033[0m";
std::string c_note = "\033[1;36m";
std::string c_warn = "\033[1;31m";
std::string c_error = "\033[41;37m";
std::string c_debug = "\033[1;35m";
std::string c_green = "\033[1;32m";
std::string c_highlight = "\033[1;33m";

void pprint::init() {
#if OS == 0     // windows
    system("cls");
#elif OS == 1   // linux
    if (system("clear") == -1) {
        std::cout << "pprint::init() caught an error." << std::endl;
    }
#endif
};


///pprint
pprint::PPrinter::PPrinter(const int output_lv) {
    switch (output_lv) {
        case 1:
            color_mark = c_note;
            break;
        case 2:
            color_mark = c_warn;
            break;
        case 3:
            color_mark = c_error;
            break;
        case 4:
            color_mark = c_debug;
            break;
        case 5:
            color_mark = c_highlight;
            break;
        case 0:
        default:
            color_mark = c_none;
            break;
    }
}

pprint::PPrinter &pprint::PPrinter::operator<<(const int value) {
    ss << value;
    return *this;
}

pprint::PPrinter &pprint::PPrinter::operator<<(const double value) {
    ss << std::setprecision(OUT_PRECISION) << value;
    return *this;
}

pprint::PPrinter &pprint::PPrinter::operator<<(const std::string &message) {
    ss << message;
    return *this;
}

void pprint::PPrinter::operator()() {
    std::cout << color_mark << ss.str() << c_none << std::endl;
    ss.str("");
}

void pprint::PPrinter::operator()(const std::string &prefix) {
    std::cout << color_mark << "[" + prefix + "] " << ss.str() << c_none << std::endl;
    ss.str("");
}

void pprint::process_bar(int step, int interval) {
    int progress = step % interval;
    if (progress == 0) progress = interval;
    if (progress == 1) std::cout << "Process: \n";
    float percentage = static_cast<float>(progress) / interval;
    int progressBarWidth = 50;
    int progressBarPosition = static_cast<int>(percentage * progressBarWidth);

    std::cout << "     [";
    for (int i = 0; i < progressBarWidth; ++i) {
        if (i < progressBarPosition) {
            std::cout << "=";
        } else {
            std::cout << " ";
        }
    }
    std::cout << "] " << progress << "/" << interval << " " << static_cast<int>(percentage * 100) << "%" << "\r";
    std::cout.flush();
    if (progress == interval) std::cout << "\n";
}

///pprint objects
pprint::PPrinter pprint::info(0);
pprint::PPrinter pprint::note(1);
pprint::PPrinter pprint::warn(2);
pprint::PPrinter pprint::error(3);
pprint::PPrinter pprint::debug(4);
pprint::PPrinter pprint::highlight(5);
