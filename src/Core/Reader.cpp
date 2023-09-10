#include <Core.h>

Reader::Reader(const std::string &name, const std::string &file_path) {
    reader_name = name;
    std::ifstream input_file(file_path);
    if (!input_file.is_open()) {
        pprint::error << "Error: Cannot open the file " << file_path;
        pprint::error(reader_name);
        throw std::invalid_argument("Reader cannot open file.");
    }
    std::string line;
    while (std::getline(input_file, line)) {
        if (!line.empty()) lines.push_back(line);
    }

    pprint::info << "Open file: " << file_path;
    pprint::info(reader_name);
}
