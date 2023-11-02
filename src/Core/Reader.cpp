#include <Core.h>

Reader::Reader(const std::string &name, const std::string &file_path) {
    reader_name = name;
    std::ifstream input_file(file_path);
    file_opened = input_file.is_open();
    if (!file_opened) {
        pprint::warn << "Error: Cannot open the file " << file_path;
        pprint::warn(reader_name);
        return;
    }
    std::string line;
    while (std::getline(input_file, line)) {
        if (!line.empty()) lines.push_back(line);
    }

    pprint::info << "Open file: " << file_path;
    pprint::info(reader_name);
}

bool Reader::is_file_open() const {
    return file_opened;
}
