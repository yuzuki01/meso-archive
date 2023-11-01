#include <Mesh.h>


/// NEU Reader
void MESH::NEUReader::parse(MESH::Mesh &mesh) {
    int read_case = 0;
    for (int i = 0; i < int(lines.size()); i++) {
        auto data = split(lines[i]);
        // 跳过注释
        if (data.empty()) continue;
        switch (read_case) {
            case 0: // find CONTROL INFO
                if (str_vec_cmp(data, split("CONTROL INFO"))) {
                    while (true) {
                        data = split(lines[++i]);
                        if (data.empty()) continue;
                        if (data[0] == end_mark) break;
                        if (str_vec_cmp(split(lines[i]), split("NUMNP NELEM NGRPS NBSETS NDFCD NDFVL"))) {
                            mesh.parse_meshInfo(split(lines[++i]));
                        }
                    }
                    read_case = 1;
                }
                break;
            case 1: // find NODAL COORDINATES
                if (str_vec_cmp(data, split("NODAL COORDINATES"))) {
                    while (true) {
                        data = split(lines[++i]);
                        if (data.empty()) continue;
                        if (data[0] == end_mark) break;
                        mesh.NODES.emplace_back(data);
                    }
                    read_case = 2;
                }
                break;
            case 2: // find ELEMENTS/CELLS
                if (str_vec_cmp(data, split("ELEMENTS/CELLS"))) {
                    while (true) {
                        data = split(lines[++i]);
                        if (data.empty()) continue;
                        if (data[0] == end_mark) break;
                        if (int(data.size()) < 3 + stoi(data[2])) {
                            VecStr next_line = split(lines[++i]);
                            for (auto &it : next_line) data.push_back(it);
                        }
                        mesh.CELLS.emplace_back(data);
                    }
                    read_case = 3;
                }
                break;
            case 3: // find BOUNDARY CONDITIONS
                if (str_vec_cmp(data, split("BOUNDARY CONDITIONS"))) {
                    data = split(lines[++i]);
                    if (data.empty()) continue;
                    int mark_id = mesh.MARKS.size();
                    mesh.MARKS.emplace_back(data);
                    auto &mark = mesh.MARKS.back();
                    mark.id = mark_id;
                    while (true) {
                        data = split(lines[++i]);
                        if (data.empty()) continue;
                        if (data[0] == end_mark) break;
                        mark.MARK_ELEM.emplace_back(split(lines[i]));
                    }
                }
                break;
            default:
                break;
        }
    }
    mesh.status = 1;
    pprint::info << "Parsed file: " << mesh.name;
    pprint::info(reader_name);
}


/// Mesh Writer
void OutputMesh(const std::string &path, MESH::Mesh &mesh) {
    int count;
    std::ofstream fp;
    // open file
    fp.open(path, std::ios::out | std::ios::trunc);
    if (!fp.is_open()) {
        pprint::error << "Cannot write file: " << path;
        pprint::error("Case");
        return;
    }
    // write file
    fp << GEOM::tecplot_file_header(mesh, {"volume"});
    /// X
    count = 0;
    fp << std::endl << "# X" << std::endl;
    for (auto &node : mesh.NODES) {
        fp << " " << node.pos.x;
        count++;
        if (count == LINE_DATA_NUM) {
            fp << std::endl;
            count = 0;
        }
    }
    /// Y
    count = 0;
    fp << std::endl << "# Y" << std::endl;
    for (auto &node : mesh.NODES) {
        fp << " " << node.pos.y;
        count++;
        if (count == LINE_DATA_NUM) {
            fp << std::endl;
            count = 0;
        }
    }
    /// Z
    if (mesh.dimension() == 3) {
        count = 0;
        fp << std::endl << "# Z" << std::endl;
        for (auto &node : mesh.NODES) {
            fp << " " << node.pos.z;
            count++;
            if (count == LINE_DATA_NUM) {
                fp << std::endl;
                count = 0;
            }
        }
    }
    /// volume
    count = 0;
    fp << std::endl << "# volume" << std::endl;
    for (auto &cell : mesh.CELLS) {
        fp << " " << cell.volume;
        count++;
        if (count == LINE_DATA_NUM) {
            fp << std::endl;
            count = 0;
        }
    }
    /// Geom
    fp << std::endl << "# Geom" << std::endl;
    for (auto &cell : mesh.CELLS) {
        fp << GEOM::tecplot_cell_format(cell);
    }
    // close file
    fp.close();

    pprint::note << "Write mesh file to: " << path;
    pprint::note("Mesh");
}


/// Config Reader
ConfigReader::ConfigReader(const std::string &case_name) : Reader("ConfigReader",
                                                                  "./case/" + case_name + "/config.txt") {
    int read_case = 0;
    Vec(Vec(VecStr)) bc_info;
    for (int i = 0; i < int(lines.size()); i++) {
        auto data = split(lines[i]);
        // 跳过注释
        if (data.empty()) continue;
        if (data[0] == "%") continue;
        switch (read_case) {
            case 0: // find CASE-INFO
                if (data[0] == "CASE-INFO") {
                    while (true) {
                        data = split(lines[++i]);
                        if (data.empty()) continue;
                        if (data[0] == "%") continue;
                        if (data[0] == end_mark) break;
                        auto it = case_info.find(data[0]);
                        if (it != case_info.end()) {
                            case_info[data[0]] = data[1];
                        }
                    }
                    read_case = 1;
                }
                break;
            case 1: // find PARAM-INFO
                if (data[0] == "PARAM-INFO") {
                    while (true) {
                        data = split(lines[++i]);
                        if (data.empty()) continue;
                        if (data[0] == "%") continue;
                        if (data[0] == end_mark) break;
                        if (data.empty() || data[0] == "%") continue;
                        param[data[0]] = data[1];
                    }
                    read_case = 2;
                }
                break;
            case 2: // find BOUNDARY-INFO
                if (data[0] == "BOUNDARY-INFO") {
                    while (true) {
                        data = split(lines[++i]);
                        if (data.empty()) continue;
                        if (data[0] == "%") continue;
                        if (data[0] == end_mark) break;
                        if (data[0] == def_mark) {
                            Vec(VecStr) bc_item;
                            while (true) {
                                data = split(lines[++i]);
                                if (data.empty()) continue;
                                if (data[0] == "%") continue;
                                if (data[0] == def_end_mark) break;
                                bc_item.push_back(data);
                            }
                            bc_info.push_back(bc_item);
                        }
                    }
                    read_case = 3;
                }
            default:
                break;
        }
    }
    // bc
    for (auto &bc_it : bc_info) {
        BoundaryParam bc;
        for (auto &it : bc_it) {
            if (it.size() < 2) continue;
            if (it[0] == "NAME") {
                bc.name = it[1];
                continue;
            }
            if (it[0] == "TYPE") {
                auto p = MESH::MarkType.find(it[1]);
                if (p == MESH::MarkType.end()) {
                    pprint::error << "Config Reader caught unsupported boundary type.";
                    pprint::error("ConfigReader");
                    throw std::invalid_argument("Config Reader caught unsupported boundary type.");
                }
                bc.type = MESH::MarkType[it[1]];
                bc.type_name = it[1];
                continue;
            }
            if (it[0] == "DENSITY") {
                bc.density = stod(it[1]);
                continue;
            }
            if (it[0] == "TEMPERATURE") {
                bc.temperature = stod(it[1]);
                continue;
            }
            if (it[0] == "VELOCITY-X") {
                bc.velocity.x = stod(it[1]);
                continue;
            }
            if (it[0] == "VELOCITY-Y") {
                bc.velocity.y = stod(it[1]);
                continue;
            }
            if (it[0] == "VELOCITY-Z") {
                bc.velocity.z = stod(it[1]);
                continue;
            }
        }
        marks.push_back(bc);
    }

    pprint::note << "Parse config file(case=" << case_info["CASE_NAME"] << ") successfully.";
    pprint::note("ConfigReader");
}


void ConfigReader::info() {
    pprint::note << "    Config info:";
    pprint::note("Config Reader");
    pprint::info << "  case name: " << case_info["CASE_NAME"] << "\n"
                 << "  phy mesh:  " << case_info["PHY_MESH"] << "\n"
                 << "  dvs mesh:  " << case_info["DVS_MESH"] << "\n"
                 << "  solver:    " << case_info["SOLVER"] << "\n"
                 << "  thread:    " << case_info["THREAD_NUM"];
    pprint::info();
}

std::string ConfigReader::operator[](const std::string &key) {
    /// find in case-info
    {
        auto it = case_info.find(key);
        if (it != case_info.end()) return case_info[key];
    }
    {
        auto it = param.find(key);
        if (it != param.end()) {
            return param[key];
        } else {
            pprint::warn << "Cannot find var:" << key << " in config context, return \"NULL\".";
            pprint::warn("Config Reader");
            return MeshReaderReturn_NULL;
        }
    }
}
