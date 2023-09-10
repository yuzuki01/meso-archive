#include "main.h"


int main(int argc, char **argv) {
    // init
    meso_init();
    // arg parse
    {
        if (arg_parser_switch(argc, argv, "-debug")) {
            try {
                handle_debug();
            } catch (std::exception &ex) {
                pprint::error << "Caught error: " << ex.what();
                pprint::error("Main");
                exit(-1);
            }
        }
        // help
        if (arg_parser_switch(argc, argv, "-h")) {
            return handle_help();
        }
        if (arg_parser_switch(argc, argv, "-help")) {
            return handle_help();
        }
        std::string parsed_param;
        // create case
        parsed_param = arg_parser_param(argc, argv, "--create", "");
        if (!parsed_param.empty()) return handle_create_case(parsed_param);
        // max_step
        parsed_param = arg_parser_param(argc, argv, "--max_step", "");
        if (!parsed_param.empty()) max_step = stoi(parsed_param);
        // save_interval
        parsed_param = arg_parser_param(argc, argv, "--save_interval", "");
        if (!parsed_param.empty()) save_interval = stoi(parsed_param);
        // residual_interval
        parsed_param = arg_parser_param(argc, argv, "--residual_interval", "");
        if (!parsed_param.empty()) residual_interval = stoi(parsed_param);
        // residual_limit
        parsed_param = arg_parser_param(argc, argv, "--esp", "");
        if (!parsed_param.empty()) residual_limit = stod(parsed_param);
        // run case
        parsed_param = arg_parser_param(argc, argv, "--case", "");
        if (!parsed_param.empty()) {
            return handle_case(parsed_param);
        }
        // parse a mesh
        parsed_param = arg_parser_param(argc, argv, "--parse_mesh", "");
        if (!parsed_param.empty()) {
            return handle_mesh(parsed_param);
        }
    }

    // arg parse failed
    if (debug_mode) {
        return 1;
    } else {
        return handle_parse_failed();
    }
}

int handle_parse_failed() {
    pprint::info << R"(Use command "-h" or "-help" to get help info.)";
    pprint::info();
    return 0;
}

int handle_case(const std::string &path) {
    // 实例化算例对象
    ConfigReader reader(path);
    reader.info();
    std::string solver_name = reader["SOLVER"];
    if (solver_name == "dugks@incompressible") {
        return handle_solver<DUGKS_INCOMPRESSIBLE>(reader);
    } else if (solver_name == "dugks@aoki") {
        return handle_solver<DUGKS_AOKI>(reader);
    }
    return 0;
}

int handle_mesh(const std::string &path) {
    try {
        MESH::Mesh mesh(0, "parsed_mesh");
        mesh.load(path);
        mesh.BuildMesh();
        mesh.info();
        /// output mesh info
        if (debug_mode) {
            pprint::debug << "Interface info:";
            pprint::debug("Main");
            int count = 0;
            for (auto &face : mesh.FACES) {
                pprint::debug << "   " << face.info_all();
                pprint::debug();
                if (++count >= 100) {
                    pprint::debug << "   ...... " << int(mesh.FACES.size() - count) << " objects left.";
                    pprint::debug();
                    break;
                }
            }
        }
        OutputMesh("parsed_mesh.dat", mesh);
        return 0;
    } catch (std::exception &ex) {
        pprint::error << "Caught error when handling parsing mesh: " << ex.what();
        pprint::error("Main");
        return -1;
    }
}

int handle_create_case(const std::string &name) {
    std::stringstream ss;
    ss << "./case/" << name;
    if (create_dir(ss.str())) {
        create_dir(ss.str() + "/result");
        create_dir(ss.str() + "/cache");
        output_default_config(name);
        pprint::note << "Create default case: " << name;
        pprint::note("Main");
        return 0;
    } else {
        pprint::error << "Cannot create case: " << name;
        pprint::error("Main");
        return -1;
    }
}

int handle_help() {
    std::cout << "Mesoscopic Kinetic Solver\n"
                 "  meso [<command>]\n"
                 "    -h / -help                      Get help info\n"
                 "    --case <NAME>                   Run case <NAME>\n"
                 "    --create <NAME>                 Create a case named <NAME> with \n"
                 "                                    a default config file\n"
                 "    --max_step <VALUE>              Set max-step to <VALUE>,\n"
                 "                                    default = 1000000\n"
                 "    --save_interval <VALUE>         Set save-interval to <VALUE>,\n"
                 "                                    default = 1000\n"
                 "    --residual_interval <VALUE>     Set residual-interval to <VALUE>,\n"
                 "                                    default = 1000\n"
                 "    --esp <VALUE>                   Set residual-limit to <VALUE>,\n"
                 "                                    default = 1e-6\n"
                 "    --parse_mesh <FILE>             Parse mesh file <FILE>, output\n"
                 "                                    tecplot file as parsed_mesh.dat\n";
    return 0;
}

void handle_debug() {
    pprint::warn << "Enter debug-mode.";
    pprint::warn("DEBUG");
    debug_mode = true;
}
