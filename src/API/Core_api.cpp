#include <API.h>


const std::string prefix = "API";


void api_init() {
    meso_init();
    pprint::note << "Initialized from api.";
    pprint::note(prefix);
}

void api_create_case(const char *name) {
    std::stringstream ss;
    ss << "./case/" << name;
    if (create_dir(ss.str())) {
        create_dir(ss.str() + "/result");
        create_dir(ss.str() + "/cache");
        output_default_config(name);
        pprint::note << "Create default case: " << name;
        pprint::note(prefix);
    } else {
        pprint::error << "Cannot create case: " << name;
        pprint::error(prefix);
    }
}

void api_set_core_params(bool api_debug_mode, int api_save_interval,
                         int api_residual_interval, int api_max_step,
                         double api_residual_limit) {
    debug_mode = api_debug_mode;
    save_interval = api_save_interval;
    residual_interval = api_residual_interval;
    max_step = api_max_step;
    residual_limit = api_residual_limit;
}

int api_parse_mesh(const char *path) {
    try {
        MESH::Mesh mesh(0, "parsed_mesh");
        mesh.load(path);
        mesh.BuildMesh();
        mesh.info();
        /// output mesh info
        if (debug_mode) {
            pprint::debug << "Interface info:";
            pprint::debug(prefix);
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
        pprint::error(prefix);
        return -1;
    }
}

ConfigReader *api_new_config_reader(const char *path) {
    auto reader = new ConfigReader(path);
    reader->info();
    return reader;
}

const char *api_config_get_value(ConfigReader *reader, const char *key) {
    const char *p;
    p = strdup ((*reader)[key].c_str());    // cstring - deep copy
    if (debug_mode) {
        pprint::debug << "api_config_get_value() return:" << p;
        pprint::debug(prefix);
    }
    return p;
}
