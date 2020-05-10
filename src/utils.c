#include "utils.h"


enum VERBOSITY_LEVELS VERBOSITY = QUIET;


json_t* read_json_file(const char *path)
{
    json_error_t j_error;
    json_t* root = json_load_file(path, 0, &j_error);
    if (root == NULL) {
        JSON_ERR_MSG(j_error, "Can't load file %s", path);
        return NULL;
    }
    return root;
}
