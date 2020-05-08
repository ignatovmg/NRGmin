#include "utils.h"


enum VERBOSITY_LEVELS VERBOSITY = QUIET;


json_t* read_json_file(const char *path)
{
    json_error_t error;
    json_t* root = json_load_file(path, 0, &error);
    if (root == NULL) {
        ERR_MSG("Can't load file %s", path);
        return NULL;
    }
    return root;
}
