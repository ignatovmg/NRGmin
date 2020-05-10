#ifndef ENERGYMIN_UTILS_H
#define ENERGYMIN_UTILS_H

#include <stdio.h>
#include <jansson.h>


enum VERBOSITY_LEVELS {
    QUIET = 0,
    ERROR = 1,
    WARNING = 2,
    INFO = 3,
    DEBUG = 4
};


extern enum VERBOSITY_LEVELS VERBOSITY;


#define JSON_ERR_MSG(json_error, fmt, ...) do {                         \
    if (VERBOSITY > QUIET) {                                            \
        fprintf(stderr,                                                 \
            "[Json parse error] (%s:%i, %s):\n\n"                       \
            "                    " fmt "\n\n"                           \
            "                    Source: %s\n"                          \
            "                    Message: %s\n"                         \
            "                    Line: %i\n"                            \
            "                    Column: %i\n\n",                       \
            __FILE__, __LINE__, __func__,                               \
            ##__VA_ARGS__,                                              \
            json_error.source,                                          \
            json_error.text,                                            \
            json_error.line,                                            \
            json_error.column);                                         \
    }                                                                   \
} while(0)

#define ERR_MSG(fmt, ...) do {                                          \
    if (VERBOSITY > QUIET) {                                            \
        fprintf(stderr,                                                 \
            "[Error] (%s:%i, %s):\n\n"                                  \
            "         " fmt "\n\n",                                     \
            __FILE__, __LINE__, __func__, ##__VA_ARGS__);               \
    }                                                                   \
} while(0)

#define WRN_MSG(fmt, ...) do {                                          \
    if (VERBOSITY > ERROR) {                                            \
        fprintf(stderr,                                                 \
             "[Warning] (%s:%i, %s):\n\n"                               \
             "           " fmt "\n\n",                                  \
             __FILE__, __LINE__, __func__, ##__VA_ARGS__);              \
    }                                                                   \
} while(0)

#define INFO_MSG(fmt, ...) do {                                         \
    if (VERBOSITY > WARNING) {                                          \
        fprintf(stderr, "[Info]:  " fmt "\n", ##__VA_ARGS__);           \
    }                                                                   \
} while(0)

#define DEBUG_MSG(fmt, ...) do {                                        \
    if (VERBOSITY > INFO) {                                             \
        fprintf(stderr,"[Debug]: " fmt "\n", ##__VA_ARGS__);            \
    }                                                                   \
} while(0)

#define FOPEN_ELSE(handle, path, mode)                                  \
    handle = fopen((path), (mode));                                     \
    if ((handle) == NULL) {                                             \
        ERR_MSG("Can't open %s", (path));                               \
    }                                                                   \
    if ((handle) == NULL)


json_t* read_json_file(const char *path);


#endif //ENERGYMIN_UTILS_H