#ifndef ENERGYMIN_UTILS_H
#define ENERGYMIN_UTILS_H

#include <stdio.h>

#define ERR_MSG(fmt, ...) do {                                              \
    fprintf(stderr,"[Error] (file %s, %s, line %i):\n" fmt "\n",            \
            __FILE__, __func__, __LINE__, ##__VA_ARGS__);                   \
} while(0)

#define WRN_MSG(fmt, ...) do {                                         \
    fprintf(stderr,"[Warning] (file %s, %s, line %i):\n" fmt "\n",     \
            __FILE__, __func__, __LINE__, ##__VA_ARGS__);              \
} while(0)

#define INFO_MSG(fmt, ...) do {                               \
    fprintf(stderr,"[Info]: " fmt "\n", ##__VA_ARGS__);            \
} while(0)

#define FOPEN_ELSE(handle, path, mode) \
    handle = fopen((path), (mode)); \
    if ((handle) == NULL) { \
        ERR_MSG("Can't open %s", (path)); \
    } \
    if ((handle) == NULL) \

#endif //ENERGYMIN_UTILS_H
