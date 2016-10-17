
#define LOG_ERROR(fmt, ...)                                 \
    do {                                                    \
        fprintf(stderr, "%s:%d:%s(): " fmt "\n",            \
            __FILE__, __LINE__, __func__, ##__VA_ARGS__);   \
    } while (0)

#define MIN(a, b)                   \
    ({                              \
        __typeof__ (a) _a = (a);    \
        __typeof__ (b) _b = (b);    \
        _a < _b ? _a : _b;          \
    })

enum
{
    RANK_MASTER     = 0,
};
