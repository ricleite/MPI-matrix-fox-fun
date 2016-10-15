
#define LOG_ERROR(fmt, ...)                                 \
    do {                                                    \
        fprintf(stderr, "%s:%d:%s(): " fmt,                 \
            __FILE__, __LINE__, __func__, ##__VA_ARGS__);   \
    } while (0)


enum
{
    RANK_MASTER     = 0,
};
