#include "tng_io.h"

/*run: gcc TestTNG.c -Wall -I ${TNG_INCLUDE_DIR} -ltng_io && ./a.out ; echo $?*/
int main()
{
//     int v;
    tng_trajectory_t traj = 0;
    char version_str[TNG_MAX_STR_LEN];

    tng_version(traj, version_str, TNG_MAX_STR_LEN);

    fprintf(stderr, "%s", version_str);

//     if (argc >= 2)
//     {
//         tng_version_major(traj, &v);
//         if(v < argv[1])
//         {
//             return 1;
//         }
//     }
//     if (argc >= 3)
//     {
//         tng_version_minor(traj, &v);
//         if(v < argv[2])
//         {
//             return 1;
//         }
//     }
//     if (argc >= 4)
//     {
//         tng_version_patchlevel(traj, &v);
//         if(v < argv[3])
//         {
//             return 1;
//         }
//     }
    return 0;
}
