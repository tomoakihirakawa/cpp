//
// create cos table
//
// (c)Copyright Spacesoft corp., 2018 All rights reserved.
//                                    Kitayama, Hiroyuki
//
#include <stdio.h>
#include <math.h>
#include <time.h>

#define PI  3.14159265358979323846

//--------------------------------------------------------------------
// main
int
main(int argc, char *argv[])
{
    double *tbl;
    int size = 4096, i, x, y, centerX, centerY;

    clock_t start, stop;

    if (argc > 1)
    {
        size = atoi(argv[1]);
    }

    tbl = (double *)malloc(sizeof(double) * size * size);

    centerX = centerY = size / 2;



    start = clock();

    double radius = sqrt(pow(centerX, 2) + pow(centerY, 2));

    const int nQueue = 4;
    const int nBloking = 64;
    int blockSize = size / nBloking;
    int blockY;

    if ((size % nBloking) != 0)
    {
        fprintf(stderr, "error: size = %d, nBloking = %d,blockSize = %d, (size %% nBloking) = %d\n",
            size, nBloking, blockSize, (size % nBloking));
        free(tbl);
        return -1;
    }

    #pragma acc data create(tbl[:size*size])
    for (blockY = 0; blockY < nBloking; blockY++)
    {
        int startY = blockY * blockSize;
        int endY = startY + blockSize;

        #pragma acc parallel loop collapse(2) gang vector async(blockY % nQueue +1)
        for (y = startY; y < endY; y++) for (x = 0; x < size; x++)
        {
            // distance from center
            double distance = sqrt(pow(centerY - y, 2) + pow(centerX - x, 2));
            // radius=ƒÎ, current radian
            double radian = (distance / radius) * (double)PI;
            // cosƒÆ, normalize -1.0`1.0 to  0`1.0
            double Y = (cos(radian) + 1.0) / 2.0;
            // normalize (Y) 0`1.0 to 0.0`255.0
            tbl[y*size+x] = Y * 255.0f;
        }
        #pragma acc update self(tbl[startY*size:blockSize*size]) async(blockY % nQueue +1)
    }
    #pragma acc wait

    stop = clock();



    fprintf(stderr, "elapsed time = %.20f [sec]\n",
        (float)(stop - start) / CLOCKS_PER_SEC);


    // print result
    if (argc < 3)
    {
        fprintf(stdout, "%d %d 1\n", size, size);
        for (y = 0; y < size; y++)
        {
            for (x = 0; x < size; x++)
            {
                fprintf(stdout, "%3d\n", (int)tbl[y*size+x]);
            }
        }
    }

    free(tbl);

    return 0;
}
