//
// negative
//
// (c)Copyright Spacesoft corp., 2018 All rights reserved.
//                               Hiro KITAYAMA
//
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

//-------------------------------------------------------------------
// read data
float*
readImgData(char* fname, int* W, int* H)
{
    FILE *fp;
    int width, height, ch;
    float data;

    if ((fp = fopen(fname, "rt")) == NULL)
    {
        fprintf(stderr, "faild file open %s\n", fname);
        return NULL;
    }

    if (fscanf(fp, "%d %d %d", &width, &height, &ch) != 3)
    {
        fprintf(stderr, "failed read file %s\n", fname);
        return NULL;
    }
    *W = width;
    *H = height;

    if (ch != 1)
    {
        fprintf(stderr, "ch != 1.");
        return NULL;
    }
    unsigned int datasize = sizeof(float)*width*height;
    float * in = (float *)malloc(datasize);

    for (int y = 0; y < height; y++)
    {
        for (int x = 0; x < width; x++)
        {
            if (fscanf(fp, "%f", &data) != 1)
            {
                fprintf(stderr, "failed read data %s\n", fname);
                return NULL;
            }
            in[y*width + x] = data;
        }
    }

    fclose(fp);

    return in;
}

//-------------------------------------------------------------------
// effect
float*
effect(const float* in, const int width, const int height,
                                const int nQueue, const int nBloking)
{
    int y, x;
    clock_t start, stop;

    int size = height * width;
    unsigned int datasize = sizeof(float)*size;
    float* out = (float *)malloc(datasize);
    int blockSize = height / nBloking;
    int blockY;

    start = clock();

    #pragma acc data create(in[:size], out[:size])
    for (blockY = 0; blockY < nBloking; blockY++)
    {
        int startY = blockY * blockSize;
        int endY = startY + blockSize;
        int hAsync = (blockY % nQueue) + 1;
        #pragma acc update device(in[startY*width:blockSize*width]) /*async(hAsync)*/
        #pragma acc parallel loop collapse(2) gang vector /*async(hAsync)*/
        for (y = startY; y < endY; y++) for (x = 0; x < width; x++)
        {
            out[y*width + x] = 255.0f - in[y*width + x];
        }
        #pragma acc update self(out[startY*width:blockSize*width]) /*async(hAsync)*/
    }
    /*#pragma acc wait*/

    stop = clock();

    fprintf(stderr, "elapsed time = %.20f [sec]\n",
        (float)(stop - start) / CLOCKS_PER_SEC);

    return out;
}

//-------------------------------------------------------------------
// main
int
main(int argc, char* argv[])
{
    int width, height;
    float* in = NULL;

    if (argc < 2)
    {
        fprintf(stderr, "no <input>");
        return -1;
    }

    if ((in = readImgData(argv[1], &width, &height)) == NULL)
    {
        return -1;
    }

    int nQueue = 4;
    if (argc > 2)
    {
        nQueue = atoi(argv[2]);
    }
    int nBloking = 64;
    if (argc > 3)
    {
        nBloking = atoi(argv[3]);
    }

    if ((height % nBloking) != 0)
    {
        free(in);
        fprintf(stderr, "error: height = %d, nBloking = %d, (height %% nBloking) = %d\n",
            height, nBloking, (height % nBloking));
        return -1;
    }

    float* out = effect(in, width, height, nQueue, nBloking);
    if (out == NULL)
    {
        free(in);
        fprintf(stderr, "error: effect!\n");
        return -1;
    }

    if (argc < 5)
    {
        fprintf(stdout, "%d %d 1\n", width, height);
        for (int y = 0; y < height; y++)
        {
            for (int x = 0; x < width; x++)
            {
                fprintf(stdout, "%d\n", (int)out[y*width + x]);
            }
        }
    }

    if (in != NULL)
        free(in);

    if (out != NULL)
        free(out);

    return 0;
}
