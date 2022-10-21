//
// rotate async
//
// (c)Copyright Spacesoft corp., 2018 All rights reserved.
//                               Hiro KITAYAMA
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#ifndef M_PI
#define M_PI 3.14159265358979f              // pi
#endif //M_PI

#define radian2degree(a) ((a)/M_PI*180.0)   // radian to degree
#define degree2radian(a) ((a)/180.0*M_PI)   // degree to radian

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
// rotate
float*
effect(const float* in, const int width, const int height,
                const float degree, const int nQueue, const int nBloking)
{
    int size = height * width;
    unsigned int datasize = sizeof(float)*size;
    float* out = (float *)malloc(datasize);
    int blockSize = height / nBloking;
    int blockY, outY, outX;
    clock_t start, stop;

    float radian = (float)degree2radian(degree);    // clockwise

    int yc = height / 2;    // y center
    int xc = width / 2;     // x center

    start = clock();

    #pragma acc data copyin(in[:size]) create(out[:size])
    for (blockY = 0; blockY < nBloking; blockY++)
    {
        int startY = -yc + (blockY * blockSize);
        int endY = startY + blockSize;
        int hAsync = (blockY % nQueue) + 1;

        #pragma acc parallel loop collapse(2) gang vector async(hAsync)
        for (outY = startY; outY < endY; outY++)
        for (outX = -xc; outX < width - xc; outX++)// dest x coord
        {
            float inY = (float)(outX*sin(radian) + outY * cos(radian));
            float inX = (float)(outX*cos(radian) - outY * sin(radian));

            int inFixY = inY > 0.0f ? (float)floor(inY + 0.5f) : (float)(-1.0*floor(fabs(inY) + 0.5));
            int inFixX = inX > 0.0f ? (float)floor(inX + 0.5f) : (float)(-1.0*floor(fabs(inX) + 0.5));

            float q = inY - (float)inFixY;
            float p = inX - (float)inFixX;

            inFixX += xc;
            inFixY += yc;
            int oX = outX + xc;
            int oY = outY + yc;

            int dstX = oX;
            int dstY = oY * width;
            int dst = dstY + dstX;

            if (inFixY >= 0 && inFixY < height - 1
                && inFixX >= 0 && inFixX < width - 1)
            {
                int srcX0 = inFixX;
                int srcX1 = srcX0 + 1;
                int srcY0 = inFixY * width;
                int srcY1 = srcY0 + width;

                int src00 = srcY0 + srcX0;
                int src01 = srcY0 + srcX1;
                int src10 = srcY1 + srcX0;
                int src11 = srcY1 + srcX1;

                float data = (1.0f - q)*((1.0f - p)*(float)in[src00]
                                              + p * (float)in[src01])
                                  + q * ((1.0f - p)*(float)in[src10]
                                              + p * (float)in[src11]);

                if (data > 255.0f) data = 255.0f;
                if (data < 0.0f) data = 0.0f;
                out[dst] = data;
            }
            else
            {
                out[dst] = 255.0f;
            }
        }
        #pragma acc update self(out[(startY+yc)*width:blockSize*width]) async(hAsync)
    }
    #pragma acc wait

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
    float degree = 33.3f;
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

    fprintf(stdout, "%d %d 1\n", width, height);
    if (argc > 2)
        degree = (float)atof(argv[2]);

    int nQueue = 4;
    if (argc > 3)
    {
        nQueue = atoi(argv[3]);
    }
    int nBloking = 64;
    if (argc > 4)
    {
        nBloking = atoi(argv[4]);
    }

    if ((height % nBloking) != 0)
    {
        free(in);
        fprintf(stderr, "error: height = %d, nBloking = %d, (height %% nBloking) = %d\n",
            height, nBloking, (height % nBloking));
        return -1;
    }

    float* out = effect(in, width, height, degree, nQueue, nBloking);
    if (out == NULL)
    {
        free(in);
        fprintf(stderr, "error: effect!\n");
        return -1;
    }

    if (argc < 6)
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
