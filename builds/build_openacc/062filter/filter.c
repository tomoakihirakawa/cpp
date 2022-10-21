//
// filter
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
effect(const float* in, const int width, const int height)
{
    int y, x;
    clock_t start, stop;
    float filter[][] =
    {
        { -1.0,  -1.0,  -1.0 },
        { -1.0,   8.0,  -1.0 },
        { -1.0,  -1.0,  -1.0 },
    };
    const int filtersize = sizeof(filter[0])/sizeof(filter[0][0]);

    start = clock();

    int size = width * height;
    unsigned int datasize = sizeof(float)*size;
    float* out = (float *)malloc(datasize);

    #pragma acc parallel loop collapse(2) gang vector \
                copyin(filter, in[:size]) copyout(out[:size])
    for (y = filtersize / 2; y < height - (filtersize / 2); y++)
        for (x = filtersize / 2; x < width - (filtersize / 2); x++)
        {
            float data = 0.0;
            #pragma acc loop independent
            for (int fy = 0; fy < filtersize; fy++)
            {
                long iy = y - (filtersize / 2) + fy;
                #pragma acc loop independent reduction(+:data)
                for (int fx = 0; fx < filtersize; fx++)
                {
                    long ix = x - (filtersize / 2) + fx;

                    data  += filter[fy][fx] * in[iy*width + ix];
                }
            }
            data  = data  <   0.0 ?   0.0: data;
            data  = data  > 255.0 ? 255.0: data;
            out[y * width + x] = data;
        }

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

    float* out = effect(in, width, height);
    if (out == NULL)
    {
        fprintf(stderr, "error: effect!\n");
        return -1;
    }

    if (argc < 3)
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
