//
// sharpen
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
readImgData(char* fname, int* W, int* H, int *C)
{
    FILE *fp;
    int width, height, ch;
    float r, g, b;

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
    *C = ch;

    if (ch != 3)
    {
        fprintf(stderr, "ch != 3.");
        return NULL;
    }
    unsigned int datasize = sizeof(float)*width*height * ch;
    float * in = (float *)malloc(datasize);

    for (int y = 0; y < height; y++)
    {
        for (int x = 0; x < width; x++)
        {
            if (fscanf(fp, "%f %f %f", &b, &g, &r) != 3)
            {
                fprintf(stderr, "failed read data %s\n", fname);
                return NULL;
            }
            in[(y*width + x)*ch + 0] = b;
            in[(y*width + x)*ch + 1] = g;
            in[(y*width + x)*ch + 2] = r;
        }
    }

    fclose(fp);

    return in;
}

//-------------------------------------------------------------------
// effect
float*
effect(const float* imgData, const int width, const int height, const int ch)
{
    int y, x;
    clock_t start, stop;
    float filter[][] =
    {
        { -1.0,  -1.0,  -1.0 },
        { -1.0,   9.0,  -1.0 },
        { -1.0,  -1.0,  -1.0 },
    };
    const int filtersize = sizeof(filter[0])/sizeof(filter[0][0]);

    start = clock();

    int size = width * height*ch;
    unsigned int datasize = sizeof(float)*size;
    float* out = (float *)malloc(datasize);

    int step = width * ch;
    float blue, green, red;

    #pragma acc parallel loop collapse(2) gang vector \
                copyin(filter, imgData[:size]) copyout(out[:size])
    for (y = filtersize / 2; y < height - (filtersize / 2); y++)
    {
        for (x = filtersize / 2; x < width - (filtersize / 2); x++)
        {
            float blue = 0.0, green = 0.0, red = 0.0;
            #pragma acc loop independent
            for (int fy = 0; fy < filtersize; fy++)
            {
                long iy = y - (filtersize / 2) + fy;
                #pragma acc loop independent \
                    reduction(+:blue) reduction(+:green) reduction(+:red)
                for (int fx = 0; fx < filtersize; fx++)
                {
                    long ix = x - (filtersize / 2) + fx;

                    blue  += filter[fy][fx] * imgData[iy*step + ix * ch + 0];
                    green += filter[fy][fx] * imgData[iy*step + ix * ch + 1];
                    red   += filter[fy][fx] * imgData[iy*step + ix * ch + 2];
                }
            }
            blue  = blue  <   0.0 ?   0.0: blue;
            blue  = blue  > 255.0 ? 255.0: blue;
            green = green <   0.0 ?   0.0: green;
            green = green > 255.0 ? 255.0: green;
            red   = red   <   0.0 ?   0.0: red;
            red   = red   > 255.0 ? 255.0: red;
            out[y * step + x * ch + 0] = blue;
            out[y * step + x * ch + 1] = green;
            out[y * step + x * ch + 2] = red;
        }
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
    int width, height, ch;
    float* in = NULL;

    if (argc < 2)
    {
        fprintf(stderr, "no <input>");
        return -1;
    }

    if ((in = readImgData(argv[1], &width, &height, &ch)) == NULL)
    {
        return -1;
    }

    float* out = effect(in, width, height, ch);
    if (out == NULL)
    {
        fprintf(stderr, "error: effect!\n");
        return -1;
    }

    if (argc < 3)
    {
        fprintf(stdout, "%d %d 3\n", width, height);
        for (int y = 0; y < height; y++)
        {
            for (int x = 0; x < width; x++)
            {
                fprintf(stdout, "%d %d %d\n",
                    (int)out[(y*width + x)*ch + 0],
                        (int)out[(y*width + x)*ch + 1],
                            (int)out[(y*width + x)*ch + 2]);
            }
        }
    }

    if (in != NULL)
        free(in);

    if (out != NULL)
        free(out);

    return 0;
}
