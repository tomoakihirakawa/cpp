//
// fir
//
// program <data file> <k file>
// >fir wav.txt k.txt > out.txt
//
// (c)Copyright Spacesoft corp., 2018 All rights reserved.
//                                    Kitayama, Hiroyuki
//
#include <stdio.h>
#include <stdlib.h>
#include <time.h>


//--------------------------------------------------------------------
//countLines
size_t
countLines(const char* fname)
{
    FILE  *fp;
    float data;

    if ((fp = fopen(fname, "rt")) == NULL)
        return 0;

    int count = 0;
    while (fscanf(fp, "%f", &data) == 1)
        count++;

    fclose(fp);

    if (count <= 0)
        return 0;

    return count;
}

//--------------------------------------------------------------------
//readData
void
readData(const char* fname, float * buf, const size_t length)
{
    FILE *fp;

    if ((fp = fopen(fname, "rt")) == NULL)
    {
        fprintf(stderr, "open faild: %s!", fname);
        return;
    }

    for (int i = 0; i < length; i++)
    {
        if (fscanf(fp, "%f", &buf[i]) != 1)
        {
            fprintf(stderr, "read faild: %s!", fname);
            break;
        }
    }
    fclose(fp);
}

//--------------------------------------------------------------------
// main
int
main(int argc, char *argv[])
{
    float *d = NULL, *k = NULL, *z = NULL;
    int dLength, kLength;
    clock_t start, stop;

    if (argc != 3)
    {
        fprintf(stderr, "missing: <data> or <parameter>.\n");
        return -1;
    }

    dLength = (int)countLines(argv[1]);
    kLength = (int)countLines(argv[2]);
    if (dLength <= 0 || kLength <= 0)
    {
        fprintf(stderr, "read faild:%s or %s.\n", argv[1], argv[2]);
        return -1;
    }

    d = (float *)malloc(sizeof(float) * (dLength + kLength - 1));
    k = (float *)malloc(sizeof(float) * kLength);
    z = (float *)malloc(sizeof(float) * dLength);

    readData(argv[1], d, dLength);          // read data
    readData(argv[2], k, kLength);          // read param

    for (int i = dLength; i < dLength + kLength; i++)
        d[i] = (float)0.0;


    start = clock();

    // do fir
    #pragma acc data copyin(d[:dLength+kLength-1], k[:kLength]) copyout(z[:dLength])
    {
        #pragma acc kernels
        #pragma acc loop independent
        for (int n = 0; n < dLength; n++)
        {
            float zz = (float)0.0;
            #pragma acc loop independent
            for (int m = 0; m < kLength; m++)
            {
                zz += (k[m] * d[n + m]);
            }
            z[n] = zz;
        }
    }
    stop = clock();

    fprintf(stderr, "elapsed time = %.20f [sec]\n",
        (float)(stop - start) / CLOCKS_PER_SEC);


    // print result
    for (size_t n = 0; n < (size_t)dLength; n++)
    {
        fprintf(stdout, "%12.4f\n", z[n]);
    }

    free(d);
    free(k);
    free(z);

    return 0;
}
