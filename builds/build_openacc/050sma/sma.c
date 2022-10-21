//
// sma:Simple Moving Average
//
// program <data file> <n:length>
// >sma wav.txt 16 > out.txt
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
    float *d = NULL, *z = NULL;
    int dLength, smaLength;
    clock_t start, stop;

    if (argc != 3)
    {
        fprintf(stderr, "missing: <data> or <n>.\n");
        return -1;
    }

    smaLength = atoi(argv[2]);
    if (smaLength <= 0)
    {
        fprintf(stderr, "[n] must be grater than 0.\n");
        return -1;
    }
    dLength = (int)countLines(argv[1]);
    if (dLength <= 0)
    {
        fprintf(stderr, "read faild:%s.\n", argv[2]);
        return -1;
    }

    d = (float *)malloc(sizeof(float) * (dLength + smaLength - 1));
    z = (float *)malloc(sizeof(float) * dLength);

    readData(argv[1], d, dLength);          // read data

    for (int i = dLength; i < dLength + smaLength; i++)
        d[i] = (float)0.0;


    start = clock();

    // do sma : simple moving sverage
    #pragma acc data copyin(d[:dLength+smaLength-1]) copyout(z[:dLength])
    {
        #pragma acc kernels
        #pragma acc loop independent
        for (int n = 0; n < dLength; n++)
        {
            float zz = (float)0.0;
            #pragma acc loop independent reduction(+:zz)
            for (int m = 0; m < smaLength; m++)
            {
                zz += d[n + m];
            }
            z[n] = zz / (float)smaLength;
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
    free(z);

    return 0;
}
