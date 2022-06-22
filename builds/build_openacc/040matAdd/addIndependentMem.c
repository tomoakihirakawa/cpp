//
//  add two matrixs and store the result in another matrix
//
//  (c)Copyright Spacesoft corp., 2018 All rights reserved.
//                                Hiro KITAYAMA
//
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// main
int
main(int argc, char* argv[])
{
    float *a, *b, *c;
    clock_t start, stop;
    int i, j, n = 4096;

    a = (float*)malloc(sizeof(float *) * n * n);
    b = (float*)malloc(sizeof(float *) * n * n);
    c = (float*)malloc(sizeof(float *) * n * n);

    // initialize array
    for (j = 0; j < n; j++)
    {
        for (i = 0; i < n; i++)
        {
            a[j * n + i] = (float)(i + 1000);
            b[j * n + i] = (float)i / 10.f;
        }
    }


    start = clock();

    // calc.
    for (j = 0; j < n; j++)
    {
        for (i = 0; i < n; i++)
        {
            c[j * n + i] = a[j * n + i] + b[j * n + i];
        }
    }

    stop = clock();

    fprintf(stdout, "      C: ");
    fprintf(stdout, "elapsed time = %.20f [sec]\n",
        (float)(stop - start) / CLOCKS_PER_SEC);


    start = clock();

    // calc.
    #pragma acc data copyin(a[:n*n], b[:n*n]) copyout(c[:n*n])
    {
        #pragma acc kernels loop independent
        for (j = 0; j < n; j++)
        {
            #pragma acc loop independent
            for (i = 0; i < n; i++)
            {
                c[j * n + i] = a[j * n + i] + b[j * n + i];
            }
        }
    }

    stop = clock();

    fprintf(stdout, "OpenACC: ");
    fprintf(stdout, "elapsed time = %.20f [sec]\n",
        (float)(stop - start) / CLOCKS_PER_SEC);

    free(a);
    free(b);
    free(c);

    return 0;
}


/*============================== Windows ==============================*/

/******************************************** pgcc / GTX 650

ši5 4570

PGI$ pgcc -acc -Minfo=all -o addIndependentMem addIndependentMem.c
main:
     55, Generating copyin(a[:n*n])
         Generating copyout(c[:n*n])
         Generating copyin(b[:n*n])
     58, Loop is parallelizable
     61, Loop is parallelizable
         Accelerator kernel generated
         Generating Tesla code
         58, #pragma acc loop gang, vector(4) /* blockIdx.y threadIdx.y * /
         61, #pragma acc loop gang, vector(32) /* blockIdx.x threadIdx.x * /

PGI$ ./addIndependentMem
      C: elapsed time = 0.04100000113248825073 [sec]
OpenACC: elapsed time = 0.13899999856948852539 [sec]


*********************/

/******************************************** pgcc / GTX 750 ti

ši5 6600

PGI$ pgcc -acc -ta=tesla:cc50 -Minfo=accel -o addIndependentMem addIndependentMem.c
main:
     55, Generating copyin(a[:n*n])
         Generating copyout(c[:n*n])
         Generating copyin(b[:n*n])
     58, Loop is parallelizable
     61, Loop is parallelizable
         Accelerator kernel generated
         Generating Tesla code
         58, #pragma acc loop gang, vector(4) /* blockIdx.y threadIdx.y * /
         61, #pragma acc loop gang, vector(32) /* blockIdx.x threadIdx.x * /

PGI$ ./addIndependentMem
      C: elapsed time = 0.02999999932944774628 [sec]
OpenACC: elapsed time = 0.29600000381469726563 [sec]

// OpenACC‚ðŽg‚¤•û‚ª’x‚­‚È‚éA•À—ñ‰»‚Å‚«‚Ä‚¢‚È‚¢B

*********************/




/**************************************************** pgcc / no GPU

ši5 4570

ši5 6600

*********************/



/*============================== Ubuntu ==============================*/

/******************************************** pgcc, pg++

ši5 4570

ši5 6600

*********************/

/******************************************** pgcc, pg++ / GTX 650

ši5 4570

ši5 6600

*********************/



/*========================= Ubuntu on VMware =========================*/

/******************************************** pgcc, pg++

ši5 4570

ši5 6600

*********************/















/*******************
PGI Community Edition 17.10
*******************/
/*

test@ubuntu:/media/test/portableHDD/GenkouPending/OpenACC/accsrc/040matAdd$ g++ -fopenacc -o addIndependent addIndependent.c
test@ubuntu:/media/test/portableHDD/GenkouPending/OpenACC/accsrc/040matAdd$ ./addIndependent
      C: elapsed time = 0.06451799720525741577 [sec]
OpenACC: elapsed time = 0.06378500163555145264 [sec]

*/
//
// --------------------------
//
// PGI$ pgcc -acc -Minfo=accel -o addIndependent.exe addIndependent.c
// main:
//      61, Generating implicit copyin(a[:4096][:4096])
//          Generating implicit copyout(c[:4096][:4096])
//          Generating implicit copyin(b[:4096][:4096])
//      63, Loop is parallelizable
//      66, Loop is parallelizable
//          Accelerator kernel generated
//          Generating Tesla code
//          63, #pragma acc loop gang, vector(4) /* blockIdx.y threadIdx.y */
//          66, #pragma acc loop gang, vector(32) /* blockIdx.x threadIdx.x */
// PGI$ ./addIndependent
//       C: elapsed time = 0.03500000014901161194 [sec]
// OpenACC: elapsed time = 0.51599997282028198242 [sec]
//
//-----
// ­‚µ‘‚­‚È‚Á‚½‚ªAOpenACC‚ª‚Ü‚¾’x‚¢A‰‰ŽZ—Ê‚ª­‚È‚¢‚©H
//
//
