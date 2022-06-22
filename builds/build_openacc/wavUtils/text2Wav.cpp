//====================================================================
// convert float text to wav faile
//
// g++ -------------------------------------------
//  g++ -o text2Wav text2Wav.cpp Cwav.cpp
//  ./text2Wav input.txt
//
// cl --------------------------------------------
//  cl /EHs /Fe:text2Wav.exe text2Wav.cpp Cwav.cpp
//  text2Wav input.txt
//
// -----------------------------------------------
// program <input.txt> <output.wav> [<m|s>]
//
// cannot be distinguished stereo or monaural from format
//     STEREO:  L, R, L, R, L, R, L, R, L, R, L, R, . . .
//     MONARAL: D, D, D, D, D, D, D, D, D, D, D, D, . . .
//
// (c)Copyright Spacesoft corp., 2018 All rights reserved.
//                                    Kitayama, Hiroyuki
//====================================================================
#include <stdio.h>
#include "Cwav.h"

//--------------------------------------------------------------------
//countLines
size_t
countLines(const char* fname)
{
    FILE  *fp;
    float data;

    if ((fp = fopen(fname, "rt")) == NULL)
        throw "input file open failed.";

    int count = 0;
    while (fscanf(fp, "%f", &data) == 1)
        count++;

    fclose(fp);

    if (count <= 0)
        throw "input file read failed.";

    return count;
}

//--------------------------------------------------------------------
//readData
void
readData(const char* fname, float data[], const size_t length)
{
    FILE *fp;

    if ((fp = fopen(fname, "rt")) == NULL)
        throw "input file open failed.";

    for (size_t i = 0; i < length; i++)
        if (fscanf(fp, "%f", &data[i]) != 1)
            throw "input file read failed.";

    for (size_t i = 0; i < length; i++)
    {
        if (data[i] > 32767.0f || data[i] < -32768.0f)
        {
            fprintf(stderr, "%8d = %10.2f\n", (int)i, data[i]);

            data[i] = min(data[i], 32767.0f);
            data[i] = max(data[i], -32768.0f);
        }
    }
    fclose(fp);
}

//--------------------------------------------------------------------
// main
int
main(int argc, char *argv[])
{
    Cwav cwav;
    float *wav = NULL;
    short *sWav = NULL;
    unsigned int len = 2u;                          // monaural

    try
    {
        if (argc < 3)                               // check parameters
            throw "missing parameters, need <output.wav> and [<m | s>].";

        if (argc == 4)
            if (argv[3][0] == 's' || argv[3][0] == 'S')
            {
                len = 4u;                           // stereo
                fprintf(stdout, "input: stereo.\n");
            }
            else
                fprintf(stdout, "input: monaural.\n");

        size_t wavLength = countLines(argv[1]);
        wav = new float[wavLength];
        sWav = new short[wavLength];

        readData(argv[1], wav, wavLength);          // read text
        for (size_t i = 0; i < wavLength; i++)
            sWav[i] = (short)wav[i];

        cwav.to16bit();
        if (len == 4u)
            cwav.toStereo();
        else
            cwav.toMonaural();
        cwav.setSamplesPerSec(44100);
        cwav.setBytesPerSec(len);
        cwav.setSizeOfData((long)wavLength * 2u);
        cwav.setBitsPerSample(16u);
        cwav.setPWav(sWav);
        cwav.setBlockAlign(len);

        cwav.SaveToFile(argv[2]);               // write wav file

        fprintf(stdout, "convert [%s] to [%s].\n", argv[1], argv[2]);
    }
    catch (char *str)
    {
        fputs(str, stderr);
        fprintf(stderr, "\n");
    }
    if (wav != NULL)
        delete[] wav;

    return 0;
}
