//====================================================================
// dump wav file, stereo -> monaural, int format
//
// converts stereo to monaural, and outputs one data per line.
//
// g++ -------------------------------------------
//  g++ -o dumpWav2M dumpWav2M.cpp Cwav.cpp
//  ./dumpWav2M input.wav
//
// cl --------------------------------------------
//  cl /EHs /Fe:dumpWav2M.exe dumpWav2M.cpp Cwav.cpp
//  dumpWav2M input.wav
//
// -----------------------------------------------
// program <input.wav>
//
// input: stereo wav file.
//
// output:
//     D, D, D, D, D, D, D, D, D, D, D, D, . . .
//
// (c)Copyright Spacesoft corp., 2018 All rights reserved.
//                                    Kitayama, Hiroyuki
//====================================================================
#include <stdio.h>
#include "Cwav.h"

//--------------------------------------------------------------------
// main
int
main(int argc, char *argv[])
{
    try
    {
        if (argc != 2)                                  // check parameters
            throw "missing input file name.";

        Cwav cwav;

        cwav.LoadFromFile(argv[1]);                     // read WAV file

        if(cwav.isMonaural())
            throw "input av not stereo.";

        unsigned int numOfUnits = cwav.getNumOfUnits();
        short *pMem = (short *)cwav.getPWav();

        for (unsigned int i = 0; i < numOfUnits; i+=2)  // dump wav
        {
            int l = (int)pMem[i];
            int r = (int)pMem[i+1];
            fprintf(stdout, "%8d\n", (l + r) / 2);
        }
    }
    catch (char const *str)
    {
        fputs(str, stderr);
        fprintf(stderr, "\n");
    }
    return 0;
}
