//====================================================================
// dump wav file, stereo/monaural ,int format
//
// Output one data in one line for even stereo or monaural.
//
// g++ -------------------------------------------
//  g++ -o dumpWav dumpWav.cpp Cwav.cpp
//  ./dumpWav input.wav
//
// cl --------------------------------------------
//  cl /EHs /Fe:dumpWav.exe dumpWav.cpp Cwav.cpp
//  dumpWav input.wav
//
// -----------------------------------------------
// program <input.wav>
//
// cannot be distinguished stereo or monaural  from format
//     STEREO:  L, R, L, R, L, R, L, R, L, R, L, R, . . .
//     MONARAL: D, D, D, D, D, D, D, D, D, D, D, D, . . .
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

        unsigned int numOfUnits = cwav.getNumOfUnits();
        short *pMem = (short *)cwav.getPWav();

        for (unsigned int i = 0; i < numOfUnits; i++)   // dump wav
            fprintf(stdout, "%8d\n", (int)pMem[i]);
    }
    catch (char const *str)
    {
        fputs(str, stderr);
        fprintf(stderr, "\n");
    }
    return 0;
}
