//===================================================================
// dump bitmap file
//
// g++ -------------------------------------------
//  g++ -o dumpBmpGray dumpBmpGray.cpp Cbmp.cpp
//  ./dumpBmpGray input.bmp
//
// cl --------------------------------------------
//  cl /EHs /Fe:dumpBmpGray.exe dumpBmpGray.cpp Cbmp.cpp
//  dumpBmpGray input.bmp
//
// -----------------------------------------------
// dumpBmpGray <input.bmp>
//
// output format:
//  width height ch.=1
//    Data
//    Data
//    Data
//    :     :     :
//    Data
//
// (c)Copyright Spacesoft corp., 2018 All rights reserved.
//                               Hiro KITAYAMA
//===================================================================
#include <iostream>
#include "Cbmp.h"

using namespace std;

//-------------------------------------------------------------------
// main
int
main(int argc, char* argv[])
{
    try
    {
        if (argc < 2)
            throw "no <input>";

        Cbmp bmp;
        bmp.loadFromFile(argv[1]);                  // load bitmap

        if (bmp.getBitsPerPixcel() != 24)
            throw "bmp must be 24bits per pixcel.";

        int ch = bmp.getBitsPerPixcel() / 8;

        int bmpWidth = bmp.getWidth();
        int bmpHeight = bmp.getAbsHeight();

        cerr << bmpWidth << " x " << bmpHeight << ", " << "1 ch." << endl;
        cout << bmpWidth << " " << bmpHeight << "  1" << endl;

        unsigned char* gs = (unsigned char*)malloc(bmpWidth*bmpHeight);
        bmp.getGSData(gs);

        unsigned char* pData = gs;
        for (int y = 0; y < bmpHeight; y++)
        {
            for (int x = 0; x < bmpWidth; x++)
            {
                cout << " " << (int)*pData++ << endl;
            }
        }
        SAFE_FREE(gs);
    }
    catch (char const *str)
    {
        cerr << str << endl;
        return -1;
    }
    return 0;
}
