//===================================================================
// dump bitmap file
//
// g++ -------------------------------------------
//  g++ -o dumpBmp dumpBmp.cpp Cbmp.cpp
//  ./dumpBmp input.bmp
//
// cl --------------------------------------------
//  cl /EHs /Fe:dumpBmp.exe dumpBmp.cpp Cbmp.cpp
//  dumpBmp input.bmp
//
// -----------------------------------------------
// dumpBmp <input.bmp>
//
// output format:
//  width height ch.=3
//    B     G     R
//    B     G     R
//    B     G     R
//    :     :     :
//    B     G     R
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

        cerr << bmpWidth << " x " << bmpHeight << ", " << ch << " ch." << endl;
        cout << bmpWidth << " " << bmpHeight << " " << ch << endl;

        for (int y = 0; y < bmpHeight; y++)
        {
            unsigned char* pData = bmp.getScanRow(y);
            for (int x = 0; x < bmpWidth; x++)
            {
                for (int rgb = 0; rgb < ch; rgb++)
                {
                    cout << " " << (int)pData[x * 3 + rgb];
                }
                cout << endl;
            }
        }
    }
    catch (char const *str)
    {
        cerr << str << endl;
        return -1;
    }
    return 0;
}
