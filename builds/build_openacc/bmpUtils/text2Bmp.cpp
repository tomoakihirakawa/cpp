//===================================================================
// text to bitmap file
//
// プログラム名 <入力ファイル>
//
// g++ -------------------------------------------
//  g++ -o text2Bmp text2Bmp.cpp Cbmp.cpp
//  ./text2Bmp input.txt input.bmp
//
// cl --------------------------------------------
//  cl /EHs /Fe:text2Bmp.exe text2Bmp.cpp Cbmp.cpp
//  dumpBmp input.txt input.bmp
//
// -----------------------------------------------
// text2Bmp <input.txt> <input.bmp>
//
// input format:
//  width height ch.=3 or 1
//    B     G     R            data
//    B     G     R            data
//    B     G     R    or      data
//    :     :     :              :
//    B     G     R            data
//
// (c)Copyright Spacesoft corp., 2018 All rights reserved.
//                               Hiro KITAYAMA
//===================================================================
#include <stdio.h>
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
        FILE *fp;
        int width, height, ch, data, r, g, b;

        if (argc != 3)
            throw "no <input> <output>";

        if ((fp = fopen(argv[1], "rt")) == NULL)
            throw "input file open failed.";

        if (fscanf(fp, "%d %d %d", &width, &height, &ch) != 3)
            throw "input file read failed.";

        Cbmp bmp;
        bmp.create24Dib(width, height);

        for (int y = 0; y < height; y++)
        {
            unsigned char* pData = bmp.getScanRow(y);
            for (int x = 0; x < width; x++)
            {
                if (ch == 1)
                {
                    if (fscanf(fp, "%d", &b) != 1)
                        throw "input file read failed.";
                    r = g = b;
                }
                else
                {
                    if (fscanf(fp, "%d %d %d", &b, &g, &r) != 3)
                        throw "input file read failed.";
                }
                pData[x * 3 + 0] = (unsigned char)b;
                pData[x * 3 + 1] = (unsigned char)g;
                pData[x * 3 + 2] = (unsigned char)r;
            }
        }
        fclose(fp);

        bmp.saveToFile(argv[2]);                    // save bitmap
    }
    catch (char const *str)
    {
        cerr << str << endl;
        return -1;
    }
    return 0;
}
