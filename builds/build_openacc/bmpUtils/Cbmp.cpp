//
// Cbmp.cpp
//
// (c)Copyright Spacesoft corp., 2018 All rights reserved.
//                               Hiro KITAYAMA
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <sys/stat.h>           // for SIDBA
#include "bitmapStruct.h"
#include "Cbmp.h"


//-------------------------------------------------------------------
// コンストラクタ
Cbmp::Cbmp()
: mPdib(NULL), mPbitmap(NULL), mDibSize(0), mRowPitch(0),
                        mPixelPitch(0), mImageSize(0), mAbsHeight(0)
{
    assert(sizeof(char) ==1);
    assert(sizeof(short)==2);
    assert(sizeof(int)  ==4);
}


//-------------------------------------------------------------------
// デストラクタ
Cbmp::~Cbmp()
{
    SAFE_FREE(mPdib);                                // free bmp
}


//================ vvvvvv private vvvvvv ============================

//-------------------------------------------------------------------
// read bitmap file header
//
// return true :0
//        false:!0=error #
int
Cbmp::readHeader(FILE* fp)
{
    if(fread(&mBmpFileHdr, sizeof(bmpFileHdr), 1, fp)!=1)
        return -1;

    if(mBmpFileHdr.bfType!='B'+'M'*256)
        return -2;                                  // not bitmap file

    return 0;
}


//-------------------------------------------------------------------
// read bitmap body
int
Cbmp::readDib(FILE* fp)
{
    if(fread(mPdib , mDibSize, 1, fp)!=1)           // read body
        return -1;

    if(mPdib->biBitCount!=16
                    && mPdib->biBitCount!=24
                                && mPdib->biBitCount!=32)
        return -2;                                  // not 16/24/32bpp

    return 0;
}


//-------------------------------------------------------------------
// write bitmap file header
int
Cbmp::writeHeader(FILE* fp)
{
    if(fwrite(&mBmpFileHdr, sizeof(bmpFileHdr), 1, fp)!=1)
        return -1;

    return 0;
}


//-------------------------------------------------------------------
// write bitmap file body
int
Cbmp::writeDib(FILE* fp)
{
    if(fwrite(mPdib , mDibSize, 1, fp)!=1)          // write bitmap body
        return -1;

    return 0;
}


//-------------------------------------------------------------------
// set bitmap file header
void
Cbmp::setBmpInfoHdr(const int width, const int height)
{
    mPdib->biSize         =sizeof(bmpInfoHdr);
    mPdib->biWidth        =width;
    mPdib->biHeight       =height;
    mPdib->biPlanes       =1;
    mPdib->biBitCount     =24;                      // 24 bpp
    mPdib->biCompression  =0;
    mPdib->biSizeImage    =0;
    mPdib->biXPelsPerMeter=0;
    mPdib->biYPelsPerMeter=0;
    mPdib->biClrUsed      =0;
    mPdib->biClrImportant =0;
}


//-------------------------------------------------------------------
// set bitmap info header
//
// set bitmap file header
// set mAbsHeight
// set mPixelPitch
// set mRowPitch
//
void
Cbmp::setBmpFileHdr(const int width, const int height)
{
    mAbsHeight=height>0 ? height : -(height);       //abs

    mPixelPitch=3;                                  // 24 bpp

    mRowPitch=width*mPixelPitch;                    // to 4byte boundary
    if(mRowPitch%4)
        mRowPitch=mRowPitch+(4-(mRowPitch%4));

    mBmpFileHdr.bfType='B'+'M'*256;
    mBmpFileHdr.bfSize=(mRowPitch*mAbsHeight)+sizeof(bmpFileHdr)+sizeof(bmpInfoHdr);
    mBmpFileHdr.bfReserved1=0;
    mBmpFileHdr.bfReserved2=0;
    mBmpFileHdr.bfOffBits=sizeof(bmpFileHdr)+sizeof(bmpInfoHdr);
}


//================ ^^^^^^ private ^^^^^^ ====================================




//-------------------------------------------------------------------
// load bitmap image from file
void
Cbmp::loadFromFile(const char* bmpFName)
{
    FILE* fp;
    struct stat statbuf;                            // for SIDBA

    SAFE_FREE(mPdib);                               // delete image

    if ((fp = fopen(bmpFName, "rb")) == 0)          // open bitmap file
        throw "input file open failed.";

    if (stat(bmpFName, &statbuf) != 0)              // for SIDBA
        throw "function stat() failed.";            // for SIDBA

    if (readHeader(fp) != 0)                        // read file header
    {
        fclose(fp);
        throw "failed to read bitmap file header.";
    }

    //mDibSize=mBmpFileHdr.bfSize-sizeof(bmpFileHdr); // size of dib
    mDibSize = statbuf.st_size - sizeof(bmpFileHdr);    // for SIDBA
    mPdib = (bmpInfoHdr *)malloc(mDibSize);         // alloc dib memory

    if (readDib(fp) != 0)                           // read dib
    {
        SAFE_FREE(mPdib);
        fclose(fp);
        throw "failed to read bitmap file body.";
    }
    fclose(fp);                                     // close bitmap file

    mPbitmap = (unsigned char *)(mPdib)             // move pos. to body
        +mBmpFileHdr.bfOffBits
        - sizeof(bmpFileHdr);

    mPixelPitch = mPdib->biBitCount / 8;

    mRowPitch = (mPdib->biWidth*mPixelPitch);       // clac. row pitch by bytes
    if (mRowPitch % 4 != 0)
        mRowPitch += (4 - (mRowPitch % 4));

    mAbsHeight = mPdib->biHeight > 0 ? mPdib->biHeight : -(mPdib->biHeight);   //abs
    mImageSize = mRowPitch*mAbsHeight;
}


//-------------------------------------------------------------------
// get mem addr of specified scanrow#
unsigned char*
Cbmp::getScanRow(const int rowNo) const
{
    int absrowNo;

    if(mPdib==0)
        return 0;

    absrowNo=rowNo;
    if(mPdib->biHeight<0)
        absrowNo=mPdib->biHeight-rowNo-1;

    return (mPbitmap+(absrowNo*mRowPitch));
}


//-------------------------------------------------------------------
// save to bitmap file
void
Cbmp::saveToFile(const char* bmpFName)
{
    FILE* fp;

    if((fp=fopen(bmpFName, "wb"))!=0)               // open file
    {
        if(writeHeader(fp)==0)                      // write header
        {
            if(writeDib(fp)!=0)                     // write dib
                throw "failed to write dib.";
        }
        else
            throw "failed to write header.";
    }
    else
        throw "failed to open file.";

    fclose(fp);
    SAFE_FREE(mPdib);
}


//-------------------------------------------------------------------
// convert color to gray scale
void
Cbmp::getGSData(unsigned char* gs) const
{
    unsigned char* pRow=mPbitmap;
    unsigned char* pDest=gs;

    for(int y=0 ; y<mAbsHeight ; y++)
    {
        for(int x=0; x<getWidth(); x++)
        {
            float m=  (float)pRow[(x*mPixelPitch)+0]*0.114478f  // blue
                    + (float)pRow[(x*mPixelPitch)+1]*0.586611f  // green
                    + (float)pRow[(x*mPixelPitch)+2]*0.298912f; // red

            *pDest=(unsigned char)m;                            // gray scale
            pDest++;
        }
        pRow+=mRowPitch;
    }
}


//-------------------------------------------------------------------
// get BGRA
void
Cbmp::getBgraData(unsigned char* dataBgra) const
{
    if(mPdib->biBitCount==32)
        memcpy(dataBgra, mPbitmap, mImageSize);         // copy BGRA to dest.
    else
    {
        int index=0;                                    // cnvert 24bpp to BGRA
        for(int i=0; i<mImageSize; i+=3)
        {
            dataBgra[index++]=mPbitmap[i+0];    // B
            dataBgra[index++]=mPbitmap[i+1];    // G
            dataBgra[index++]=mPbitmap[i+2];    // R
            dataBgra[index++]=255;              // A
        }
    }
}


//-------------------------------------------------------------------
// gray scale to 24bpp RGB or 32bpp BGRA
void
Cbmp::gs2bgra(unsigned char* gs) const
{
    for(int y=0; y<mAbsHeight; y++)
    {
        int rowOffset=y*mRowPitch;
        for(int x=0; x<mPdib->biWidth; x++)
        {
            mPbitmap[rowOffset+(mPixelPitch*x)+0]=gs[(y*mPdib->biWidth)+x];  // B
            mPbitmap[rowOffset+(mPixelPitch*x)+1]=gs[(y*mPdib->biWidth)+x];  // G
            mPbitmap[rowOffset+(mPixelPitch*x)+2]=gs[(y*mPdib->biWidth)+x];  // R
                                                                             // A
        }
    }
}


//-------------------------------------------------------------------
// create 24 bit DIB
int
Cbmp::create24Dib(const int width, const int height)
{
    setBmpFileHdr(width, height);

    SAFE_FREE(mPdib);                               // delete bmp
    mDibSize=mBmpFileHdr.bfSize-sizeof(bmpFileHdr); // size of dib
    mPdib=(bmpInfoHdr *)malloc(mDibSize);           // alloc dib memory

    setBmpInfoHdr(width, height);

    mPbitmap=(unsigned char *)(mPdib)               // move pos. to body
                            +mBmpFileHdr.bfOffBits
                                    -sizeof(bmpFileHdr);

    mImageSize=mRowPitch*mAbsHeight;

    memset(mPbitmap, 0xFF, mImageSize);             // init. image data

    return 0;
}


//-------------------------------------------------------------------
// easy format analyzer
int
Cbmp::easyFmtAna(void) const
{
    if(mPdib==NULL)
        return BMTFMT_UNKNOWN;

    if(mPdib->biBitCount<16)
        return BMTFMT_UNKNOWN;

    if(mPdib->biBitCount==16 )
    {
        if(mPdib->biCompression==0)
            return BMTFMT_16_555;
        else
            return BMTFMT_16_565;
    }

    if(mPdib->biBitCount==24 && mPdib->biCompression==0)
        return BMTFMT_24_888;

    if(mPdib->biBitCount==32)
    {
        if(mPdib->biCompression==0)
            return BMTFMT_32_BGRA;
        else
            return BMTFMT_32_BGRX;
    }

    return BMTFMT_UNKNOWN;
}
