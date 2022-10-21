//
//             _________________________
//            |                         |
//            |       mBmpFileHdr       | Bitmap File Header
//            |                         |
//            |_________________________|
// mPdib----->|                         | }
//            |    Bitmap Info Header   | }
//            |                         | }
//            |_________________________| }
// mPbitmap-->|                         | } -> DIB
//            |                         | }
//            |      image  data        | }
//            |        RGB(BGRA)        | }
//            |                         | }
//            |_________________________| }
//
// (c)Copyright Spacesoft corp., 2018 All rights reserved.
//                               Hiro KITAYAMA
#ifndef __CBMPH__
#define __CBMPH__

#include "bitmapStruct.h"

//-------------------------------------------------------------------
class Cbmp
{
private:
    // ----- Methods -----
    int readHeader(FILE* fp);
    int readDib(FILE* fp);
    int writeHeader(FILE* fp);
    int writeDib(FILE* fp);
    void setBmpInfoHdr(const int width, const int height);
    void setBmpFileHdr(const int width, const int height);


    // ----- Members -------------------------------
    bmpFileHdr  mBmpFileHdr;                            // ヘッダ



public:
    // ----- Constructor/Destructor ----------------
    Cbmp();                                             // コンストラクタ
    virtual ~Cbmp();                                    // デストラクタ

    // ----- Methods -----
    void loadFromFile(const char* bmpFName);
    int getWidth(void) const { return (mPdib==0 ? 0 : mPdib->biWidth); }
    int getHeight(void) const { return (mPdib==0 ? 0 : mPdib->biHeight); }
    int getAbsHeight(void) const { return (mPdib==0 ? 0 : mAbsHeight); }
    pBmpInfoHdr getPdib(void) const { return mPdib;}
    unsigned char* getPbitmapBody(void) const { return (unsigned char*)(mPdib==0 ? 0 : mPbitmap); }
    unsigned char* getScanRow(const int rowNo) const;
    int getBitsPerPixcel(void) const { return mPdib->biBitCount; }
    void saveToFile(const char* bmpFName);
    void getGSData(unsigned char* gs) const;
    void getBgraData(unsigned char* dataBgra) const;
    void gs2bgra(unsigned char* gs) const;
    int create24Dib(const int width, const int height);
    int easyFmtAna(void) const;


    // ----- Members -------------------------------
    pBmpInfoHdr mPdib;                                  // pointer to BITMAP(DIB)
    unsigned char* mPbitmap;                            // pointer to image
    int mDibSize;                                       // size of BITMAP(DIB)
    int mRowPitch;                                      // row per bytes
    int mPixelPitch;                                    // pixel per bytes
    int mImageSize;                                     // size of image
    int mAbsHeight;                                     // absolute height
};
//-------------------------------------------------------------------

#endif  /* __CBMPH__ */
