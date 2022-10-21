//
//  bitmap structs
//
// (c)Copyright Spacesoft corp., 2018 All rights reserved.
//                               Hiro KITAYAMA
#ifndef __BITMAPSTRUCT__
#define __BITMAPSTRUCT__


#pragma pack(push, 1)

typedef struct
{
    unsigned short  bfType;
    unsigned int    bfSize;
    unsigned short  bfReserved1;
    unsigned short  bfReserved2;
    unsigned int    bfOffBits;
}
bmpFileHdr, *pBmpFileHdr;

typedef struct
{
    unsigned int    biSize;
    int             biWidth;
    int             biHeight;
    unsigned short  biPlanes;
    unsigned short  biBitCount;
    unsigned int    biCompression;
    unsigned int    biSizeImage;
    int             biXPelsPerMeter;
    int             biYPelsPerMeter;
    unsigned int    biClrUsed;
    unsigned int    biClrImportant;
}
bmpInfoHdr, *pBmpInfoHdr;


#pragma pack(pop)


#define BMTFMT_16_555     1
#define BMTFMT_16_565     2
#define BMTFMT_24_888     3
#define BMTFMT_32_BGRA    4
#define BMTFMT_32_BGRX    5
#define BMTFMT_32_101010  6
#define BMTFMT_UNKNOWN    -1

#define SAFE_DELETE(p)    if(p!=NULL) { delete(p); p=NULL; }
#define SAFE_FREE(p)      if(p!=NULL) { free(p);   p=NULL; }

#endif // __BITMAPSTRUCT__
