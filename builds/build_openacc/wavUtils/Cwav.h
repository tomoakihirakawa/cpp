//==========================================================================
// Class Cwav Header
//
// (c)Copyright Spacesoft corp., 2018 All rights reserved.
//                                    Kitayama, Hiroyuki
//==========================================================================

#ifndef CwavH
#define CwavH

#include "common.h"

static const char *STR_RIFF = "RIFF";
static const char *STR_WAVE = "WAVE";
static const char *STR_fmt  = "fmt ";
static const char *STR_data = "data";

static const int WAV_MONAURAL = 1;
static const int WAV_STEREO   = 2;

//--------------------------------------------------------------------
// structures
#pragma pack(push,1)

typedef struct tagSWaveFileHeader
{
    char                hdrRiff[4];             // 'RIFF'
    unsigned int        sizeOfFile;             // file size - 8
    char                hdrWave[4];             // 'WAVE'
} SWaveFileHeader;

typedef struct tagChank
{
    char                hdr[4];                 // 'fmt ' or 'data'
    unsigned int        size;                   // sizeof(PCMWAVEFORMAT)
                                                //    or Wave data size
} tChank;

typedef struct tagWaveFormatPcm
{
    unsigned short      formatTag;              // WAVE_FORMAT_PCM
    unsigned short      channels;               // number of channels
    unsigned int        samplesPerSec;          // sampling rate
    unsigned int        bytesPerSec;            // samplesPerSec * channels
                                                //        * (bitsPerSample/8)
    unsigned short      blockAlign;             // block align
    unsigned short      bitsPerSample;          // bits per sampling
} tWaveFormatPcm;

typedef struct tagWrSWaveFileHeader
{
    SWaveFileHeader     wfh;                    // Wave File Header
    tChank              cFmt;                   // 'fmt '
    tWaveFormatPcm      wfp;                    // Wave Format Pcm
    tChank              cData;                  // 'data'
} WrSWaveFileHeader;

#pragma pack(pop)


//--------------------------------------------------------------------
// class header
class Cwav
{

private:
    // ----- private member ------------------------------------------
    SWaveFileHeader wFH;
    tWaveFormatPcm  wFP;
    void*           pMem;                       // pointer to wav data
    long            sizeOfData;

    char wavInFName[_MAX_PATH];                 // input  wav file name
    char wavOutFName[_MAX_PATH];                // output wav file name

    // ----- private method ------------------------------------------
    bool readfmtChunk(FILE *fp, tWaveFormatPcm* waveFmtPcm);
    int  wavHeaderWrite(FILE *fp);
    bool wavDataWrite(FILE *fp);


public:
    // ----- Constructor/Destructor ----------------------------------
    Cwav(void);                                 // Constructor
    virtual ~Cwav(void);                        // Destructor

    // ----- public method -------------------------------------------
    void LoadFromFile(const char *wavefile);    // read  wav file
    void SaveToFile(const char *wavefile);      // write wav file
    bool printWavInfo(void);                    // print wav info.

    //----------------------------------------------------------------
    bool isPCM(void)                            // is PCM
    { return wFP.formatTag==1 ? true: false;    }

    //----------------------------------------------------------------
    bool is16bit(void)                          // is 16bits/sample
    { return wFP.bitsPerSample==16 ? true: false;   }

    //----------------------------------------------------------------
    void to16bit(void)                          // to 16bits/sample
    { wFP.bitsPerSample=16;                     }

    //----------------------------------------------------------------
    bool isStereo(void)                         // is stereo
    { return wFP.channels==WAV_STEREO ? true: false;}

    //----------------------------------------------------------------
    void toStereo(void)                         // to stereo
    { wFP.channels=WAV_STEREO;                  }

    //----------------------------------------------------------------
    bool isMonaural(void)                       // is monaural
    { return wFP.channels==WAV_MONAURAL ? true: false;  }

    //----------------------------------------------------------------
    void toMonaural(void)                       // to monaura
    { wFP.channels=WAV_MONAURAL;                }

    //----------------------------------------------------------------
    unsigned int getSamplesPerSec(void)         // get sampling rate
    { return wFP.samplesPerSec;                 }

    //----------------------------------------------------------------
    void setSamplesPerSec(unsigned int samplesPerSec)   // set sampling rate
    { wFP.samplesPerSec=samplesPerSec;          }

    //----------------------------------------------------------------
    void setBytesPerSec(unsigned int bytesPerSec)// set bytes/second
    { wFP.bytesPerSec=bytesPerSec;              }

    //----------------------------------------------------------------
    long getSizeOfData(void)                    // get wav data size
    { return sizeOfData;                        }

    //----------------------------------------------------------------
    void setSizeOfData(long size)               // set wav data size
    { sizeOfData=size;                          }

    //----------------------------------------------------------------
    unsigned short getBitsPerSample(void)       // get bits/sample
    { return wFP.bitsPerSample;                 }

    //----------------------------------------------------------------
    void setBitsPerSample(unsigned short bitsPerSample) // set bits/sample
    { wFP.bitsPerSample=bitsPerSample;          }

    //----------------------------------------------------------------
    void* getPWav(void)                         // get addr. of wav data
    { return pMem;                              }

    //----------------------------------------------------------------
    void setPWav(void* pInMem)                  // set addr. of wav data
    { pMem=pInMem;                              }

    //----------------------------------------------------------------
    unsigned short getBlockAlign(void)          // get blockAlign
    { return wFP.blockAlign;                    }

    //----------------------------------------------------------------
    void setBlockAlign(unsigned short blockAlign)   // set blockAlign
    { wFP.blockAlign=blockAlign;                }

    //----------------------------------------------------------------
    unsigned int getNumOfUnits(void)            // get units of wav data
    { return sizeOfData/(getBitsPerSample()/8); }

    //----------------------------------------------------------------
    unsigned int getNumOfSamples(void)          // get num. of samples
    { return sizeOfData/getBlockAlign();        }

    bool stereo2monaural(void);                 // Stereo -> Monaural
    bool monaural2stereo(void);                 // Monaural -> Stereo
};

//--------------------------------------------------------------------
#endif
