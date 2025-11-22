//====================================================================
// Class Cwav
//
// (c)Copyright Spacesoft corp., 2018 All rights reserved.
//                                    Kitayama, Hiroyuki
//====================================================================
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "Cwav.h"

//--------------------------------------------------------------------
// constractor
Cwav::Cwav(void) : pMem(NULL), sizeOfData(0)
{
    memset(&wFH, 0, sizeof(wFH));           // initialization
    memset(&wFP, 0, sizeof(wFP));
    wavInFName[0] = '\0';
    wavOutFName[0] = '\0';
}


//--------------------------------------------------------------------
// destructor
Cwav:: ~Cwav(void)
{
    SP_FREE(pMem);                          // free WAV data memory
}


/****** vvv private method vvv ******/

//--------------------------------------------------------------------
// read and check fmt chank
bool Cwav::readfmtChunk(FILE *fp, tWaveFormatPcm* waveFmtPcm)
{
    if (fread(waveFmtPcm, sizeof(tWaveFormatPcm), 1, fp) != 1)
        return false;

    return true;
}


//--------------------------------------------------------------------
// write wav header
int Cwav::wavHeaderWrite(FILE *fp)
{
    unsigned short bytes;
    WrSWaveFileHeader wrWavHdr;
    int rCode = -1;

    //RIFF header
    strncpy(wrWavHdr.wfh.hdrRiff, STR_RIFF, sizeof wrWavHdr.wfh.hdrRiff);

    //file size
    wrWavHdr.wfh.sizeOfFile = sizeOfData + sizeof(wrWavHdr) - 8;

    //WAVE header
    strncpy(wrWavHdr.wfh.hdrWave, STR_WAVE, sizeof wrWavHdr.wfh.hdrWave);

    //fmt chunk
    strncpy(wrWavHdr.cFmt.hdr, STR_fmt, sizeof(wrWavHdr.cFmt.hdr));

    //fmt chunk
    wrWavHdr.cFmt.size = sizeof(wrWavHdr.wfp);

    //no compression PCM = 1
    wrWavHdr.wfp.formatTag = 1;

    //ch (mono=1, stereo=2)
    wrWavHdr.wfp.channels = wFP.channels;

    //sampleng rate(Hz)
    wrWavHdr.wfp.samplesPerSec = wFP.samplesPerSec;

    //bytes/sec
    bytes = wFP.bitsPerSample / 8;

    wrWavHdr.wfp.bytesPerSec = bytes*wFP.channels*wFP.samplesPerSec;

    //byte/sample*channels
    wrWavHdr.wfp.blockAlign = bytes*wFP.channels;

    //bit/samples
    wrWavHdr.wfp.bitsPerSample = wFP.bitsPerSample;

    //data chunk
    strncpy(wrWavHdr.cData.hdr, STR_data, sizeof(wrWavHdr.cData.hdr));

    //data length(byte)
    wrWavHdr.cData.size = sizeOfData;

    //write header
    if (fwrite(&wrWavHdr, sizeof(wrWavHdr), 1, fp) == 1)
        rCode = ftell(fp);
    else
        rCode = -1;

    return rCode;
}


//--------------------------------------------------------------------
// write wav content to file
bool Cwav::wavDataWrite(FILE *fp)
{
    if (fwrite(pMem, sizeOfData, 1, fp) != 1)
        return false;

    return true;
}

/****** ^^^ private method ^^^ ******/


//--------------------------------------------------------------------
// read wav file
void Cwav::LoadFromFile(const char* wavefile)
{
    tChank chank;
    long   cursor, len;
    FILE   *fp = NULL;

    try
    {
        wavInFName[0] = '\0';                               // wav file name

        if ((fp = fopen(wavefile, "rb")) == NULL)
            throw "input file open failed.";

        if (fread(&wFH, sizeof(wFH), 1, fp) != 1)           // file header
            throw "error in wav header.";

        if (memcmp(wFH.hdrWave, STR_WAVE, 4) != 0)          // wav header
            throw "error in wav header.";

        if (memcmp(wFH.hdrRiff, STR_RIFF, 4) != 0)
            throw "error in wav header.";

        // 4 byte, bytes after this = (file size - 8)(Byte)
        len = wFH.sizeOfFile;

        while (fread(&chank, sizeof chank, 1, fp) == 1)     // chunk
        {
            if (memcmp(chank.hdr, STR_fmt, sizeof chank.hdr) == 0)
            {
                len = chank.size;
                cursor = ftell(fp);
                if (!readfmtChunk(fp, &wFP))
                    throw "error in wav file format.";
                fseek(fp, cursor + len, SEEK_SET);
            }
            else if (memcmp(chank.hdr, STR_data, 4) == 0)
            {
                sizeOfData = chank.size;
                if ((pMem = malloc(sizeOfData)) == NULL)
                    throw "failed malloc.";

                if (fread(pMem, sizeOfData, 1, fp) != 1)    // read whole
                    throw "failed wav file read.";
            }
            else
            {
                len = chank.size;
                cursor = ftell(fp);
                fseek(fp, cursor + len, SEEK_SET);
            }
        }
        fclose(fp);

        if (!isPCM())                                       // not PCM
            throw "not PCM format.";

        strcpy(wavInFName, wavefile);                       // input wav file name
    }
    catch (char *str)
    {
        SP_FREE(pMem);
        if (fp != NULL)
            fclose(fp);

        throw str;
    }
}

//--------------------------------------------------------------------
// write wav file
void Cwav::SaveToFile(const char *outFile)
{
    FILE *fp = NULL;
    int rCode = 0;

    try
    {
        if ((fp = fopen(outFile, "wb")) == NULL)
            throw "output file open failed.";

        // write wav header
        if (wavHeaderWrite(fp) != sizeof(WrSWaveFileHeader))
            throw "wav header write failed.";

        if (!wavDataWrite(fp))                       // write wav data
            throw "wav data write failed.";

        fclose(fp);

        strcpy(wavOutFName, outFile);               // output wav file name
    }
    catch (char *str)
    {
        SP_FREE(pMem);
        if (fp != NULL)
            fclose(fp);

        throw str;
    }
}

//--------------------------------------------------------------------
// print WAV info
bool Cwav::printWavInfo(void)
{
    printf("         data format: %u (1 = PCM)\n", wFP.formatTag);
    printf("  number of channels: %u\n", wFP.channels);
    printf(" frequency of sample: %u [Hz]\n", wFP.samplesPerSec);
    printf("   bytes per seconds: %u [bytes/sec]\n", wFP.bytesPerSec);
    printf("    bytes x channels: %u [bytes]\n", wFP.blockAlign);
    printf("    bits per samples: %u [bits/sample]\n", wFP.bitsPerSample);
    printf("      wav data size = %lu\n\n", sizeOfData);
    printf("       playback time=%.3f\n", (float)sizeOfData / (float)wFP.bytesPerSec);

    return true;
}

//--------------------------------------------------------------------
// stereo to monaural
bool Cwav::stereo2monaural(void)
{
    setSizeOfData(getSizeOfData() >> 1);
    setBlockAlign(getBlockAlign() >> 1);
    toMonaural();

    return true;
}

//--------------------------------------------------------------------
// monaural to stereo
bool Cwav::monaural2stereo(void)
{
    setSizeOfData(getSizeOfData() << 1);
    setBlockAlign(getBlockAlign() << 1);
    toStereo();

    return true;
}
