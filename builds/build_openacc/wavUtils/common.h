//====================================================================
// common.h
//
// (c)Copyright Spacesoft corp., 2018 All rights reserved.
//                                    Kitayama, Hiroyuki
//====================================================================

#ifndef COMMONH__
#define COMMONH__

//--------------------------------------------------------------------
//  macros
#define SP_FREE(p)          if(p) {free(p);     p=NULL;}

#ifndef min
#define min(a,b) (((a)<(b))?(a):(b))
#endif
#ifndef max
#define max(a,b) (((a)>(b))?(a):(b))
#endif

#ifndef _MAX_PATH
#define _MAX_PATH   1024
#endif

//--------------------------------------------------------------------
#endif  /* COMMONH__ */
