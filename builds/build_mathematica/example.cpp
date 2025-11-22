
#include "/Applications/Mathematica.app/Contents/SystemFiles/IncludeFiles/C/WolframLibrary.h"

EXTERN_C DLLEXPORT int addTwoIntegers(WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res) {
   mint i1 = MArgument_getInteger(Args[0]);
   mint i2 = MArgument_getInteger(Args[1]);
   MArgument_setInteger(Res, i1 + i2);
   return LIBRARY_NO_ERROR;
}
