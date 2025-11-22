#include "basic_arithmetic_array_operations.hpp"
#include "kernelFunctions.hpp"

extern "C" {
double exported_w_Bspline3(double r, double h) {
   double result = w_Bspline3(r, h);
   return result;
}
}
