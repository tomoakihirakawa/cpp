#include "fundamental.hpp"

using V_double = std::vector<double>;
using VV_double = std::vector<std::vector<double>>;
using VVV_double = std::vector<std::vector<std::vector<double>>>;

int main(){
  
  VV_double mat = {{1.,2.,3.},
		   {4.,1.,6.},
		   {7.,8.,1.}};

  VV_double mat2 = {{1.,2.,4.},
		    {2.,1.,6.},
		    {4.,4.,1.}};
  
  V_double vec = {1.,2.,3.};

  Print(Inverse(mat),red);
  Print(Dot(mat,vec),blue);
  Print(Dot(vec,mat),green);  
  Print(Dot(mat,mat2),red);
  Print(Dot(mat2,mat),blue);
  
  return 0;      
};
