#include "fundamental.hpp"

#define check1

#if defined(check0)

int main(){

  std::vector<std::vector<double>> ans;
  
  Print("3d cases",Red);
  for(auto p=-4; p<=4; p++){
    std::vector<double> ans_;  
    for(auto q=-4; q<=4; q++){
      std::vector<double> v0{1.,1.,1.};
      std::vector<double> v1{1.,2.,3.};
      std::vector<double> v2{2.,cos(M_PI*q/4.),sin(M_PI*p/4.)};  
      ans_.push_back(VectorAngle(v1,v2,v0));
    }
    ans.emplace_back(ans_);
  }
  Print(ans,Green);  
  
  // mathematicaと比較
  /*
    ClearAll["Global`*"];
    Table[Module[{},
    v0={1.,1.,1.};
    v1={1.,2.,3.};
    v2={2.,Cos[p/4*Pi],Sin[q/4*Pi]};
    N@VectorAngle[v1-v0,v2-v0]
    ],{q,-4,4},{p,-4,4}]
  */

  ans.clear();
  
  Print("3d cases",Red);
  for(auto p=-4; p<=4; p++){
    std::vector<double> ans_;  
    for(auto q=-4; q<=4; q++){
      std::vector<double> v0{1.,1.,1.};
      std::vector<double> v1{1.,2.,3.};
      std::vector<double> v2{2.,cos(M_PI*q/4.),sin(M_PI*p/4.)};  
      ans_.push_back(MathematicaVectorAngle(v1-v0,v2-v0));
    }
    ans.emplace_back(ans_);
  }

  Print(ans,Green);  

  

};

#elif defined(check1)

int main(){

  
  
  std::cout << "-45deg = " << MyVectorAngle({0.,1.,0.}, {1.,0.,0.}, {0.,0.,1.}) << std::endl;

  std::cout << "45deg = " << MyVectorAngle({1.,0.,0.}, {0.,1.,0.}, {0.,0.,1.}) << std::endl;
};


#endif
