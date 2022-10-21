#include "fundamental_vectors.hpp"

#include <iostream>
#include <string>
#include <sstream>
#include <cmath>

#include <fstream>
#include <functional>

using V_d = std::vector<double>;
using VV_d = std::vector<std::vector<double>>;
using VV_i = std::vector<std::vector<int>>;

class WavefrontObj{
public:
  VV_d v;
  VV_i f_v;
  V_d scale;
  std::function<double(V_d)> zfunc;
  
  WavefrontObj(const V_d& scale_IN={1.,1.,1.},
               const std::function<double(V_d)>& zfunc_IN=[](const V_d& a){return 0.;})
    :scale(scale_IN),zfunc(zfunc_IN){};
  
  void clear(){
    v.clear();
    f_v.clear();    
  };
  
  void setType0(const int Nx, const int Ny){
    this->clear();
    double x(0), y(0), z(0);
    int c=0, k=0;
    for(auto i=-Nx; i<=Nx; i++){
      x=double(i);
      for(auto j=-Ny; j<=Ny; j++){
        y=double(j);        
        v.push_back({x, y, z});
        v.push_back({-.5 +x, -.5 +y, z});
        v.push_back({ .5 +x, -.5 +y, z}); 
        f_v.push_back({c=k++,k++,k++});
        //
        v.push_back({.5 +x, -.5 +y, z});
        v.push_back({.5 +x,  .5 +y, z});            
        f_v.push_back({c,k++,k++});
        //
        v.push_back({ .5 +x, .5 +y, z});
        v.push_back({-.5 +x, .5 +y, z});      
        f_v.push_back({c,k++,k++});
        //
        v.push_back({-.5 +x,  .5 +y, z});
        v.push_back({-.5 +x, -.5 +y, z});       
        f_v.push_back({c,k++,k++});
      }
    }

    for(auto& w:v){
      w[0] *= scale[0]/(2.*Nx+1);
      w[1] *= scale[1]/(2.*Ny+1);
      w[2] = scale[2]*this->zfunc({w[0],w[1]});
    };
    
  };

  void setType1(const int Nx, const int Ny){
    this->clear();
    double x, y, z=0;    
    int c=0, k=0;
    for(auto i=-Nx; i<=Nx; i++){
      x=double(i);
      for(auto j=-Ny; j<=Ny; j++){
        y=double(j);
        v.push_back({ .5+x,  .5+y, z});
        v.push_back({-.5+x, -.5+y, z});
        v.push_back({ .5+x, -.5+y, z}); 
        f_v.push_back({c=k++,k++,k++});
        //
        v.push_back({-.5 +x,  .5 +y, z});
        v.push_back({-.5 +x, -.5 +y, z});
        f_v.push_back({c,k++,k++});
      }
    }
    for(auto& w:v){
      w[0] = w[0]*scale[0]/(2.*Nx+1);
      w[1] = w[1]*scale[1]/(2.*Ny+1);
      w[2] = scale[2]*this->zfunc({w[0],w[1]});
    }
  };

  void setType2(const int Nx, const int Ny){
    this->clear();
    double X, Y, x, y, z=0;
    int c=0, k=0;
    for(auto i=-Nx; i<=Nx; i++){
      x=double(i);
      for(auto j=-Ny; j<=Ny; j++){
        y=double(j);

        X = .5 + x;
        Y = .5 + y;
        v.push_back({(X-Y)/sqrt(2),(X+Y)/sqrt(6),z});

        X = -.5 + x;
        Y = -.5 + y;
        v.push_back({(X-Y)/sqrt(2),(X+Y)/sqrt(6),z});

        X = .5 + x;
        Y = -.5 + y;
        v.push_back({(X-Y)/sqrt(2),(X+Y)/sqrt(6),z});
          
        f_v.push_back({c=k++,k++,k++});
        //

        X = -.5 + x;
        Y = .5 + y;
        v.push_back({(X-Y)/sqrt(2),(X+Y)/sqrt(6),z});

        X = -.5 + x;
        Y = -.5 + y;
        v.push_back({(X-Y)/sqrt(2),(X+Y)/sqrt(6),z});

        f_v.push_back({c,k++,k++});
      }
    }
    for(auto& w:v){
      w[0] = w[0]*scale[0]/(2.*Nx+1);
      w[1] = w[1]*scale[1]/(2.*Ny+1);
      w[2] = scale[2]*this->zfunc({w[0],w[1]});
    }
  };

  std::string show(){
    std::stringstream strs;

    strs << "#vertices " << v.size() << std::endl;
    strs << "#faces " << f_v.size() << std::endl;
    
    strs << std::endl;
    
    for(const auto& s:this->v){
      strs << "v" << " ";     
      for(auto it=s.begin(); it<s.end()-1; it++)      
        strs << *it << " ";
      strs << *s.rbegin() << std::endl;    
    }

    strs << std::endl;    
    
    for(const auto& s:this->f_v){
      strs << "f" << " ";
      for(auto it=s.begin(); it<s.end()-1; it++)      
        strs << *it + 1 << " ";
      strs << *s.rbegin() + 1 << std::endl;    
    }

    auto ret = strs.str();
    std::cout << ret << std::endl;
    return ret;
    
  };
  
};


int main(){

  
  WavefrontObj obj({50.,50.,1.},[](const V_d& xy){return sin(std::abs(xy[0])+std::abs(xy[1]));});
  obj.setType2(10,10);
  obj.show();

  std::ofstream ofs("./test.obj");
  ofs << obj.show();
  
}
