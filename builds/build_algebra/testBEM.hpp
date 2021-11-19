using V_double = std::vector<double>;
using VV_double = std::vector<std::vector<double>>;
using VVV_double = std::vector<std::vector<std::vector<double>>>;

VV_double gw;
virtual V_double N(const double& a, const double& b){
  return {a, b, 1.-(a+b)};
};
virtual V_double X(networkFace* f, const double& a, const double& b){
  return Dot({a, b, 1.-(a+b)},f->getLocations());
};
virtual V_double dXds(networkFace* f, const double& a, const double& b, const int& i){
  switch(i){
  case 0:
    return Dot({1, b, -(1+b)}, f->getLocations());
  case 1:
    return Dot({a, 1, -(a+1)}, f->getLocations());
  default:
    Print(ERROR+" dxds takes 0 or 1 ",Red);
    return {0.};
  }
};

virtual double IGn(networkFace* f, networkPoint* origin){
  double ret=0;
  V_double v, r;
  for(const auto& tw0:gw)
    for(const auto& tw1:gw){
      v = tw1*(1-tw0[0]);
      r = X(f, tw0[0],v[0])-origin->xyz;
      for(const auto& n:N(tw0[0],v[0]))
	ret -= Dot(r/pow(Norm(r),3), Cross(dXds(f, tw0[0],v[0],0),dXds(f, tw0[0],v[0],1))) * n * tw0[1] * v[1];
    }
  return ret;
};
virtual double IG(networkFace* f, networkPoint* origin){
  double ret=0, detJ;
  V_double v, r;
  for(const auto& tw0:gw)
    for(const auto& tw1:gw){
      v = tw1*(1-tw0[0]);
      r = X(f, tw0[0],v[0]) - origin->xyz;
      detJ = Dot({0.,0.,1.},Cross(dXds(f, tw0[0],v[0],0),dXds(f, tw0[0],v[0],1)));
      for(const auto& n:N(tw0[0],v[0]))
	ret += 1./Norm(r) * detJ * n * tw0[1] * v[1];
    }
  return ret;
};
void calculateIG(networkPoint* p, networkPoint* origin){
  p->IG = 0;
  for(const auto& f:p->getFaces()){
    p->IG += this->IG(f, origin);
  };    
};
void calculateIGn(networkPoint* p, networkPoint* origin){
  p->IGn = (p==origin) ? -2*M_PI : 0;
  for(const auto& f:p->getFaces()){
    p->IGn += this->IGn(f, origin);
  };    
};




class testBEM{
public:
  
  double length(const VV_double& vv){  
    return Norm(vv[0]-vv[1]);
  };
  double length(const networkLine* l){  
    return length(l->getLocations());
  };
  
  testBEM(){
    Print("**************** testBEM ****************", Red);
    NetworkObj obj("./obj/tank.obj");
    NetworkW water({6,6},{1.133,1.13,1/2.},.2);
    GNUPLOT plot;

    for(auto i=0; i<16; i++){
      plot.Set({{"style","arrow 1 nohead lc \"blue\" lw .5"}});
      plot.Set({{"style","arrow 2 nohead lc \"magenta\" lw .5 dt 2"}});    
      plot.Set({{"style","arrow 3 nohead lc \"red\" lw .5 dt 2"}});
      plot.Set({{"style","arrow 4 nohead lc \"green\" lw 2 dt 2"}});
      plot.Set({{"style","arrow 5 nohead lc \"orange\" lw 3 dt 2"}});
      plot.Set({{"style","arrow 10 nohead lc \"red\" lw .5"}});            
      plot.Set({{"key",""}});
      
      Print(obj.Points.size(),Red);
      Network corssNetwork(obj,water);
      network::setStatus(corssNetwork.Points,true);
      
      plot.plot3d();
      std::cin.ignore();
      plot.Clear();      
    }
  };
  ~testBEM(){Print(CHECK+" testBEM",Red);};
};
