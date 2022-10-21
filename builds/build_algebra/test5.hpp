class test5{
public:
  using V_double = std::vector<double>;
  using VV_double = std::vector<V_double>;
  using VVV_double = std::vector<VV_double>;
  
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
  
  test5(){
    gw = GaussianQuadratureWeights(5,-1.,1.);    
    Print("**************** test5 **************** ",Red);
    NetworkObj obj("./obj/tank.obj");
    NetworkW water({2,2},{1.12,1.15,1/2.},.2);  
    GNUPLOT plot;
    //=========== plot ===========
    plot.Set({{"xrange","[-.8:.8]"},{"yrange","[-.8:.8]"},{"zrange","[0:1.2]"}});
    plot.Set({{"style","arrow 1 nohead lc \"blue\" lw 1"}});
    plot.Set({{"style","arrow 2 nohead lc 2 lw 1"}});
    plot.Set({{"style","arrow 3 nohead lc \"red\" lw 2"}});
    for(auto i=0; i<50; i++){//ここが増えるとより多くcrossをチェックするため，crossの周辺だけをより細かくする
    
      plot.Set({{"key",""},{"title","\"divide "+std::to_string(i)+"\""}});
      Network corssNetwork(obj,water);    
      int counter=0;

      // {
      //   std::vector<std::vector<std::vector<double>>> vecvec; 
      //   for(const auto& p:water.getPoints())
      // 	vecvec.push_back({p->xyz, p->getNormal()});
      //   plot.SaveVectorData(vecvec,{{"notitle",""}});      
      // }  
      
      {
      	V_double cross;
      	VVV_double tmp;
      	double s = 4;
        for(const auto& f:water.Faces){
      	  for(auto i=0; i<(int)(s+1); i++){
      	    for(auto j=0; j<(int)(s+1)-i; j++){
      	      cross = Cross(dXds(f, i/s, j/s, 0),dXds(f, i/s, j/s, 1));
      	      tmp = {{ X(f, i/s, j/s), cross/Norm(cross)/5. }};
      	      plot.SaveVectorData(tmp, {{"notitle",""},{"lc", plot.rgb(235.*V_double{i/s, j/s, 1.-(i+j)/s})}});
      	    }
      	  }
      	}
      }
      
      // {
      //   std::vector<std::vector<std::vector<double>>> vecvec; 
      //   for(const auto& f:water.Faces[5])
      // 	  vecvec.push_back({f->X(.5,.5), f->dXds(.5,.5,1)});
      //   plot.SaveVectorData(vecvec,{{"notitle",""}});      
      // }

      for(const auto& p:water.getPoints()){
	double x = p->xyz[0];
	double y = p->xyz[1];
	p->xyz += {0,0,cos(M_PI*(x*y)+i/4.)/30.};
      }
      
      plot.SaveVectorData(getVectorData(water.getFaces()),{{"arrowstyle","1"},{"title","water"}});      
      plot.SaveVectorData(getVectorData(obj.getLines()),{{"arrowstyle","2"},{"title","obj"}});          
      plot.SaveVectorData(getVectorData(corssNetwork.getLines()),{{"arrowstyle","3"},{"title","corssNetwork"}});
      plot.SaveData(corssNetwork.getLocations(),{{"w","p"},{"lc","\"magenta\""},{"pt","7"},{"title","cross points"}});    

      network::setStatus(water.getNeighborsIncludeCross(water.getNearestPoints(water.getMeanLocation())),false);
      plot.SaveData(getLocations(water.getPoints(true)),
		    {{"w","p"},{"lc","\"black\""},{"pt","7"},{"title","cross points"}});
    
      plot.Plot3D_All();
      std::cin.ignore();
      plot.Clear();
    }
  };
  ~test5(){Print(CHECK+" test5",Red);};
};
