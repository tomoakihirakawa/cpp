namespace BEM
{
  
  using V_d = std::vector<double>;
  using VV_d = std::vector<V_d>;
  using VVV_d = std::vector<VV_d>;

  using netP = networkPoint;  
  using V_netPp = std::vector<networkPoint*>;  
  using VV_netPp = std::vector<V_netPp>;

  using netF = networkFace;
  
  double N(const double& a, const double& b, const int& i){
    switch(i){
    case 0:
      return a;
    case 1:
      return b;
    case 2:
      return 1.-a-b;
    default:
      Print(__PRETTY_FUNCTION__,Red);
      return 0;      
    }       
  };
  
  V_d N(const double& a, const double& b){
    return {a, b, 1.-(a+b)};
  };

  V_d dNd_(const int& i, const double& a, const double& b){
    switch(i){
    case 0:
      return {1, b, -(1+b)};
    case 1:
      return {a, 1, -(a+1)};
    default:
      Print(__PRETTY_FUNCTION__,Red);
      Print(ERROR+" dxds takes 0 or 1 ",Red);
      return {0.,0.,0.};
    }
  };
  
  V_d X(const V_netPp ps, const double& a, const double& b){
    return Dot(N(a,b), obj3D::extractX(ps));
  };

  V_d X(const netF* f, const double& a, const double& b){
    return Dot(N(a,b),f->getLocations());
  };
  
  V_d X(const VV_d& xyzs, const double& a, const double& b){
    return Dot(N(a,b),xyzs);
  };

  V_d dXd_(const V_netPp ps, const double& a, const double& b, const int& i){
    return Dot(dNd_(i,a,b), obj3D::extractX(ps));
  };
  
  V_d dXd_(const netF* f, const double& a, const double& b, const int& i){
    return Dot(dNd_(i,a,b), f->getLocations());
  };

  V_d dXd_(const VV_d& xyzs, const double& a, const double& b, const int& i){
    return Dot(dNd_(i,a,b), xyzs);
  };
  
  VV_d dNd_(const double& a, const double& b){
    return {{1, b, -(1+b)}, {a, 1, -(a+1)}};
  };

  // V_d nabla_phi(const V_d& phi, const double& phi_n, const VV_d& xyzs, const V_d& normal, const double& a, const double& b){
  //   V_d U_Local = Dot(phi, dNd_(a,b));// nabra_phi_local_horizon
  //   U_Local.insert(U_Local.begin(),phi_n);// local velocity

  //   V_d dxd_0 = dXd_(xyzs,a,b,0);
  //   V_d dxd_1 = dXd_(xyzs,a,b,1);

  //   return Dot({dxd_0,dxd_1,normal}, U_Local);
  // };
  
  using map_P_Vd = std::map<netP*, V_double>;
  using map_P_d = std::map<netP*, double>;

  /// THIS ONLY USES ACTIVATED FACES
  map_P_d DphiDt(map_P_Vd& P_phiphin, const bool TorF=true){
    map_P_d ret;
    for(const auto& [p, phiphin]:P_phiphin){
      V_double dphidt;


      
      for(const auto& f:p->getFaces(TorF)){
	V_d Phi(3);
	int i=0;
	double a,b;    
	for(const auto& q:f->getPoints()){
	  if(q==p){if(i==0){a=1.;b=0.;}else if(i==1){a=0.;b=1.;}else if(i==2){a=0.;b=0.;}else{Print(__PRETTY_FUNCTION__);abort();}}
	  Phi[i++] = P_phiphin[q][0];
	}    
	V_d u = Dot(Phi, dNd_(a,b));// nabra_phi_local_horizon
        dphidt.push_back( ( u[0]*u[0] + u[1]*u[1] + phiphin[1]*phiphin[1] ) / 2. - p->getX()[2] );	
      }

      ret[p] = Mean(dphidt);
    }

    return ret;
  };

  /// THIS ONLY USES ACTIVATED FACES
  map_P_Vd nablaPhi(map_P_Vd& P_phiphin, const bool TorF=true){
    map_P_Vd ret;
    for(const auto& [p, phiphin]:P_phiphin){
      VV_double nablaPhi;
      double phi_n = phiphin[1], z = p->getX()[2];
      for(const auto& f:p->getFaces(TorF)){
	V_d Phi(3);
	int i=0;
	double a,b;
	
	for(const auto& q:f->getPoints()){
	  if(q==p){if(i==0){a=1.;b=0.;}else if(i==1){a=0.;b=1.;}else if(i==2){a=0.;b=0.;}else{Print(__PRETTY_FUNCTION__);abort();}}
	  Phi[i++] = P_phiphin[q][0];
	}
	
	VV_d M = {dXd_(f,a,b,0)/Norm(dXd_(f,a,b,0)),
		  dXd_(f,a,b,1)/Norm(dXd_(f,a,b,1)),
		  f->getNormal()/Norm(f->getNormal())};	
	V_d u = Dot(Phi,dNd_(a,b));// nabra_phi_local_horizon
	
        nablaPhi.push_back(Dot(Inverse(M),{u[0],u[1],phi_n}));
      }
      ret[p] = Mean(nablaPhi);
    }

    return ret;
  };
  
  double IGn(const netF* f, const netP* origin, const VV_d& gw = GaussianQuadratureWeights(5,-1.,1.)){
    double ret=0;
    V_d v, r;
    for(const auto& tw0:gw)
      for(const auto& tw1:gw){
	v = tw1*(1-tw0[0]);
	r = X(f, tw0[0],v[0])-origin->getX();
	for(const auto& n:N(tw0[0],v[0]))
	  ret -= Dot(r/pow(Norm(r),3), Cross(dXd_(f, tw0[0],v[0],0),dXd_(f, tw0[0],v[0],1))) * n * tw0[1] * v[1];
      }
    return ret;
  };
  
  double IG(const netF* f, const netP* origin, const VV_d& gw = GaussianQuadratureWeights(5,-1.,1.)){
    double ret=0, detJ;
    V_d v, r;
    for(const auto& tw0:gw)
      for(const auto& tw1:gw){
	v = tw1*(1-tw0[0]);
	r = X(f, tw0[0],v[0]) - origin->getX();
	detJ = Dot({0.,0.,1.},Cross(dXd_(f, tw0[0],v[0],0),dXd_(f, tw0[0],v[0],1)));
	for(const auto& n:N(tw0[0],v[0]))
	  ret += 1./Norm(r) * detJ * n * tw0[1] * v[1];
      }
    return ret;
  };

  V_d IGIGn(const netF* f, const netP* origin, const VV_d& gw = GaussianQuadratureWeights(5,-1.,1.)){
    double ig=0, ign=0, detJ;
    V_d v, r;
    for(const auto& tw0:gw)
      for(const auto& tw1:gw){
	v = tw1*(1-tw0[0]);
	r = X(f, tw0[0],v[0]) - origin->getX();
	detJ = Dot({0.,0.,1.},Cross(dXd_(f, tw0[0],v[0],0),dXd_(f, tw0[0],v[0],1)));
	for(const auto& n:N(tw0[0],v[0])){
	  ig += 1./Norm(r) * detJ * n * tw0[1] * v[1];
	  ign -= Dot(r/pow(Norm(r),3), Cross(dXd_(f, tw0[0],v[0],0),dXd_(f, tw0[0],v[0],1))) * n * tw0[1] * v[1];
	}
      }
    return {ig,ign};
  };
  
  std::map<netP*, V_d> sum(const std::map<netF*, std::map<netP*, V_d>>& IGIGn){
    std::map<netP*, V_d> ret;

    for( auto const& [f, p_igign] : IGIGn){
      for( auto const& [p, igign] : p_igign){      
	ret[p] += igign;
      }
    }      
    return ret;
  };

  VV_netPp triangulate(const V_netPp& pointsOfPoly, const V_d& normal){
    VV_netPp ret;
    for(const auto& index:(geometry::polygon(obj3D::extractX(pointsOfPoly))).triangulate(normal))
      ret.emplace_back(V_netPp{pointsOfPoly[index[0]], pointsOfPoly[index[1]], pointsOfPoly[index[2]]});

    return ret;
  };

  std::map<netP*, V_d> calc_P_IGIGn(const netF* f, const netP* origin, const VV_d& gw = GaussianQuadratureWeights(5,-1.,1.)){
    double ig_=0, ign_=0, ig=0, ign=0, detJ, a, b, norm_r, weight;
    V_d v, r, igign={0,0}, cross={0,0,0}, A = origin->getX();
    V_netPp ps = f->getPoints();
    std::map<netP*, V_d> ret;
    ret[ps[0]] = {0,0};
    ret[ps[1]] = {0,0};
    ret[ps[2]] = {0,0};    

    if(f->intersectQ()){

      VV_d M;
      VV_d invX = Inverse(Transpose(f->getLocations()));
      
      for(const auto& routeP:f->getPointsCutFacesBehind()){
  	for(const auto& tri_ps:triangulate(routeP,f->getNormal())){
	  
	  M = Dot(obj3D::extractX(tri_ps), invX);
	  
	  for(const auto& tw0:gw){
	    a = tw0[0];
	    for(const auto& tw1:gw){	
	      v = tw1*(1-a);
	      weight = tw0[1] * v[1];
	      b = v[0];
	      
	      norm_r = Norm(r = X(f, a, b) - A);
	      cross = Cross(dXd_(tri_ps, a, b, 0),dXd_(tri_ps, a, b, 1));

	      detJ = Dot({0.,0.,1.}, cross);
	      igign = {1./norm_r * detJ * weight, -Dot(r/pow(norm_r,3), cross) * weight};
		
	      for(const auto& shape:N(a,b)){
	 
		ret[ps[0]] += igign * shape * Dot(M[0], N(a, b));
		ret[ps[1]] += igign * shape * Dot(M[1], N(a, b));
		ret[ps[2]] += igign * shape * Dot(M[2], N(a, b));
	  
	      }
	    }	  
	  }
	}
      }
      
    }else{

      for(const auto& tw0:gw){
	a = tw0[0];	  
	for(const auto& tw1:gw){	
	  v = tw1*(1-a);
	  weight = tw0[1] * v[1];
	  b = v[0];

	  norm_r = Norm(r = X(f, a, b) - A);
	  cross = Cross(dXd_(f, a, b,0),dXd_(f, a, b, 1));

	  detJ = Dot({0.,0.,1.}, cross);

	  igign = {1./norm_r * detJ * weight, -Dot(r/pow(norm_r,3), cross) * weight};

	  for(const auto& shape:N(a,b)){
	  
	    ret[ps[0]] += igign * shape * N(a, b, 0);
	    ret[ps[1]] += igign * shape * N(a, b, 1);
	    ret[ps[2]] += igign * shape * N(a, b, 2);
	    
	  }
	}
      }
      
    }

    return ret;
  };  
  
}
