namespace BEM
{
  
  using V_d = std::vector<double>;
  using VV_d = std::vector<V_d>;
  using VVV_d = std::vector<VV_d>;

  using netFp = networkFace*;
  
  using netP = networkPoint;
  using netPp = networkPoint*;    
  using V_netPp = std::vector<networkPoint*>;  
  using VV_netPp = std::vector<V_netPp>;

  using netF = networkFace;

  /// USE LIKE Dot( N(ab), sample )
  V_d N(const V_d& ab){return {ab[0], ab[1], 1.-(ab[0]+ab[1])};};

  V_d dNd_(const int& i){
    switch(i){
    case 0:
      return {1., 0., -1.};
    case 1:
      return {0., 1., -1.};
    default:
      Print(__PRETTY_FUNCTION__);
      abort();
      return {0.,0.,0.};
    }
  };

  VV_d dNd_(){return {{1., 0., -1.},{0., 1., -1.}};};
  
  using map_P_Vd = std::map<netPp, V_d>;
  using map_P_d = std::map<netPp, double>;

  V_d parameterize(const V_netPp& points, const netPp& p){
    int i=0;
    for(const auto& q:points){
      if(q==p){
	if(i==0){
	  return {1.,0.};
	}else if(i==1){
	  return {0.,1.};
	}else if(i==2){
	  return {0.,0.};
	}
      }
    }
    Print(__PRETTY_FUNCTION__);
    abort();    
  };

  /// global to local: {s,m,n}
  VV_d M(const netFp& f, const V_d& ab){
    VV_d xyzs = f->getLocations();
    auto dXds0 = Dot(dNd_(0), xyzs);
    auto dXds1 = Dot(dNd_(1), xyzs);
    return {dXds0/Norm(dXds0)/*s*/,
	    dXds1/Norm(dXds1)/*m*/,
	    f->getNormal()/Norm(f->getNormal())/*n*/};
  };
  
  /// THIS ONLY USES ACTIVATED FACES
  map_P_d DphiDt(map_P_Vd& P_phiphin, const bool TorF=true){
    // GNUPLOT plot_;
    map_P_d ret;
    double norm;
    V_d Phi(3), u, v;    
    for(const auto& [p, phiphin]:P_phiphin){
      V_d dphidt;
      // p周りの各faceにおけるdphidtを計算
      for(const auto& f:p->getFaces(TorF)){
	auto ps = f->getPoints();
	Phi[0] = P_phiphin[ps[0]][0];
	Phi[1] = P_phiphin[ps[1]][0];
	Phi[2] = P_phiphin[ps[2]][0];

	v = Dot(Inverse(M(f,parameterize(ps, p))), { Dot(dNd_(0),Phi), Dot(dNd_(1),Phi), phiphin[1] }/*local velocity*/);// global velocity

	// plot_.SaveVectorData({{p->getX(),v}});
	
	dphidt.emplace_back(  pow(Norm(v),2)/2. - p->getX()[2] );
      }
      Print(dphidt, blue);
      ret[p] = Mean(dphidt);
    }

    // plot_.plot3d();
    // std::cin.ignore();
    
    return ret;
  };
  
  /// THIS ONLY USES ACTIVATED FACES
  ///   *---*
  ///  / \ / \
  /// *---@---*
  ///  \ / \ /
  ///   *---*
  map_P_Vd nablaPhi(map_P_Vd& P_phiphin, const bool TorF=true){
    map_P_Vd ret;
    V_d Phi(3), u, ab;

    // GNUPLOT plot;
    // VV_d vv;      
    // V_netFp tmpF;      

    for(auto& [p, phiphin]:P_phiphin){
      VV_d nablaPhi;

      for(const auto& f:p->getFaces(TorF)){
        // tmpF.push_back(f);
	auto ps = f->getPoints();
	Phi[0] = P_phiphin[ps[0]][0];
	Phi[1] = P_phiphin[ps[1]][0];
	Phi[2] = P_phiphin[ps[2]][0]; 
        nablaPhi.emplace_back(Dot(Inverse(M(f,parameterize(ps,p))), { Dot(dNd_(0), Phi), Dot(dNd_(1), Phi) , phiphin[1] }));
      }

      // vv.push_back(p->getX());
      // plot.SaveData(vv,{{"ps","2"},{"pt","7"},{"notitle",""}});
      // plot.SaveVectorData(getVectorData(tmpF),{{"arrowstyle","1"},{"notitle",""}});
      // plot.plot3d();
      // std::cin.ignore();
      // plot.Clear();
      
      ret[p] = Mean(nablaPhi);
      
    }
    
    return ret;
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

  V_d F(const VV_d& xyzs, const V_d& ab, const V_d& A){
    /// ローカルな面上の座標に沿って積分
    /// |detJ|の代わりに，z=0として，cross = Cross(Dot(dNd_(0), xyzs), Dot(dNd_(1), xyzs))のabs(cross[2])を使用，
    V_d r = Dot(N(ab),xyzs)-A, cross = Cross(Dot(dNd_(0), xyzs), Dot(dNd_(1), xyzs));
    double norm_r = Norm(r);
    return {1./norm_r * abs(cross[2]), - Dot(r,cross)/pow(norm_r,3)};
    /// calculate IG and IGn that satisfy
    /// IG *  phi_n  = (IGn - c * delta)  * phi
  };
      
  std::map<netP*, V_d> calc_P_IGIGn(const netFp& f, const netPp& origin, const VV_d& gwgw){
    V_d v, r, igign={0,0}, cross={0,0,0}, A = origin->getX();    
    V_netPp ps = f->getPoints();
    std::map<netP*, V_d> ret;
    V_d shapeF;
    ret[ps[0]] = {0,0};
    ret[ps[1]] = {0,0};
    ret[ps[2]] = {0,0};    
    
    VV_d xyzs_big = f->getLocations();

    if(f->intersectQ()){

      VV_d invTrX = Inverse(xyzs_big);
      VV_d xyzs_small, trT;
      
      for(const auto& routeP:f->getPointsCutFacesBehind()){
    	for(const auto& tri_ps:network::triangulate(routeP,f->getNormal())){
	  
    	  trT = Transpose(Dot( xyzs_small=obj3D::extractX(tri_ps), invTrX));

    	  for(const auto& abw:gwgw){    
	    igign = F(xyzs_small, abw, A) * abw[2];
    	    shapeF = Dot(trT, N(abw));
    	    ret[ps[0]] += igign * shapeF[0];
    	    ret[ps[1]] += igign * shapeF[1];
    	    ret[ps[2]] += igign * shapeF[2];
    	  }
	  
    	}
      }
      
    }
    else    
      {
	for(const auto& abw:gwgw){
	  igign = F(xyzs_big, abw, A) * abw[2];
	  shapeF = N(abw);
	  ret[ps[0]] += igign * shapeF[0];
	  ret[ps[1]] += igign * shapeF[1];
	  ret[ps[2]] += igign * shapeF[2];
	}
      }
    
    return ret;
  };
  
}
