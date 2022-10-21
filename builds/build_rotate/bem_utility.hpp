
// できるだけカプセル化することが，今後の変更を最小に抑えてくれる．
// この観点から幾何学関連の操作（面の抽出や三角分割など）は，bemには持ち込まない．
#ifndef INCL_BEM
#define INCL_BEM
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
  using V_netFp = std::vector<netF*>;  
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
    for(auto i=0; i<points.size(); i++){
      if(p==points[i]){
	switch(i){
        case 0:
	  return {1.,0.,0.};
        case 1:
	  return {0.,1.,0.};
        case 2:
	  return {0.,0.,1.};
        default:
          Print(__PRETTY_FUNCTION__);
          abort();          
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
  VV_d M(const VV_d& xyzs, const V_d& ab){
    auto dXds0 = Dot(dNd_(0), xyzs);
    auto dXds1 = Dot(dNd_(1), xyzs);
    auto cross = Cross(dXds0,dXds1);
    return {dXds0/Norm(dXds0)/*s*/,
	    dXds1/Norm(dXds1)/*m*/,
	    cross/Norm(cross)/*n*/};
  };
  
  /// THIS ONLY USES ACTIVATED FACES
  map_P_d DphiDt(map_P_Vd& P_phiphin, const bool TorF=true){
    Print(__PRETTY_FUNCTION__,Blue);
    //    GNUPLOT plot_;
    map_P_d ret;
    double norm;
    V_d Phi(3), v_L, v_G;
    VVV_d vvv_L, vvv_G;
    V_netFp fs;
    
    for(const auto& [p, phiphin]:P_phiphin){
      V_d dphidt;
      // p周りの各faceにおけるdphidtを計算
      //      Print(p,red);
      for(const auto& f:p->getFaces_intersectQ(false)){
	fs.push_back(f);

	// plot_.SaveVectorData(getVectorData(fs),{{"lc","'black'"},{"title","vvv_L"}});      
	// plot_.plot3d();
	// plot_.Clear();	
	//Print(f,blue);

	auto ps = f->getPoints();//このポイントがP_phiphinにありますか?
	// Print("    p---*    ex.  ps[0]---ps[2]",Red);
	// Print(" ps \\ /            \\   /",Red);
	// Print("      *               ps[1]",Red);

	// Print("ps",Red);	
	// Print(ps,red);

	// Print("P_phiphin[ps[0]]",Red);
	// Print(P_phiphin[ps[0]],Red);
	// Print("P_phiphin[ps[1]]",Red);
	// Print(P_phiphin[ps[1]],Red);
	// Print("P_phiphin[ps[2]]",Red);
	// Print(P_phiphin[ps[2]],Red);
		

	Phi[0] = P_phiphin[ps[0]][0];
	Phi[1] = P_phiphin[ps[1]][0];
	Phi[2] = P_phiphin[ps[2]][0];

	v_L = { Dot(dNd_(0),Phi), Dot(dNd_(1),Phi), phiphin[1] };
	v_G = Dot(Inverse(M({ps[0]->getX(),ps[1]->getX(),ps[2]->getX()},parameterize(ps, p))), v_L/*local velocity*/);// global velocity
	
	vvv_L.push_back({p->getX(),v_L});
	vvv_G.push_back({p->getX(),v_G});	
	norm = Norm(v_G);
	dphidt.emplace_back(  norm*norm/2. - p->getX()[2] );
      }
      ret[p] = Mean(dphidt);
    }
    //    std::cin.ignore();
	
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
    V_d Phi(3), u, ab, v_L, v_G;

    // GNUPLOT plot;
    // VV_d vv;      
    // V_netFp tmpF;      

    for(auto& [p, phiphin]:P_phiphin){
      VV_d nablaPhi;
      for(const auto& f:p->getFaces_intersectQ(false)){
        // tmpF.push_back(f);
	auto ps = f->getPoints();
	Phi[0] = P_phiphin[ps[0]][0];
	Phi[1] = P_phiphin[ps[1]][0];
	Phi[2] = P_phiphin[ps[2]][0];
	v_L = { Dot(dNd_(0),Phi), Dot(dNd_(1),Phi), phiphin[1] };	
        nablaPhi.emplace_back(Dot(Inverse(M(f,parameterize(ps,p))), v_L));
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





  map_P_Vd nablaPhi(map_P_Vd& P_phiphin, const V_netFp& Faces ,const bool TorF=true){
    map_P_Vd ret;
    V_d Phi(3), u, ab, v_L, v_G;

    // GNUPLOT plot;
    // VV_d vv;      
    // V_netFp tmpF;      

    for(auto& [p, phiphin]:P_phiphin){
      VV_d nablaPhi;
      for(const auto& f:p->getFaces_intersectQ(false))
        if(MemberQ(Faces,f)){
          // tmpF.push_back(f);
          auto ps = f->getPoints();
          Phi[0] = P_phiphin[ps[0]][0];
          Phi[1] = P_phiphin[ps[1]][0];
          Phi[2] = P_phiphin[ps[2]][0];
          v_L = { Dot(dNd_(0),Phi), Dot(dNd_(1),Phi), phiphin[1] };	
          nablaPhi.emplace_back(Dot(Inverse(M(f,parameterize(ps,p))), v_L));
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
    
    V_d r = Dot(N(ab),xyzs)-A;
    V_d cross = Cross(xyzs[0]-xyzs[2], xyzs[1]-xyzs[2]);
    double nr = Norm(r);
    return {abs(cross[2])/nr, - Dot(r,cross)/(nr*nr*nr)};
    /// calculate IG and IGn that satisfy
    /// IG *  phi_n  = (IGn - c * delta)  * phi
  };

    //  static constexpr VV_d gwgw = GaussianQuadratureWeights(5, 0., 1.);
    
  std::map<netP*, V_d> calc_P_IGIGn(const netFp& f, const netPp& origin, const VV_d& gwgw){    
    V_d igign={0,0}, cross={0,0,0}, A = origin->getX(), shapeF;
    V_netPp ps = f->getPoints();
    VV_d xyz = f->getLocations();

    if(false && MemberQ(ps,origin)){
      int n = gwgw.size();
      V_d param = parameterize(ps, origin);
      VV_d gw0_ = GaussianQuadratureWeights(n, 0., 1., param[0], 2.);
      VV_d gw1_ = GaussianQuadratureWeights(n, 0., 1., param[1], 2.);
      // Print("singular GaussianQuadratureWeights",red);
      //std::cout << "parameterize : " << parameterize(ps, origin)<< std::endl;
      // std::cout << "gw0_ : " << gw0_ << std::endl;
      // std::cout << "gw1_ : " << gw1_ << std::endl;      
      V_d v;
      double a, b;
      VV_d gwgw_;
      for(const auto& tw0:gw0_){
        a = tw0[0];
        for(const auto& tw1:gw1_){
          v = tw1*(1-a);
          b = v[0];
          gwgw_.push_back({a, b , tw0[1] * v[1]});
        }
      }

      V_d tmp0={0,0}, tmp1={0,0}, tmp2={0,0};
      for(const auto& abw:gwgw_){
        igign = F(xyz, abw, A) * abw[2];
        shapeF = N(abw);
        tmp0 += igign * shapeF[0];
        tmp1 += igign * shapeF[1];
        tmp2 += igign * shapeF[2];
      }
      return {{ps[0],tmp0},{ps[1],tmp1},{ps[2],tmp2}};      
    }else{
      // for(auto i=0; i<3; i++)
      //   if(ps[i]!=origin){
      //     if(Norm(ps[i]->getX()-A)<1E-1)
      //       return ret;
      //   }

      V_d tmp0={0,0}, tmp1={0,0}, tmp2={0,0};
      for(const auto& abw:gwgw){
        igign = F(xyz, abw, A) * abw[2];
        shapeF = N(abw);
        tmp0 += igign * shapeF[0];
        tmp1 += igign * shapeF[1];
        tmp2 += igign * shapeF[2];
      }
      return {{ps[0],tmp0},{ps[1],tmp1},{ps[2],tmp2}};      
    }
  };

  // std::map<netP*, V_d> calc_P_IGIGn(const netFp& f, const netPp& origin, const VV_d& gwgw){
  //   V_d v, r, igign={0,0}, cross={0,0,0}, A = origin->getX();    
  //   V_netPp ps = f->getPoints();
  //   std::map<netP*, V_d> ret;
  //   V_d shapeF;
  //   ret[ps[0]] = {0,0};
  //   ret[ps[1]] = {0,0};
  //   ret[ps[2]] = {0,0};    
    
  //   VV_d xyzs_big = f->getLocations();

  //   if(f->intersectQ()){

  //     VV_d invTrX = Inverse(xyzs_big);
  //     VV_d xyzs_small, trT;
  //     V_d tmp0={0,0}, tmp1={0,0}, tmp2={0,0};
  //     for(const auto& routeP:f->getPointsCutFacesBehind()){
  //   	for(const auto& tri_ps:network::triangulate(routeP,f->getNormal())){
	  
  //   	  trT = Transpose(Dot( xyzs_small=obj3D::extractX(tri_ps), invTrX));
  // 	  tmp0={0,0};
  // 	  tmp1={0,0};
  // 	  tmp2={0,0};
  //   	  for(const auto& abw:gwgw){    
  // 	    igign = F(xyzs_small, abw, A) * abw[2];
  //   	    shapeF = Dot(trT, N(abw));
  //   	    tmp0 += igign * shapeF[0];
  //   	    tmp1 += igign * shapeF[1];
  //   	    tmp2 += igign * shapeF[2];
  //   	  }
  // 	  ret[ps[0]] += tmp0;
  // 	  ret[ps[1]] += tmp1;
  // 	  ret[ps[2]] += tmp2;	  
  //   	}
  //     }
      
  //   }
  //   else    
  //     {
  // 	V_d tmp0={0,0}, tmp1={0,0}, tmp2={0,0};
  // 	for(const auto& abw:gwgw){
  // 	  igign = F(xyzs_big, abw, A) * abw[2];
  // 	  shapeF = N(abw);
  // 	  tmp0 += igign * shapeF[0];
  // 	  tmp1 += igign * shapeF[1];
  // 	  tmp2 += igign * shapeF[2];
  // 	}
  // 	ret[ps[0]] = tmp0;
  // 	ret[ps[1]] = tmp1;
  // 	ret[ps[2]] = tmp2;		
  //     }
    
  //   return ret;
  // };
  
}
#endif
