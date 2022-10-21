std::vector<networkFace*> getIntersectingFaces(const std::vector<networkFace*>& faces,
					       std::vector<double> a, std::vector<double> b){
  std::vector<networkFace*> ret;
  std::vector<std::vector<double>> locs;
  double ratio = 1E-6;
  for(const auto& face:faces){
    locs = (1.+ratio)*(face->getLocations())-ratio*(face->getMeanLocation());//少し広げた頂点      
    if(isIntersectingSurface(locs, std::vector<std::vector<double>>{a,b}) == 3){
      ret.emplace_back(face);
    }
  }
  return ret;  
};
class test3{
public:
  test3(){   
    NetworkObj obj("./obj/teddy.obj");
    GNUPLOT plot;
    plot.Set({{"style","arrow 1 nohead lc 1 lw .1"}});

    std::vector<std::vector<std::vector<double>>> objvec;
    for(const auto& line:obj.getLines()){
      auto ab = line->getLocations();
      objvec.push_back({ab[0],ab[1]-ab[0]});
    }
  
    std::vector<double> a={0.,20.,0.}, b={0.,0.,0};
    auto face = getIntersectingFaces(obj.Faces,a,b)[0];
  
    std::vector<double> hitP = pOnSurface(face->getLocations(),{a,b});
  
    plot.SaveData(hitP,{{"lw","5"},{"lt","5"},{"lc","6"}});  

    for(auto k=1; k<100; k++){
      plot.Clear();
      plot.SaveVectorData({{a,b-a}});
      plot.SaveData(hitP,{{"lw","5"},{"lt","5"},{"lc","6"}});
      plot.SaveVectorData(objvec,{{"arrowstyle","1"}});

      plot.Set({{"cbrange","[0:1000]"}});    
      auto Rmat = RotationMatrix(2.*1*M_PI/100.,face->getNormal());
      std::vector<double> vvv = 500.*Dot(Rmat, face->getLocations()[0]-hitP);
      //    networkLine* line = face->getCrossLine(hitP, vvv);
      pathInfo info = face->getPathInfo(hitP, vvv);

      int l=0;
    
      while(l<k){
	info = info.face[1]->getPathInfo(info);
	plot.SaveVectorData({{info.xyz[0],info.xyz[2]-info.xyz[0],{double(l++)}}},{{"lw","5"},{"lc","palette"}});
	if(info.face[1]==NULL)
	  break;
      }    

      if(k%10==0){
	plot.Plot3D_All();
	//std::cin.ignore();
      }
    }
  };
  ~test3(){Print(CHECK+" test3",Red);};
};
