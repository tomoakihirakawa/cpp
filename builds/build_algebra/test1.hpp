class test1{
public:
  test1(){
    Print("**************** test1 **************** ",Red);
    NetworkObj obj("./obj/tank.obj");
    NetworkW water({5,5},{1.12,1.15,1/2.},.1);  
    GNUPLOT plot;
    //=========== plot ===========
    plot.Set({{"xrange","[-.8:.8]"},{"yrange","[-.8:.8]"},{"zrange","[0:1.2]"}});
    plot.Set({{"style","arrow 1 nohead lc \"blue\" lw 1"}});
    plot.Set({{"style","arrow 2 nohead lc 2 lw 1"}});
    plot.Set({{"style","arrow 3 nohead lc \"red\" lw 1"}});
    plot.Set({{"style","arrow 4 nohead lc \"blue\" lw 1"}});
    plot.Set({{"style","arrow 5 nohead lc \"red\" lw 1"}});    
    for(auto i=0; i<1000; i++){//ここが増えるとより多くcrossをチェックするため，crossの周辺だけをより細かくする
    
      plot.Set({{"key",""},{"title","\"divide "+std::to_string(i)+"\""}});
      Network corssNetwork(obj,water);    
      int counter=0;

      // {
      //   std::vector<std::vector<std::vector<double>>> vecvec; 
      //   for(const auto& p:water.getPoints())
      // 	  vecvec.push_back({p->xyz, p->getNormal()});
      //   plot.SaveVectorData(vecvec,{{"notitle",""}});      
      // }    

      {
        std::vector<std::vector<std::vector<double>>> vecvec; 
        for(const auto& f:water.Faces)
      	  vecvec.push_back({f->getMeanLocation(), f->getNormal()});
        plot.SaveVectorData(vecvec,{{"notitle",""}});      
      }

      for(const auto& p:water.getPoints()){
	double x = p->xyz[0];
	double y = p->xyz[1];
	p->xyz += {0,0,cos(M_PI*(x*y)+i/4.)/30.};
      }
    
      //      plot.SaveVectorData(getVectorData(water.getFaces()),{{"arrowstyle","1"},{"title","water"}});      
      plot.SaveVectorData(getVectorData(obj.getLines()),{{"arrowstyle","2"},{"title","obj"}});          
      plot.SaveVectorData(getVectorData(corssNetwork.getLines()),{{"arrowstyle","3"},{"title","corssNetwork"}});
      //      plot.SaveData(corssNetwork.getLocations(),{{"w","p"},{"lc","\"magenta\""},{"pt","7"},{"title","cross points"}});    

      network::setStatus(water.getNeighborsIncludeCross(water.getNearestPoints(water.getMeanLocation())),false);
      plot.SaveData(getLocations(water.getPoints(true)),
		    {{"w","p"},{"lc","\"black\""},{"pt","7"},{"title","cross points"}});
    

      {
	VVV_double YYY;
	network::setStatus(water.Points,false);
	auto startP = network::takeNearest(obj3D::takeInsideOfBounds(water.Points,{{-0.55,0.55},{-0.55,0.55},{0.1,1.}}),{0,0,0});      
	auto ps = startP->getNeighbors_AvoidCross_recursive();		      
	network::setStatus(ps,true);
	for(const auto& f:network::takeIfIntersect(water.Faces)){
	  for(const auto& routeP:f->getPointsCutFacesBehind()){
	    if(true){
	      for(const auto& ind: triangulate(routeP,f->getNormal()) ){
		YYY.push_back({routeP[ind[0]]->getX(),routeP[ind[1]]->getX()-routeP[ind[0]]->getX()});
		YYY.push_back({routeP[ind[1]]->getX(),routeP[ind[2]]->getX()-routeP[ind[1]]->getX()});
		YYY.push_back({routeP[ind[2]]->getX(),routeP[ind[0]]->getX()-routeP[ind[2]]->getX()});		
	      }
	    }
	  }
	}
	plot.SaveVectorData(YYY,{{"arrowstyle","4"},{"notitle",""}});      
      }
      {
	VVV_double YYY;
	network::setStatus(obj.Points,false);
	auto startP = network::takeNearest(obj3D::takeInsideOfBounds(obj.Points,{{-0.55,0.55},{-0.55,0.55},{0.1,1.}}),{0,0,0});      
	auto ps = startP->getNeighbors_AvoidCross_recursive();		      
	network::setStatus(ps,true);
	for(const auto& f:network::takeIfIntersect(obj.Faces)){
	  for(const auto& routeP:f->getPointsCutFacesBehind()){
	    if(true){
	      for(const auto& ind: triangulate(routeP,f->getNormal()) ){
		YYY.push_back({routeP[ind[0]]->getX(),routeP[ind[1]]->getX()-routeP[ind[0]]->getX()});
		YYY.push_back({routeP[ind[1]]->getX(),routeP[ind[2]]->getX()-routeP[ind[1]]->getX()});
		YYY.push_back({routeP[ind[2]]->getX(),routeP[ind[0]]->getX()-routeP[ind[2]]->getX()});		
	      }
	    }
	  }
	}
	plot.SaveVectorData(YYY,{{"arrowstyle","5"},{"notitle",""}});      
      }
      

      plot.Plot3D_All();
      //std::cin.ignore();
      plot.Clear();
    }
  };
  ~test1(){Print(CHECK+" test1",Red);};
};
