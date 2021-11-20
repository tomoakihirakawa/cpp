class test2{
public:
  test2(){
    NetworkObj obj("./obj/tank.obj");
    obj.rotate(2*M_PI, {1.,0.,0.});
    NetworkW water({1,1},{1.17,1.15,1/2.},.3);
    GNUPLOT plot;
    //=========== plot ===========
    plot.Set({{"style","arrow 1 nohead lc \"blue\" lw 1"}});
    plot.Set({{"style","arrow 2 nohead lc 2 lw 1"}});
    plot.Set({{"style","arrow 3 nohead lc \"red\" lw 2"}});
    for(auto i=0; i<8; i++){//ここが増えるとより多くcrossをチェックするため，crossの周辺だけをより細かくする
      plot.Set({{"key",""},{"title","\"divide "+std::to_string(i)+"\""}});
      Network corssNetwork(obj,water);
    
      plot.SaveVectorData(getVectorData(water.Faces),{{"arrowstyle","1"},{"title","water"}});      
      plot.SaveVectorData(getVectorData(obj.getLines()),{{"arrowstyle","2"},{"title","obj"}});          
      plot.SaveVectorData(getVectorData(corssNetwork.getLines()),{{"arrowstyle","3"},{"title","corssNetwork"}});

      plot.SaveData(corssNetwork.getLocations(),{{"w","p"},{"lc","\"magenta\""},{"pt","7"},{"title","cross points"}});    
      plot.SaveData(water.getMeanLocation(), {{"pt","20"},{"lc","0"},{"ps","3"},{"title","中心点"}});    
      plot.SaveData(water.getNearestPoints(water.getMeanLocation())->xyz,{{"pt","20"},{"lc","10"},{"ps","3"},{"title","中央付近の点"}});
    
      plot.Plot3D_All();
      //std::cin.ignore();
      plot.Clear();
      //======== divide =======
      for(const auto& line:water.getLines())
	if(line->penetrateQ())
	  water.divide(longerLine(line));

      //for(auto dummy:std::vector<int>{1,1,1,1}){
      //   for(const auto& line:obj.getLines())
      // 	if(!line->crossinfos.empty())
      // 	  obj.divide(line->longerLine());
      // }
      //-----------------------
      obj.clearXPoints();
      water.clearXPoints();    
      obj.setXPoints(water);
      water.setXPoints(obj);    
      //-----------------------
      // water.setLinkedPointsIncludeCross(water.getNearestPoints(water.getMeanLocation()),false);
      // for(const auto& p:water.Points)
      //   if(p->status)
      // 	p->Delete();

      bool found;
      network::setStatus(water.getNeighborsIncludeCross(water.getNearestPoints(water.getMeanLocation())),false);
      for(const auto& f:water.Faces){
	found = false;
	for(const auto& l:f->Lines)
	  if(!l->status){
	    found = true;
	    break;
	  }
	if(!found && !f->penetratedQ())
	  f->Delete();
      }
    
    }
  };
  ~test2(){Print(CHECK+" test2",Red);};
};
