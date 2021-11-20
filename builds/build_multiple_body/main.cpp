#define bem
#ifdef bem

#include "GNUPLOT.hpp"
#include "Network.hpp"
#include "bem.hpp"

int main()
{
  try
  {
    NetworkObj camel("../../obj/camel.obj", "camel");
    camel.scale(0.1);

    NetworkObj cow("../../obj/cow.obj", "cow");
    NetworkObj bunny("../../obj/bunny.obj", "bunny");
    bunny.scale(50);

    // NetworkObj cow("../../obj/galleon.obj", "cow");
    // NetworkObj bunny("../../obj/bunny.obj", "bunny");
    // cow.rotate(M_PI/2., {1.,0.,0.});
    // bunny.scale(2000);

    // NetworkX water(
    //     {25, 25},
    //     {50 / 3.2, 50 / 3.1, -1.},
    //     2,
    //     [](const V_d &xy) { return 0; },
    //     "water");
    
    cow.translate({RandomReal({-1.,1.}),RandomReal({-1.,1.}),RandomReal({-1.,1.})});
    bunny.translate({RandomReal({-1.,1.}),RandomReal({-1.,1.}),RandomReal({-1.,1.})});
    mk_vtu("./vtu/cow.vtu", cow.getFaces());
    mk_vtu("./vtu/bunny.vtu", bunny.getFaces());
    // mk_vtu("./vtu/camel.vtu", camel.getFaces());

 
    std::cin.ignore();

    BEM::CompGrid cow_bunny_camel({&cow, &bunny}, {{1000., 1000., 1000.}});

    // {
    //   std::cout << "------------------------------------" << std::endl;    
    //   XLoops loops(cow.Faces);
    //   std::cout << "success " << loops.success << std::endl;
    //   std::cout << "failure " << loops.failure << std::endl;
    // }
    // {
    //   std::cout << "------------------------------------" << std::endl;    
    //   XLoops loops(bunny.Faces);
    //   std::cout << "success " << loops.success << std::endl;
    //   std::cout << "failure " << loops.failure << std::endl;
    // }

    // display(takeIntxn(cow_bunny_camel.xnet->getLines()));
    // mk_vtu("./vtu/takeIntxn(cow_bunny_camel.xnet->getLines()).vtu", takeIntxn(cow_bunny_camel.xnet->getLines()));
    // std::cin.ignore();

    // display(takeIntxn(cow.getLines()));
    // mk_vtu("./vtu/takeIntxn(cow.getLines()).vtu", takeIntxn(cow.getLines()));
    // std::cin.ignore();

    // display(takeIntxn(bunny.getLines()));
    // mk_vtu("./vtu/takeIntxn(bunny.getLines()).vtu", takeIntxn(bunny.getLines()));
    // std::cin.ignore();    

    {
      Print("------------ cow ------------");
      display(&cow);
      searcherIntxn s;
      s.set(cow.getNearestFace({1000., 1000., -1000.}));
      s.addNetwork(cow_bunny_camel.xnet);
      s.search();
      mk_vtu("./vtu/cow_searcherIntxn_getObjects.vtu", s.getObjects());
      mk_vtu("./vtu/cow_searcherIntxn_getObjects_.vtu", s.getObjects_());
      mk_vtu("./vtu/cow_searcherIntxn_getObjects__.vtu", s.getObjects__());
      mk_vtu("./vtu/cow_searcherIntxn_getEnteredLines_.vtu", s.getEnteredLines_());      
      mk_vtu("./vtu/cow_searcherIntxn_intxn_.vtu", s.intxn);            
      mk_vtu("./vtu/cow_searcherIntxn_path_.vtu", s.path);                  
    }

    {
      Print("------------ bunny ------------");
      display(&bunny);
      searcherIntxn s;
      s.set(bunny.getNearestFace({1000., 1000., 0.}));
      s.addNetwork(cow_bunny_camel.xnet);
      s.search();
      mk_vtu("./vtu/bunny_searcherIntxn_getObjects.vtu", s.getObjects());
      mk_vtu("./vtu/bunny_searcherIntxn_getObjects_.vtu", s.getObjects_());
      mk_vtu("./vtu/bunny_searcherIntxn_getObjects__.vtu", s.getObjects__());
      mk_vtu("./vtu/bunny_searcherIntxn_getEnteredLines_.vtu", s.getEnteredLines_());
      mk_vtu("./vtu/bunny_searcherIntxn_intxn_.vtu", s.intxn);
      mk_vtu("./vtu/bunny_searcherIntxn_path_.vtu", s.path);      
    }

    Print(Green + "Enterを押して終了");
    std::cin.ignore();
  }
  catch (error_message e)
  {
    std::cout << e.what() << reset << std::endl;
    abort();
  };
};

#else

#include "GNUPLOT.hpp"
#include "Network.hpp"

int main()
{
  try
  {

    NetworkObj camel("./obj/camel.obj", "camel");
    camel.scale(0.1);
    V_d cardinal = {1000., 100., .2};
    NetworkObj cow("./obj/cow.obj", "cow");
    NetworkObj bunny("./obj/bunny.obj", "bunny");
    bunny.scale(50);
    Network camel_cow(camel, cow);
    Network camel_bunny(camel, bunny);
    Network cow_bunny(cow, bunny);

    mk_vtu("./vtu/bunny.vtu", bunny.Faces);
    mk_vtu("./vtu/cow.vtu", cow.Faces);
    mk_vtu("./vtu/camel.vtu", camel.Faces);
    mk_vtu("./vtu/camel_cow.vtu", camel_cow.getLines());
    mk_vtu("./vtu/camel_bunny.vtu", camel_bunny.getLines());
    mk_vtu("./vtu/cow_bunny.vtu", cow_bunny.getLines());

    Print(Green + "Enterを押して終了");
    std::cin.ignore();
  }
  catch (error_message e)
  {
    std::cout << e.what() << reset << std::endl;
    abort();
  };
};
#endif
