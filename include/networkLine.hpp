#pragma once

inline V_netFp networkLine::getSurfaces(networkFace *const f_excluded) const {
  V_netFp surfaces;
  surfaces.reserve(this->Faces.size());
  for (const auto &f : this->Faces)
    if (f != f_excluded && f->SurfaceQ())
      surfaces.emplace_back(f);
  return surfaces;
};

inline bool networkLine::SurfaceQ() const {
  return std::any_of(this->Faces.begin(), this->Faces.end(), [](const auto &f) { return f->SurfaceQ(); });
};

inline T2Tddd networkLine::getLocationsTuple() const { return {this->Point_A->getXtuple(), this->Point_B->getXtuple()}; };

inline void networkLine::setBoundsSingle() { CoordinateBounds::setBounds(getLocationsTuple()); };

inline bool networkLine::Replace(netP *oldP, netP *newP) {
  auto bool1 = this->replace(oldP, newP); // 1
  auto bool2 = oldP->erase(this);         // 2
  auto bool3 = newP->add(this);           // 3
  // このステップがdouble replace
  //  switchでないと，順番に意味のあるFaceではおかしくなるので注意
  if (bool1 && bool2 & bool3)
    return true;
  else
    return false;
};

inline networkLine::networkLine(Network *network_IN, netP *sPoint_IN, netP *ePoint_IN) : CoordinateBounds(Tddd{{0., 0., 0.}}), Faces(0), network(network_IN), Point_A(nullptr), Point_B(nullptr), Neumann(false), Dirichlet(false), CORNER(false) {
#ifdef DEM
  this->tension = 0.;
#endif
  network->Lines.emplace(this);
  set(sPoint_IN, ePoint_IN);
  sPoint_IN->add(this);
  ePoint_IN->add(this);
  setBoundsSingle();
  // setBounds();
};

inline bool networkLine::isIntxn() {
  //  return intxn; /*face-face intersection*/
  auto fs = this->Faces;
  if (fs.size() < 2)
    return false;
  else {
    for (auto i = 0; i < fs.size(); i++)
      for (auto j = i + 1; j < fs.size(); j++)
        if (fs[i]->getNetwork() != fs[j]->getNetwork())
          return true;
    return false;
  }
};

inline double networkLine::length() const { return Norm(this->Point_A->X - this->Point_B->X); };

inline Tddd networkLine::getNormal() const {
  Tddd ret = {0, 0, 0};
  for (const auto &f : this->Faces)
    ret += f->normal;
  return ret / (double)(this->Faces.size());
};

bool isLinkedDoubly(const netLp l, const netPp p) {
  try {
    auto [q0, q1] = l->getPoints();
    if (q0 && q0 == p) {
      for (const auto &m : p->getLines())
        if (m && m == l)
          return true;
    }
    if (q1 && q1 == p) {
      for (const auto &m : p->getLines())
        if (m && m == l)
          return true;
    }
    return false;
  } catch (std::exception &e) {
    std::cerr << e.what() << colorReset << std::endl;
    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
  };
};

bool isLinkedDoubly(const netLp l, const netFp f) {
  try {
    for (const auto &q : l->getFaces())
      if (q && q == f) {
        auto [l0, l1, l2] = f->getLines();
        if ((l0 && l0 == l) || (l1 && l1 == l) || (l2 && l2 == l))
          return true;
      }
    return false;
  } catch (std::exception &e) {
    std::cerr << e.what() << colorReset << std::endl;
    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
  };
};

bool isLinkedDoubly(const netPp p, const netLp l) { return isLinkedDoubly(l, p); };
bool isLinkedDoubly(const netFp f, const netLp l) { return isLinkedDoubly(l, f); };

//% ------------------------------------------------------ */
//%                          辺の分割                      */
//% ------------------------------------------------------ */

inline netPp networkLine::divide(const Tddd &midX) {

  /*
              p2
             /  \
           /      \
      l2 /          \  l1
       /      fA      \
     /                  \
p0, q1 ------ this ------ p1, q0
    \                   /
      \               /
        \     fB    /
      e1  \       / e2
            \   /
              q2

                p2
             / |   \
           /   |     \
      l2 /     |       \  l1
       /  fA  (LA)  (FA) \
     /         |           \
p0, q1--this---(P)---(LC)---- p1, q0
    \          |            /
      \       (LB)  (FB)  /
        \ fB   |        /
      e1  \    |      / e2
            \  |   /
               q2
*/

  /* ------------------------------------- */
  try {
    const auto adjacent_faces = this->getSurfaces();
    if (adjacent_faces.size() != 2) {
      std::cout << "adjacent_faces.size() != 2 but " << adjacent_faces.size() << std::endl;
      return nullptr;
    }
    auto fA = adjacent_faces[0]; // will be deleted
    auto fB = adjacent_faces[1]; // will be deleted
    std::cout << "fA: " << fA << ", fB: " << fB << std::endl;
    if (!fA || !fB) {
      std::cout << "fA == nullptr or fB == nullptr" << std::endl;
      return nullptr;
    }
    const auto [p0, this_, p1, l1, p2, l2] = fA->getPointsAndLines(this);
    const auto [q0, this__, q1, e1, q2, e2] = fB->getPointsAndLines(this);
    if (!p0 || !p1) {
      std::cout << "p0 == nullptr or p1 == nullptr" << std::endl;
      return nullptr;
    }
    /* ------------------------------------- */
    auto P = new networkPoint(this->network, (p0->X + p1->X) / 2.);
    dual_replace(this, p1, P);

    auto LC = new networkLine(this->network, P, p1);
    auto LA = new networkLine(this->network, P, p2);
    auto LB = new networkLine(this->network, P, q2);
    auto FA = new networkFace(this->network, P, LC, p1, l1, p2, LA);
    auto FB = new networkFace(this->network, P, LB, q2, e2, p1, LC);
    std::cout << "P: " << P << ", LA: " << LA << ", LB: " << LB << ", LC: " << LC << ", FA: " << FA << ", FB: " << FB << std::endl;

    dual_replace(fA, l1, LA, FA, nullptr);
    dual_replace(fB, e2, LB, FB, nullptr);
    FA->setPoints(P, LC, p1, l1, p2, LA);
    FB->setPoints(P, LB, q2, e2, p1, LC);
    dual_replace(fA, p1, P, FA);
    dual_replace(fB, p1, P, FB);
    fA->setPoints(P, LA, p2, l2, p0, this);
    fB->setPoints(P, this, p0, e1, q2, LB);
    //
    p2->add(FA);
    q2->add(FB);
    //
    p1->add(LC);
    p1->add(FA);
    p1->add(FB);
    p1->erase(this);
    p1->erase(fA);
    p1->erase(fB);
    P->Faces = {fA, fB, FA, FB};
    P->Lines = {LA, LB, LC, this};
    LA->Faces = {fA, FA};
    LB->Faces = {fB, FB};
    LC->Faces = {FA, FB};
    fA->syncPLPLPL();
    fB->syncPLPLPL();
    this->set(p0, P);

    auto p = P;
    {
      auto &phiphin = p->phiphin;
      auto &phiphin_t = p->phiphin_t;
      phiphin = Tdd{0, 0};
      phiphin_t = Tdd{0, 0};
      for (auto &q : p->getNeighborPointsOnSurfaces()) {
        phiphin += q->phiphin;
        phiphin_t += q->phiphin_t;
      }
      phiphin /= p->getNeighborPointsOnSurfaces().size();
      phiphin_t /= p->getNeighborPointsOnSurfaces().size();
    }

    return P;
  } catch (std::exception &e) {
    std::cerr << e.what() << colorReset << std::endl;
    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
  };
};

//@ ------------------------------------------------------ */
//@                          辺の削除                        */
//@ ------------------------------------------------------ */

inline bool networkLine::isMergeable() const {
  // isConsistent(this);
  // auto oldPoints = this->getPoints(); //最後にsetBoundsする必要がる
  //* ------------------------------------------------------ */
  auto [p0, p1] = this->getPoints();
  if (p0->getNeighbors().size() < 4 || p1->getNeighbors().size() < 4) {
    Print("今の所，mergeできない条件");
    return false;
  }

  auto AB = this->getFaces();
  if (AB.size() != 2) {
    Print("AB.size() != 2");
    return false;
  }
  auto A = AB[0];
  auto B = AB[1];
  auto Aps = A->getPoints(this);
  auto Bps = B->getPoints(this);
  //* ------------------------------------------------------ */
  int c = 0;
  std::ranges::for_each(Join(Bps, Aps), [&](const auto &p) {
    if (p->getLines().empty())
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "p->getLines().empty()");
  });

  /*ここが間違っていたconcaveの場合，拒否する*/
  if (isConcavePolygon({std::get<1>(Aps)->getXtuple(), std::get<2>(Aps)->getXtuple(), std::get<1>(Bps)->getXtuple(), std::get<2>(Bps)->getXtuple()})) {
    Print("Concave Polygonです．mergeしません．");
    return false;
  }
  //* ------------------------------------------------------ */
  auto Abl = A->getLineBack(this);
  auto Afl = A->getLineFront(this);
  auto Bbl = B->getLineBack(this);
  auto Bfl = B->getLineFront(this);

  auto a = (*Abl)(A);
  auto b = (*Bfl)(B);
  if (MemberQ(std::get<1>(Aps)->getFaces(), a) || MemberQ(std::get<1>(Aps)->getFaces(), b)) {
    std::cout << "面が潰れる条件なので，mergeしない" << std::endl;
    return false;
  }
  //* ------------------------------------------------------ */
  auto Abl_F_Except_A = TakeExcept(Abl->getFaces(), A);
  auto Bfl_F_Except_B = TakeExcept(Bfl->getFaces(), B);
  if (Abl_F_Except_A.empty()) {
    Print("Abl_F_Except_A.empty()");
    if (Bfl_F_Except_B.empty()) {
      Print("Bfl_F_Except_B.empty()");
      return false;
    }
    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
  }
  return true;
};

/* -------------------------------------------------------------------------- */

inline netPp networkLine::merge() {

  /*
                 p2
                /  \       fl1     /
   \          /      \            /
    \     l2/         (l1)       /
     \     /    (f0)     \      /
      \  /                 \   /
  p0, q1 <---->(this)<---->p1 q0 <----
       / \                 /   \
      /    \    (f1).     /     \
     /      \           /        \
    /       e1\      (e2)   fe2   \
   /            \   /              \
                 q2                 X <--------> p0 this case is not allowed
  */

  const auto adjacent_faces = this->getSurfaces();
  if (adjacent_faces.size() != 2) {
    std::cout << "adjacent_faces.size() != 2 but " << adjacent_faces.size() << std::endl;
    return nullptr;
  }
  auto f0 = adjacent_faces[0]; // will be deleted
  auto f1 = adjacent_faces[1]; // will be deleted

  std::cout << "f0: " << f0 << ", f1: " << f1 << std::endl;
  if (!f0 || !f1) {
    std::cout << "f0 == nullptr or f1 == nullptr" << std::endl;
    return nullptr;
  }

  const auto [p0, this_, p1, l1, p2, l2] = f0->getPointsAndLines(this);
  const auto [q0, this__, q1, e1, q2, e2] = f1->getPointsAndLines(this);

  if (this_ != this || this__ != this || p0 != q1 || p1 != q0)
    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "Topology error in flip");

  if (!p0 || !p1) {
    std::cout << "p0 == nullptr or p1 == nullptr" << std::endl;
    return nullptr;
  }

  if ((p0->getLines().size() + p1->getLines().size() - 3) <= 3 || p2->getLines().size() <= 3 || q2->getLines().size() <= 3) {
    // p2->getLines().size() <=2となるため，マージ不可
    // q2->getLines().size() <=2となるため，マージ不可
    std::cout << "p0 or p1 has less than 4 lines. merge is not allowed." << std::endl;
    return nullptr;
  }

  /* ------------------------------------- */
  auto meanX = p0->X + p1->X;
  meanX /= 2.;

  std::cout << "f0: " << f0 << ", f1: " << f1 << std::endl;
  std::cout << "p0: " << p0 << ", p1: " << p1 << ", p2: " << p2 << std::endl;
  std::cout << "p0->getLines(): " << p0->getLines() << std::endl;
  std::cout << "q0: " << q0 << ", q1: " << q1 << ", q2: " << q2 << std::endl;
  std::cout << "q0->getLines(): " << q0->getLines() << std::endl;
  std::cout << "l1: " << l1 << ", l1->Points: " << l1->getPoints() << ", l1->Faces: " << l1->Faces << std::endl;
  std::cout << "l2: " << l2 << ", l2->Points: " << l2->getPoints() << ", l2->Faces: " << l2->Faces << std::endl;
  std::cout << "e1: " << e1 << ", e1->Points: " << e1->getPoints() << ", e1->Faces: " << e1->Faces << std::endl;
  std::cout << "e2: " << e2 << ", e2->Points: " << e2->getPoints() << ", e2->Faces: " << e2->Faces << std::endl;
  std::cout << "this: " << this << ", this->Points: " << this->getPoints() << ", this->Faces: " << this->Faces << std::endl;

  auto faces_p1 = p1->getSurfaces();
  auto lines_p1 = p1->getLinesOnSurfaces();

  auto commonPoints = Intersection(p0->getNeighborPointsOnSurfaces(), p1->getNeighborPointsOnSurfaces());
  if (commonPoints.size() != 2) {
    std::cout << "inteersection : " << commonPoints << std::endl;
    std::cout << "p0 and p1 are directly connected by more than one line." << std::endl;
    return nullptr;
  }

  for (auto f : faces_p1)
    if (f != f0 && f != f1)
      dual_replace(f, p1, p0);

  for (auto l : lines_p1)
    if (l != this && l != l1 && l != e2)
      dual_replace(l, p1, p0);

  auto faces = l1->getSurfaces();
  if (faces.size() != 2)
    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "faces.size() !=2");
  auto fl1 = (faces[0] == f0 ? faces[1] : faces[0]);
  dual_replace(l2, f0, fl1, nullptr, l1);
  std::cout << "fl1: " << fl1 << ", fl1->Lines: " << fl1->Lines << std::endl;

  faces = e2->getSurfaces();
  if (faces.size() != 2)
    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "faces.size() !=2");
  auto fe2 = (faces[0] == f1 ? faces[1] : faces[0]);
  dual_replace(e1, f1, fe2, nullptr, e2);
  std::cout << "fe2: " << fe2 << ", fe2->Lines: " << fe2->Lines << std::endl;

  p2->erase(l1);
  q2->erase(e2);
  p0->erase(this);

  // この時点で，l1はまだ，f0,fl1を持っている．
  // この時点で，e2はまだ，f1,fe2を持っている．
  l1->erase(fl1); // delete l1でfl1が消されないように
  l1->erase(f0);
  e2->erase(fe2); // delete e2でfe2が消されないように
  e2->erase(f1);
  p0->erase(f0);
  p0->erase(f1);

  for (const auto &f : p0->getSurfaces())
    f->syncPLPLPL();

  p0->setXSingle(meanX);

  delete l1;
  delete e2;
  delete f0;
  delete f1;
  delete this;
  delete p1;

  auto checkLine = [&](netPp p0, netPp p1) {
    auto l0 = Line(p0, p1);
    auto l1 = Line(p1, p0);
    if (l0 && l1) {
      if (l0 != l1) {
        std::cout << "p0: " << p0 << ", p1: " << p1 << std::endl;
        std::cout << "l0: " << l0 << ", l1: " << l1 << std::endl;
        throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "Points are linked but by differrent lines");
      }
    } else
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "Points are not linked by any line");
  };
  std::cout << "p0->getLines(): " << p0->getLines() << std::endl;
  for (auto p : p0->getNeighbors())
    checkLine(p0, p);

  return p0;
}

/* -------------------------------------------------------------------------- */

inline netPp networkLine::mergeIfMergeable() {
  if (isMergeable())
    return this->merge();
  else
    return nullptr;
};

/*networkLine::divide_code*/
#include <csignal>
#include <iostream>
void signal_handler(int signal) {
  if (signal == SIGABRT) {
    std::cerr << "SIGABRT received: The program called abort()." << std::endl;
  }
}

/* -------------------------------------------------------------------------- */

inline bool networkLine::checkTopology() const {
  auto AB = this->getSurfaces();
  if (AB.size() != 2)
    return false;
  auto initial_faces = this->Faces;
  auto fA = AB[0];
  auto fB = AB[1];
  if (!fA || !fB || fA == fB)
    return false;
  auto [p0, this_, p1, l1, p2, l2] = fA->getPointsAndLines(this);
  auto [q0, this__, q1, e1, q2, e2] = fB->getPointsAndLines(this);
  if (this_ != this || this__ != this || p0 != q1 || p1 != q0) {
    std::cout << "Topology error " << std::endl;
    return false;
  }
  return true;
}

inline bool networkLine::flip(bool force = false) {
  try {

    auto AB = this->getSurfaces();
    if (AB.size() != 2)
      return false;

    auto initial_faces = this->Faces;

    auto fA = AB[0];
    auto fB = AB[1];

    if (!fA || !fB || fA == fB)
      return false;

    fA->syncPLPLPL();
    fB->syncPLPLPL();
    auto [p0, this_, p1, l1, p2, l2] = fA->getPointsAndLines(this);
    auto [q0, this__, q1, e1, q2, e2] = fB->getPointsAndLines(this);

    if (this_ != this || this__ != this || p0 != q1 || p1 != q0)
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "Topology error in flip");

    if (p0->getLines().size() <= 3 || p1->getLines().size() <= 3) {
      // p0->getLines().size() <=2となるため，フリップ不可
      // p1->getLines().size() <=2となるため，フリップ不可
      std::cout << "p0 or p1 has less than 4 lines. flip is not allowed." << std::endl;
      return false;
    }

    std::cout << "fA: " << fA << ", fB: " << fB << std::endl;
    std::cout << "p0: " << p0 << ", p1: " << p1 << ", p2: " << p2 << std::endl;
    std::cout << "p0->getLines(): " << p0->getLines() << std::endl;
    std::cout << "q0: " << q0 << ", q1: " << q1 << ", q2: " << q2 << std::endl;
    std::cout << "q0->getLines(): " << q0->getLines() << std::endl;
    std::cout << "l1: " << l1 << ", l1->Points: " << l1->getPoints() << ", l1->Faces: " << l1->Faces << std::endl;
    std::cout << "l2: " << l2 << ", l2->Points: " << l2->getPoints() << ", l2->Faces: " << l2->Faces << std::endl;
    std::cout << "e1: " << e1 << ", e1->Points: " << e1->getPoints() << ", e1->Faces: " << e1->Faces << std::endl;
    std::cout << "e2: " << e2 << ", e2->Points: " << e2->getPoints() << ", e2->Faces: " << e2->Faces << std::endl;
    std::cout << "this: " << this << ", this->Points: " << this->getPoints() << ", this->Faces: " << this->Faces << std::endl;

    if (!l1 || !l2 || !e1 || !e2) {
      std::cout << "Topology error: null lines in flip" << std::endl;
      return false;
    }

    auto n_before = 0.5 * (TriangleNormal(p0->X, p1->X, p2->X) + TriangleNormal(q0->X, q1->X, q2->X));
    auto n_after = 0.5 * (TriangleNormal(q2->X, p2->X, p0->X) + TriangleNormal(p2->X, q2->X, p1->X));

    if (force == false)
      if (!isFlat(n_before, n_after, M_PI / 180. * 20.))
        return false;

    /* ------------------------------------------------------------------- */

    auto initial_tetras = this->Tetras;
    auto this_faces = this->Faces;
    std::unordered_set<networkFace *> faces_L1_without_L0_plus_Surface = {fA, fB};
    for (const auto &t : initial_tetras)
      for (const auto &f : t->Faces)
        if (std::find(this_faces.begin(), this_faces.end(), f) == this_faces.end())
          faces_L1_without_L0_plus_Surface.emplace(f);

    // faces_L1_without_L0_plus_Surfaceのサイズが，4の場合，このflip操作は，四面体を潰すことになる．
    // 現在これには対応していないので，その場合はflipを行わない．
    // if (faces_L1_without_L0_plus_Surface.size() <= 4)
    //   return false;

    /* ------------------------------------------------------------------- */
    /*    点，線，面の接続関係の入れ替えを行う                             */
    /* ------------------------------------------------------------------- */

    /*
                   p2
                  /  \
                /      \
           l2 /          \  l1
            /      fA      \
          /                  \
    p0, q1 <----> this <----> p1, q0
         \                   /
           \               /
             \     fB    /
           e1  \       / e2
                 \   /
                   q2

                   p2
                  /A \
                /  |   \
           l2 /    |     \ l1
            /      |       \
          /        V         \
    p0, q1  fB   this  fA     p1, q0
         \         A         /
           \       |       /
             \     |     /
           e1  \   |   / e2
                 \ V /
                   q2
    */

    // #線の入れ替え
    std::cout << "線の点の入れ替え" << std::endl;
    this->set(p2, q2);
    p0->erase(this);
    p1->erase(this);
    p2->add(this);
    q2->add(this);

    // b% 線と点の接続関係については入れ替え終了

    // #線の面の入れ替え
    std::cout << "線の面の入れ替え" << std::endl;
    // l2 (p0側) を fA から fB へ移動
    l2->erase(fA);
    l2->add(fB);

    // e1 (p1側) を fB から fA へ移動 (e2ではなくe1を移動すべき)
    e2->erase(fB);
    e2->add(fA);

    // #面に点・線を入れる
    std::cout << "面の点・線の入れ替え" << std::endl;
    fA->setPoints(p2, this, q2, e2, p1, l1);
    fB->setPoints(q2, this, p2, l2, p0, e1);
    fA->syncPLPLPL();
    fB->syncPLPLPL();
    // b% 面と線の接続関係については入れ替え終了

    // #面の持つ点の入れ替え
    std::cout << "面の点の入れ替え" << std::endl;
    fA->setGeometricProperties(std::array<std::array<double, 3>, 3>{p2->X, q2->X, p1->X});
    fB->setGeometricProperties(std::array<std::array<double, 3>, 3>{q2->X, p2->X, p0->X});

    // #点の面の入れ替え
    std::cout << "点の面の入れ替え" << std::endl;
    p0->erase(fA);
    p1->erase(fB);
    p2->add(fB);
    q2->add(fA);

    // b% 点と面の接続関係については入れ替え終了
    std::cout << Green << "flip done" << colorReset << std::endl;
    return true;

  } catch (std::exception &e) {
    std::cerr << e.what() << colorReset << std::endl;
    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
    return false;
  };
};

/* -------------------------------------------------------------------------- */

inline bool networkLine::islegal() const {
  try {
    auto fs = this->getFaces();
    if (fs.size() < 2)
      return true;
    auto [Ap0, Ap1, Ap2] = ToX(fs[0]->getPoints(this));
    auto [Bp0, Bp1, Bp2] = ToX(fs[1]->getPoints(this));
    double sumangle = std::abs(VectorAngle(Bp1 - Bp2, Bp0 - Bp2)) + std::abs(VectorAngle(Ap1 - Ap2, Ap0 - Ap2));
    if (M_PI /*180.0.....1 若干大きい場合はOKとする*/ > sumangle)
      return true /*正*/;
    else
      return false /*不正*/;
  } catch (std::exception &e) {
    std::cerr << e.what() << colorReset << std::endl;
    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
  };
};
inline bool networkLine::flipIfIllegal() {
  if (!islegal() && !isIntxn())
    return (this->flip()); // flipはflipが成功した場合trueを返す．convexでない場合flipされない場合がある
  else
    return false;
};

inline bool networkLine::isAdjacentFacesFlat(const double minangle = M_PI / 180.) const {
  auto fs = this->getFaces();
  auto [p0, p1, p2] = fs[0]->getPoints();
  auto [P0, P1, P2] = fs[1]->getPoints();
  if (fs.size() != 2)
    return false;
  else if (isFlat(Cross(p1->X - p0->X, p2->X - p0->X), Cross(P1->X - P0->X, P2->X - P0->X), minangle))
    return true;
  else
    return false;
};

/*DOC_EXTRACT flip

### `flip`可能かどうかの判定

\ref{canFlip}{`canFlip`}でフリップ可能かどうかを判定する．直感的に次のような条件の場合，境界面が崩れるため，フリップさせたくない．

* フリップ前後で，辺に隣接する面の面積の和が大きく変化する場合，フリップさせない
* フリップ前後で，辺に隣接する面の法線ベクトルが大きく変換する場合，フリップさせない

しかし，これの判定において必要となる計算：三角形の内角や法線方向，ベクトルの成す角度の計算は，精確に判定できない領域があるようだ．
なので，その領域をおおよそ実験的に調べて，まずはその領域に入らせない条件を設ける（信頼できる三角形）．
次のような三角形は信頼しない：

* 三角形の内角が小さすぎる，または大きすぎる場合
* 内角の和が$\pi$にならない場合

信頼できる三角形の判定には，\ref{isValidTriangle}{`isValidTriangle`}を用いる．

*/

// \label{canFlip}
inline bool networkLine::canFlip(const double acceptable_n_diff_before_after = 1E-3 * M_PI / 180.) const {
  try {
    auto f_and_F = this->getFaces();
    if (f_and_F.size() != 2)
      return false;
    auto [f0, f1, f2] = f_and_F[0]->getPoints(this);
    auto [F0, F1, F2] = f_and_F[1]->getPoints(this);
    auto tri0_now = T3Tddd{f0->X, f1->X, f2->X};
    auto tri1_now = T3Tddd{F0->X, F1->X, F2->X};

    auto tri0 = T3Tddd{f0->X, F2->X, f2->X};
    auto tri1 = T3Tddd{F0->X, f2->X, F2->X};

    if (!isValidTriangle(tri0, M_PI / 180.) || !isValidTriangle(tri1, M_PI / 180.))
      return false;

    // if (isFlat(tri0[1] - tri0[0], tri0[2] - tri0[0], 1E-1))
    //    return false;
    // if (isFlat(tri0[2] - tri0[1], tri0[0] - tri0[1], 1E-1))
    //    return false;
    // if (isFlat(tri1[1] - tri1[0], tri1[2] - tri1[0], 1E-1))
    //    return false;
    // if (isFlat(tri1[2] - tri1[1], tri1[0] - tri1[1], 1E-1))
    //    return false;
    //$ large difference of normal vector after and before flip
    if (!isFlat(tri0[0], tri0[1], tri0[2], tri0_now[0], tri0_now[1], tri0_now[2], acceptable_n_diff_before_after) || !isFlat(tri1[0], tri1[1], tri1[2], tri1_now[0], tri1_now[1], tri1_now[2], acceptable_n_diff_before_after))
      return false;

    if (!isFlat(TriangleNormal(tri0), TriangleNormal(tri0_now), acceptable_n_diff_before_after) || !isFlat(TriangleNormal(tri1), TriangleNormal(tri1_now), acceptable_n_diff_before_after))
      return false;

    //$ area conservation
    return true; // TriangleArea(tri0) + TriangleArea(tri1) == TriangleArea(tri0_now) + TriangleArea(tri1_now);

  } catch (const std::exception &e) {
    std::cerr << e.what() << std::endl;
    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
  };
}

inline bool networkLine::flipIfBetter(const double n_diff_tagert_face, const double acceptable_n_diff_before_after, const int min_n) {
  try {
    //@ flipを実行するには面の法線方向が成す全ての角度かこれよりも小さくなければならない
    //@ フリップ前後の両方で不正な辺と判定された場合，
    //@ 線の数と面の面積の差をチェックし，差が少ない方を選択する．
    if (!canFlip(acceptable_n_diff_before_after))
      return false;
    // else if (this->isAdjacentFacesFlat(n_diff_tagert_face /*ここで引っかかってしまいフリップされないことがよくある*/) && !islegal() && !isIntxn()) {
    else if (this->isAdjacentFacesFlat(n_diff_tagert_face /*ここで引っかかってしまいフリップされないことがよくある*/) && !islegal()) {
      auto [p0, p1] = this->getPoints();

      auto f0f1_ = this->getFaces();
      int s0 = p0->getLines().size();
      int s1 = p1->getLines().size();
      auto p2 = f0f1_[0]->getPointOpposite(this);
      auto p3 = f0f1_[1]->getPointOpposite(this);
      int s2 = p2->getLines().size();
      int s3 = p3->getLines().size();
      // if (s0 > 3 || s1 > 3 || s2 > 3 || s3 > 3)
      //    if (s0 - 1 < 4 || s1 - 1 < 4 || s2 + 1 < 4 || s3 + 1 < 4)
      //       return false;  // 3以下はつくらない

      /* -------------------------------------------------------------------------- */
      // auto [p0, p2] = this->getPoints();
      auto f0f1 = this->getFaces();
      auto [f0pb, f0pf, f0po] = f0f1[0]->getPoints(this);
      auto [f1pb, f1pf, f1po] = f0f1[1]->getPoints(this);
      // 面積の減少
      // double diffAinit = std::abs(TriangleArea(T3Tddd{f0pb->getXtuple(), f0pf->getXtuple(), f0po->getXtuple()}) - TriangleArea(T3Tddd{f1pb->getXtuple(), f1pf->getXtuple(), f1po->getXtuple()}));
      // double diffAnext = std::abs(TriangleArea(T3Tddd{f0pb->getXtuple(), f1po->getXtuple(), f0po->getXtuple()}) - TriangleArea(T3Tddd{f0pf->getXtuple(), f0po->getXtuple(), f1po->getXtuple()}));
      // double diff = (diffAnext - diffAinit) / (f0f1[0]->getArea() + f0f1[1]->getArea());
      /*
      //     f0po *------* f0pf,f1pb
      //          |   /  |
      //f0pb,f1pf *------* f1po
      */

      double min_init = std::min(Min(TriangleAngles(T3Tddd{f0pb->X, f0pf->X, f0po->X})), Min(TriangleAngles(T3Tddd{f1pb->X, f1pf->X, f1po->X})));
      double min_later = std::min(Min(TriangleAngles(T3Tddd{f0pb->X, f1po->X, f0po->X})), Min(TriangleAngles(T3Tddd{f1pb->X, f0po->X, f1po->X})));

      int next_s0 = s0 - 1;
      int next_s1 = s1 - 1;
      int next_s2 = s2 + 1;
      int next_s3 = s3 + 1;
      /*
      @ <-------------(-c)------------ (c) --------------
      @ <---- flip -----|--- topology --|----- none -----
      */
      // if (min_init < min_later)
      // if (std::abs(min_init - min_later) < M_PI / 180. * 5)
      // {
      // 	this->flipIfTopologicalyBetter(min_degree_to_flat);
      // 	return true;
      // }
      // else
      if (min_init <= min_later && (next_s0 >= min_n || p0->CORNER /*min_nよりも小さくても，角ならOK*/) && (next_s1 >= min_n || p1->CORNER /*min_nよりも小さくても，角ならOK*/) && (next_s2 >= min_n || p2->CORNER /*min_nよりも小さくても，角ならOK*/) && (next_s3 >= min_n || p3->CORNER /*min_nよりも小さくても，角ならOK*/)) {
        this->flip();
        return true;
      } else
        return false;
    } else
      return false;
  } catch (std::exception &e) {
    std::cerr << e.what() << colorReset << std::endl;
    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
  };
};

inline bool networkLine::flipIfTopologicallyBetter(const double n_diff_tagert_face, const double acceptable_n_diff_before_after, const int s_meanIN) {
  try {
    // Check if the flip is allowed based on the minimum angle of the triangle after flipping
    if (!canFlip(acceptable_n_diff_before_after))
      return false;
    if (this->isAdjacentFacesFlat(n_diff_tagert_face) && !isIntxn()) {
      auto [p0, p1] = this->getPoints();
      auto f0f1 = this->getFaces();
      int s0 = p0->getLines().size();
      int s1 = p1->getLines().size();
      auto p2 = f0f1[0]->getPointOpposite(this);
      auto p3 = f0f1[1]->getPointOpposite(this);
      int s2 = p2->getLines().size();
      int s3 = p3->getLines().size();
      double s_mean = s_meanIN;
      double v_init = Norm(std::array<int, 4>{s0, s1, s2, s3} - s_mean);
      double v_next = Norm(std::array<int, 4>{s0 - 1, s1 - 1, s2 + 1, s3 + 1} - s_mean);
      // ただし，例えばs0が角であった場合，角が多くの線を持つことは問題ないため，考慮に入れない．つまり辺の数は変わっても変わらないものとして，v_nextを考える．
      // v_next = Norm(std::array<int, 4>{p0->CORNER || p0->isMultipleNode ? s0 : s0 - 1, p1->CORNER || p1->isMultipleNode ? s1 : s1 - 1,
      //                                  p2->CORNER || p2->isMultipleNode ? s2 : s2 + 1, p3->CORNER || p3->isMultipleNode ? s3 : s3 + 1} -
      //               s_mean);
      if (v_init > v_next || s0 >= 8 || s1 >= 8 || s2 <= 4 || s3 <= 4) {
        this->flip();
        return true;
      } else
        return false;
    } else
      return false;
  } catch (std::exception &e) {
    std::cerr << e.what() << std::endl;
    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
  }
};

inline void networkLine::divideIfIllegal() {
  if (!islegal())
    this->divide();
};
