#ifndef BEM_calculateVelocities_H
#define BEM_calculateVelocities_H

#include "BEM_utilities.hpp"
#include "Network.hpp"

//$ -------------------------------------------------------------------------- */
//$                         calculateVecToSurface                              */
//$ -------------------------------------------------------------------------- */

Tddd RK_without_Ubuff(const networkPoint *p) { return p->RK_X.getX(p->U_update_BEM); };
T3Tddd RK_without_Ubuff(const networkPoint *p0, const networkPoint *p1, const networkPoint *p2) { return {RK_without_Ubuff(p0), RK_without_Ubuff(p1), RK_without_Ubuff(p2)}; };
T3Tddd RK_without_Ubuff(const T_PPP &p012) { return RK_without_Ubuff(std::get<0>(p012), std::get<1>(p012), std::get<2>(p012)); };
T3Tddd RK_without_Ubuff(const networkFace *f) { return RK_without_Ubuff(f->getPoints()); };
T2Tddd RK_without_Ubuff(const networkPoint *p0, const networkPoint *p1) { return {RK_without_Ubuff(p0), RK_without_Ubuff(p1)}; };
Tddd RK_without_Ubuff_Normal(const networkFace *f) { return TriangleNormal(RK_without_Ubuff(f->getPoints())); };
Tddd RK_without_Ubuff_Normal(const networkPoint *p) {
  Tddd normal = {0., 0., 0.};
  double a = 0, total = 0;
  for (const auto &f : p->getSurfaces()) {
    a = TriangleArea(RK_without_Ubuff(f));
    FusedMultiplyIncrement(a, RK_without_Ubuff_Normal(f), normal);
    total += a;
  }
  return Normalize(normal / total);
};
Tddd getNextNormalDirichlet_BEM(const networkPoint *p) {
  Tddd normal = {0., 0., 0.};
  double a = 0, total = 0;
  for (const auto &f : p->getFacesDirichlet()) {
    a = TriangleArea(RK_without_Ubuff(f));
    // normal += a * RK_without_Ubuff_Normal(f);
    FusedMultiplyIncrement(a, RK_without_Ubuff_Normal(f), normal);
    total += a;
  }
  return Normalize(normal / total);
};
Tddd getNextNormalNeumann_BEM(const networkPoint *p) {
  Tddd normal = {0., 0., 0.};
  double a = 0, total = 0;
  for (const auto &f : p->getFacesNeumann()) {
    a = TriangleArea(RK_without_Ubuff(f));
    // normal += a * RK_without_Ubuff_Normal(f);
    FusedMultiplyIncrement(a, RK_without_Ubuff_Normal(f), normal);
    total += a;
  }
  return Normalize(normal / total);
};

// Tddd RK_with_Ubuff(const networkPoint *p) { return p->RK_X.getX(p->U_update_BEM + p->vecToSurface / p->RK_X.getdt()); };
Tddd RK_with_Ubuff(const networkPoint *p) {
  // return p->RK_X.toReachAtNextTimeQ(p->RK_X.getX(p->U_update_BEM) + p->vecToSurface);
  const auto U = p->RK_X.getVectorToReachAtNextTimeQ(RK_without_Ubuff(p) + p->vecToSurface);
  return p->RK_X.getX(U);
};
Tddd RK_with_Ubuff(const networkPoint *p, const Tddd &vecToSurface) {
  // return p->RK_X.getX(p->U_update_BEM + vecToSurface / p->RK_X.getdt());
  const auto U = p->RK_X.getVectorToReachAtNextTimeQ(RK_without_Ubuff(p) + vecToSurface);
  return p->RK_X.getX(U);
};
T3Tddd RK_with_Ubuff(const networkPoint *p0, const networkPoint *p1, const networkPoint *p2) { return {RK_with_Ubuff(p0), RK_with_Ubuff(p1), RK_with_Ubuff(p2)}; };
T3Tddd RK_with_Ubuff(const T_PPP &p012) { return {RK_with_Ubuff(p012[0]), RK_with_Ubuff(p012[1]), RK_with_Ubuff(p012[2])}; };
T3Tddd RK_with_Ubuff(const networkFace *f) { return RK_with_Ubuff(f->getPoints()); };
Tddd RK_with_Ubuff_Normal(const networkFace *f) { return TriangleNormal(RK_with_Ubuff(f->getPoints())); };
double RK_with_Ubuff_Area(const networkFace *f) { return TriangleArea(RK_with_Ubuff(f->getPoints())); };
Tddd RK_with_Ubuff_Normal(const networkPoint *p) {
  Tddd normal = {0., 0., 0.};
  double a = 0, total = 0;
  for (const auto &f : p->getSurfaces()) {
    a = RK_with_Ubuff_Area(f);
    FusedMultiplyIncrement(a, RK_with_Ubuff_Normal(f), normal);
    total += a;
  }
  return Normalize(normal / total);
};

void add_vecToSurface_BUFFER_to_vecToSurface(const auto &p, const double a = 1.) {
  if (isFinite(p->vecToSurface_BUFFER))
    p->vecToSurface += p->vecToSurface_BUFFER * a;
};

std::array<double, 3> nextPositionOnBody(Network *net, networkPoint *p) {
  auto velocity = net->velocityTranslational();
  auto rotation = net->velocityRotational();
  if (net->isRigidBody) {
    for (auto i = 0; i < 3; ++i)
      if (net->isFixed[i])
        velocity[i] = 0;
    for (auto i = 3; i < 6; ++i)
      if (net->isFixed[i])
        rotation[i - 3] = 0;
  }
  auto COM_next = net->RK_COM.getX(velocity);
  auto Q_next = Quaternion(net->RK_Q.getX(AngularVelocityTodQdt(rotation, net->RK_Q.getX())));
  // return rigidTransformation(net->ICOM, COM_next, Q_next.Rv(), p->initialX);
  return rigidTransformation(net->ICOM, COM_next, Q_next.Rv(), p->initialX);
  // else return p->X;
}

std::vector<std::tuple<T3Tddd, T3Tddd>> nextBodyVertices(const auto &Fs) {
  std::vector<std::tuple<T3Tddd, T3Tddd>> ret(Fs.size());
  int i = -1;
  for (auto &f : Fs) {
    i++;
    auto net = f->getNetwork();
    auto [p0, p1, p2] = f->getPoints();
    // if (net->isFixed)
    //    ret[i] = ToX(f);
    // else
    if (net->isRigidBody) {
      // 現在のv^nを使って問題ない．加速度はいらない．
      // Quaternion q;
      // q = q.d_dt(net->velocityRotational());
      // Quaternion next_Q(net->RK_Q.getX(net->Q.AngularVelocityTodQdt(net->velocityRotational())));
      // auto next_COM = net->RK_COM.getX(net->velocityTranslational());

      auto velocity = net->velocityTranslational();
      auto rotation = net->velocityRotational();

      if (std::ranges::all_of(velocity, [](double v) { return v == 0; }) && std::ranges::all_of(rotation, [](double v) { return v == 0; })) {
        ret[i] = {{p0->initialX, p1->initialX, p2->initialX}, {p0->initialX, p1->initialX, p2->initialX}};
        continue;
      }

      for (auto i = 0; i < 3; ++i) {
        if (net->isFixed[i])
          velocity[i] = 0;
        if (net->isFixed[i + 3])
          rotation[i] = 0;
      }

      auto next_COM = net->RK_COM.getX(velocity);
      auto next_Q = Quaternion(net->RK_Q.getX(AngularVelocityTodQdt(rotation, net->RK_Q.getX())));
      auto X_next = [&](const auto &p) { return rigidTransformation(net->ICOM, next_COM, next_Q.Rv(), p->initialX); };

      ret[i] = {{p0->X, p1->X, p2->X}, {X_next(p0), X_next(p1), X_next(p2)}};
    } else if (net->isSoftBody) {
      auto v0 = p0->velocityTranslational();
      auto v1 = p1->velocityTranslational();
      auto v2 = p2->velocityTranslational();
      if (std::ranges::all_of(net->isFixed, [](bool isfixed) { return isfixed; })) {
        v0.fill(0.);
        v1.fill(0.);
        v2.fill(0.);
      }
      auto X0 = p0->RK_X.getX(v0);
      auto X1 = p1->RK_X.getX(v1);
      auto X2 = p2->RK_X.getX(v2);
      ret[i] = {{p0->X, p1->X, p2->X}, {X0, X1, X2}};
    }
  }
  return ret;
};

// \label{BEM:vectorTangentialShift}
// Tddd vectorTangentialShift(const networkPoint *p, double scale = 1.) {
//    Tddd V = {0., 0., 0.};
//    // V += scale * ArithmeticWeightedSmoothingVector(p, [](const networkPoint *p) -> Tddd { return RK_with_Ubuff(p); });
//    // V += scale * AreaWeightedSmoothingVector(p, [](const networkPoint *p) -> Tddd { return RK_with_Ubuff(p); });
//    V += scale * DistorsionMeasureWeightedSmoothingVector2(p, [](const networkPoint *p) -> Tddd { return RK_with_Ubuff(p); });
//    return condition_Ua(V, p);
// };

// やはりこの辺りなのだろう

// \label{BEM:vectorToNextSurface}
const double dxi = 0.3;
const T3Tdd t0t1_vertices = {{{1., 0.}, {1. - dxi, dxi}, {1. - dxi, 0.}}};
const auto t0t1 = SymmetricSubdivisionOfTriangle(t0t1_vertices, 5);
const auto t0t1_high = SymmetricSubdivisionOfTriangle(t0t1_vertices, 15);

Tddd vectorToNextSurface(const networkPoint *p) {
  Tddd X_shifted = RK_with_Ubuff(p);
  Tddd X_not_shifted = RK_without_Ubuff(p);
  Tddd X, ret = {1E+20, 1E+20, 1E+20};

  struct DirectionInfo {
    double distance;
    Tddd direction;
  };
  std::vector<DirectionInfo> direction_infos;

  auto add_vector = [&](const Tddd &V, Tddd n) {
    n = Normalize(n);
    double new_dist = Dot(V, n);

    // delete if the new vector is better in a similar direction
    std::erase_if(direction_infos, [&](const auto &info) { return Dot(n, info.direction) > std::cos(M_PI / 180. * 20.) && info.distance > new_dist; });

    // Add the new vector if no existing better vector is in the same direction
    if (std::ranges::none_of(direction_infos, [&](const auto &info) { return Dot(n, info.direction) > std::cos(M_PI / 180. * 20.); })) {
      direction_infos.push_back({new_dist, n});
    }
  };

  auto faces = p->getSurfaces();

  if (p->Dirichlet || p->CORNER) {

    auto Nearest_pseudo = [&](const auto X_shifted, const networkFace *f, const auto &t0t1) {
      int index = -1;
      for (auto i = 0; i < 3; ++i)
        if (f->Points[i] == p) {
          index = i;
          break;
        }
      if (index == -1)
        throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "error");
      return Nearest(X_shifted, f->dodecaPoints[index]->interpolate(t0t1, [&](networkPoint *p) -> Tddd { return RK_without_Ubuff(p); }));
    };

    networkFace *closest_face = nullptr;

    for (const auto &f : faces) {
      if (f->Dirichlet) {
        if (_ALE_ON_LINEAR_ELEMENT_)
          X = Nearest(X_shifted, RK_without_Ubuff(f));
        else
          X = Nearest_pseudo(X_shifted, f, t0t1);

        if (Norm(ret) >= Norm(X - X_shifted)) {
          ret = X - X_shifted;
          closest_face = f;
        }
      }
    }

    if (closest_face != nullptr) {
      X = Nearest(X_shifted, RK_without_Ubuff(closest_face));
      add_vector(X - X_shifted, Normalize(X - X_shifted));
    }

    if (p->Dirichlet)
      return ret;
  }

  if (/*p->Neumann || */ p->CORNER) {

    const double short_range = 0.01; //@ short_range * radius is the detection range of the structure face

    Tddd to_structure_face = {0., 0., 0.};
    auto Vec_CurrentX012_NextX012 = nextBodyVertices(bfs(p->getContactFaces(), 1));

    //! facesはDirichet面も含んでいる．これはCORNERの張り付きにおいて重要だと思われる．これがないと，コーナーは無限に上下方向に移動しても条件を満たすことになるだろう．
    //! Dirichet面の法線方向距離も小さくする様にしている．

    auto contactAngle = [p](double distance) {
      auto max_deg = 70.;
      auto min_deg = 30.;
      auto r = distance / p->contact_range;
      auto deg = max_deg - (max_deg - min_deg) * r;
      return M_PI * deg / 180.;
    };

    /*cornerのnextBodyVerticesへの張り付きは，思わぬ歪みを生じることがある*/

    //! 面p_fが干渉する，最も近い構造物面を抽出
    const double angle = 60 * M_PI / 180.; //! 少なくとも60度は必要のようだ．変更するとしたら，0.5 * radius >= distanceの値を変更する．
    if (!Vec_CurrentX012_NextX012.empty())
      for (const auto &f : faces /*p->getFacesNeumann()*/) {
        for (const auto &CurrentX012_NextX012 : Vec_CurrentX012_NextX012) {
          auto [struct_vertex, next_struct_vertex] = CurrentX012_NextX012;
          if (isInContact(X_not_shifted /*シフトなし*/, f->normal, struct_vertex, p)) {
            /*以下はシフトした場所からのベクトルを計算*/
            auto [t0, t1, X_nearest, n] = Nearest_(X_shifted, struct_vertex);
            std::array<double, 3> To_nearest = X_nearest - X_shifted;
            double distance = Norm(To_nearest);
            if ((isFlat(n, To_nearest, contactAngle(distance)) || isFlat(n, -To_nearest, contactAngle(distance))) && (p->contact_range >= distance))
              add_vector(To_nearest, n);
            else if ((short_range * p->contact_range >= distance))
              add_vector(To_nearest, n);
          }
        }
      }

    std::vector<double> distances;
    std::vector<Tddd> directions;
    for (const auto &info : direction_infos) {
      distances.push_back(info.distance);
      directions.push_back(info.direction);
    }
    to_structure_face = optimalVector(distances, directions, {0., 0., 0.});

    return to_structure_face;
  }
  return {0., 0., 0.};
};

/*DOC_EXTRACT 0_3_1_INITIAL_VALUE_PROBLEM

### Arbitrary Lagrangian–Eulerian Methods (ALE)

水面を移動する物体が存在する場合，格子が潰れてしまうため，何らかの方法で格子を綺麗に保たなければ，長時間の計算は不可能である．
その方法の一つが，節点の位置を流速$\nabla \phi$ではなく任意のベクトル${\bf v}$で移動させる，ALEである．
${\bf v}$で移動する節点位置を${\mathbfit\chi}(t)$と置くと，
$\frac{D\phi}{Dt}=\frac{\partial\phi}{\partial t}+\nabla\phi\cdot\nabla\phi$の代わりに，
$\frac{D\phi}{Dt}=\frac{\partial\phi}{\partial t}+\frac{d{\mathbfit\chi}}{dt} \cdot\nabla\phi,\frac{d{\mathbfit\chi}}{dt} = \bf v$
を使って$\phi$を時間発展させると次時刻の節点位置${\mathbfit\chi}(t+\delta t)$での$\phi$が得られる．

ディリクレ節点（水面）：

求めた流速から，次の時刻の境界面$\Omega(t+\Delta t)$を見積もり，その面上で節点を移動させ歪さを解消する．
修正ベクトルは，$\Delta t$で割り，求めた流速$\nabla \phi$に足し合わせて，節点を時間発展させる．

ノイマン節点：

ノイマン節点も修正流速を加え時間発展させる．
ただし，ノイマン節点の修正流速に対しては，節点が水槽の角から離れないように，工夫を施している．

\ref{BEM:calculateVecToSurface}{`calculateVecToSurface`}で$\Omega(t+\Delta t)$上へのベクトルを計算する．

1. まず，\ref{BEM:vectorTangentialShift}{`vectorTangentialShift`}で接線方向にシフトし，
2. \ref{BEM:vectorToNextSurface}{`vectorToNextSurface`}で近くの$\Omega(t+\Delta t)$上へのベクトルを計算する．

#### 使っているALEの手法

`CircumradiusToInradius`をもとに，節点に隣接する面の歪みを評価し，その歪みに応じて節点を修正する．
面の歪みを修正するための節点の移動方向は，節点の向かい側にある辺の中点から辺が作る正三角形の高さの方向としている．
歪みを小さくするための方向を歪みの`CircumradiusToInradius`の微分から求めることもできるだろうが，プログラムが複雑になるためこの方法を採用している．

* 重みの最大値と最小値を設定している
* よりディリクレ面の歪みを緩和するように重みを大きくしている
* より喫水線の歪みを緩和するように重みを大きくしている

*/

// \label{BEM:calculateVecToSurface}
Tddd DistorsionMeasureWeightedSmoothingVector_modified(const networkPoint *p, const std::array<double, 3> &current_pX, std::function<Tddd(const networkPoint *)> position) {
  const int max_sum_depth = 20;
  //   auto faces = p->CORNER ? p->getFacesDirichlet() : p->getSurfaces();
  auto faces = p->getSurfaces();
  std::vector<double> weights;
  weights.reserve(faces.size());
  std::vector<Tddd> positions;
  positions.reserve(faces.size());
  Tddd X0, X1, X2, X3, Xmid, vertical, X_ideal, To_ideal;
  double W, height, normalized_discrepancy;

  // tetraがある場合
  std::vector<networkTetra *> tetras = p->getTetras();

  auto isSurface = [](const networkTetra *tet) -> bool {
    for (const auto &f : tet->Faces)
      if (f->Tetras[0] == nullptr || f->Tetras[1] == nullptr)
        return true;
    return false;
  };

  if (!tetras.empty()) {
    for (const auto &tet : tetras) {
      auto [p0, p1, p2, p3] = tet->getPoints(p);
      X0 = current_pX;
      X1 = position(p1);
      X2 = position(p2);
      X3 = position(p3);
      W = CircumradiusToInradius(X0, X1, X2, X3) - 3; // CircumradiusToInradius(X0, X1, X2)は最小で3
                                                      //   if (std::ranges::any_of(tet->Faces, [](const auto &f) { return f->SurfaceQ(); }))
                                                      //     W *= 2.;
      //! 　面それぞれに対する，修正方向は正三角形の高さの方向としている
      auto Xmid = (X1 + X2 + X3) / 3.;
      auto cross = Cross(X2 - X1, X3 - X1);
      auto vertical = Normalize(cross);
      if (Dot(vertical, X0 - Xmid) < 0)
        vertical = -vertical;
      auto area = Norm(cross) * 0.5;
      auto approx_edge_length = (Norm(X1 - X2) + Norm(X2 - X3) + Norm(X3 - X1)) / 3.;
      auto height = approx_edge_length * std::sqrt(6.) / 3.; // 正四面体の高さ
      if (Dot(vertical, current_pX - Xmid) < 0)
        vertical = -vertical;
      auto To_ideal = height * vertical + (Xmid - current_pX);
      auto normalized_discrepancy = Norm(To_ideal) / height;
      double a = 1.;
      weights.push_back(a * (normalized_discrepancy + W));
      positions.emplace_back(To_ideal);
    }
  }
  //   else
  //     throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "Tetrahedra are required for DistorsionMeasureWeightedSmoothingVector_modified.");

  {
    for (const auto &f : faces) {
      auto [p0, p1, p2] = f->getPoints(p);
      X0 = current_pX;
      X1 = position(p1);
      X2 = position(p2);
      //! 重みの計算
      W = CircumradiusToInradius(X0, X1, X2) - 2.; // CircumradiusToInradius(X0, X1, X2)は最小で2
      //! 　面それぞれに対する，修正方向は正三角形の高さの方向としている
      Xmid = (X2 + X1) * 0.5;
      vertical = Normalize(Chop(X0 - Xmid, X2 - X1));
      height = Norm(X2 - X1) * std::sqrt(3.) * 0.5;
      X_ideal = height * vertical + Xmid;
      To_ideal = X_ideal - current_pX;
      normalized_discrepancy = Norm(To_ideal) / height;
      //! 重みと位置を保存
      double a = 5.;
      //   a *= std::pow(1.1, std::ranges::count_if(f->getPoints(), [](const auto &p) { return p->isMultipleNode; })); // 重みは，面に含まれるマルチプルノードの数に応じて増加する
      weights.push_back(a * (normalized_discrepancy + W));
      positions.emplace_back(To_ideal);
    }
  }

  double sum = Sum(weights);
  if (sum < 1E-12)
    weights /= static_cast<double>(weights.size());
  else
    weights /= sum;
  Tddd V = {0., 0., 0.};
  for (int i = 0; i < weights.size(); ++i)
    V += weights[i] * positions[i];
  return V;
};

/* -------------------------------------------------------------------------- */

void calculateVecToSurface(const Network &net, const int loop, const double coefIN) {
  auto points = ToVector(net.getPoints());
  //! 初期化
  for (const auto &p : points) {
    p->temporary_bool = true;
    p->vecToSurface.fill(0.);
    p->vecToSurface_BUFFER.fill(0.);
    p->vecToSurface_BUFFER_BUFFER.fill(0.);
  }

  points = RandomSample(points);
  std::vector<networkPoint *> points_multiple_node, points_inner;
  points_multiple_node.reserve(points.size());
  points_inner.reserve(points.size());
  for (const auto &p : points) {
    if (p->isMultipleNode)
      points_multiple_node.push_back(p);
    else
      points_inner.push_back(p);
  }

  auto shiftable_distance = [&](const networkPoint *p) -> double {
    double ret = 1E+20;
    auto tet = p->getTetras();
    // もし，隣接tetraがあれば，tetraの内接円半径の平均がシフト距離の上限
    // tetraがない場合は，隣接面の内接円半径の平均をシフト距離の上限とする
    if (!tet.empty()) {
      for (const auto &t : tet) {
        auto [p0, p1, p2, p3] = t->Points;
        double dist = Inradius(RK_with_Ubuff(p0), RK_with_Ubuff(p1), RK_with_Ubuff(p2), RK_with_Ubuff(p3));
        ret = std::min(ret, dist);
      }
      return ret;
    } else {
      for (const auto &f : p->getSurfaces()) {
        auto [p0, p1, p2] = f->getPoints();
        double dist = Inradius(RK_with_Ubuff(p0), RK_with_Ubuff(p1), RK_with_Ubuff(p2));
        ret = std::min(ret, dist);
      }
      return ret;
    }
  };

  /* ------------------------------ 要素を整えるためのベクトル ----------------------------- */

  TimeWatch watch;
  double coef = coefIN;

  auto noThroughCondition = [&](const networkPoint *p, Tddd vec) {
    if (p->Neumann || p->CORNER) {
      auto Vec_CurrentX012_NextX012 = nextBodyVertices(p->getContactFaces());
      for (const auto &CurrentX012_NextX012 : Vec_CurrentX012_NextX012) {
        auto [CurrentX012, NextX012] = CurrentX012_NextX012;
        auto n_current = Normalize(Cross(CurrentX012[1] - CurrentX012[0], CurrentX012[2] - CurrentX012[0]));
        auto n_next = Normalize(Cross(NextX012[1] - NextX012[0], NextX012[2] - NextX012[0]));
        vec = Chop(vec, (n_current + n_next) / 2.);
      }
    }
    return vec;
  };

  for (auto steps = 0; steps < loop; ++steps) {

    if (steps < 3 || steps >= loop - 3)
      coef = 1E-2 * coefIN;
    else if (steps < 10)
      coef = 1E-1 * coefIN;
    else
      coef = coefIN;

    /* -------------------------------------------------------------------------- */
    for (auto &points : {points_multiple_node, points_inner}) {
#pragma omp parallel for
      for (const auto &p : points) {
        auto current_pX = RK_with_Ubuff(p);
        Tddd V = DistorsionMeasureWeightedSmoothingVector_modified(p, current_pX, [&](const networkPoint *p) -> Tddd { return RK_with_Ubuff(p); });
        if (isFinite(V)) {
          V = noThroughCondition(p, V); // 不透過条件 No through flow boundary condition
          p->vecToSurface_BUFFER = std::min(Norm(V), coef * shiftable_distance(p)) * Normalize(V);
        }
      }
      for (const auto &p : points) {
        p->vecToSurface += p->vecToSurface_BUFFER;
        p->vecToSurface_BUFFER.fill(0.);
      }
#pragma omp parallel for
      for (const auto &p : points) {
        p->clungSurface = vectorToNextSurface(p);
        // p->clungSurface = noThroughCondition(p, p->clungSurface);
        if (!isFinite(p->clungSurface))
          p->clungSurface.fill(0.);
      }
      for (const auto &p : points) {
        p->vecToSurface += std::min(Norm(p->clungSurface), shiftable_distance(p)) * Normalize(p->clungSurface);
        p->clungSurface.fill(0.);
      }
    }
    /* -------------------------------------------------------------------------- */
  }
};

// ! -------------------------------------------------------------------------- */
// !                               calculateVelocities                          */
// ! -------------------------------------------------------------------------- */

void calculateCurrentVelocities(const Network &net) {
#pragma omp parallel for
  for (const auto &p : ToVector(net.getPoints())) {
    p->U_update_BEM = p->U_BEM = gradPhi(p);
    if (p->Neumann)
      p->U_update_BEM = contactNormalVelocity(p);
  }
}

void calculateCurrentUpdateVelocities(const Network &net, const int loop, const double coef) {

  for (const auto &p : net.getPoints()) {
    if (!isFinite(p->U_update_BEM)) {
      std::cout << "p->X = " << p->X << std::endl;
      std::cout << "p->U_BEM = " << p->U_BEM << ", gradPhi(p) = " << gradPhi(p) << std::endl;
      std::cout << "p->U_update_BEM = " << p->U_update_BEM << std::endl;
      std::cout << "p->vecToSurface = " << p->vecToSurface << std::endl;
      std::cout << "p->Dirichlet = " << p->Dirichlet << std::endl;
      std::cout << "p->Neumann = " << p->Neumann << std::endl;
      std::cout << "p->CORNER = " << p->CORNER << std::endl;
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
    }
  }
  calculateVecToSurface(net, loop, coef);

  for (const auto &p : net.getPoints())
    p->U_update_BEM = p->RK_X.getVectorToReachAtNextTimeQ(p->vecToSurface + RK_without_Ubuff(p));
}

/*DOC_EXTRACT 0_7_OTHERS

### エネルギー保存則（計算精度のチェックに利用できる）

流体全体の運動エネルギーは，ラプラス方程式と発散定理を使うと，次のように境界面に沿った積分で表される．

$$
E_K =\frac{\rho}{2} \iint_\Gamma \phi\nabla\phi\cdot {\bf n} d\Gamma
$$

また，流体の位置エネルギーは，次のように表される．

$$
E_P = \frac{\rho}{2} \iint_\Gamma (0,0,g(z - z_0)^2) \cdot {\bf n} d\Gamma
$$


???   NOTE: なぜか？

      テンソルを使って考えてみると

      $$
      \begin{align*}
      \nabla \cdot (\phi\nabla\phi) &= \frac{\partial\phi}{\partial x_i} \frac{\partial\phi}{\partial x_i} + \phi \frac{\partial^2\phi}{\partial x_i \partial x_i}\\
      &= \nabla \phi \cdot \nabla \phi + \phi \nabla^2 \phi\\
      &= \nabla \phi \cdot \nabla \phi
      \end{align*}
      $$

      よって，

      $$
      \iiint_\Omega \nabla\phi\cdot\nabla\phi d\Omega = \iiint_\Omega \nabla \cdot (\phi\nabla\phi) d\Omega = \iint_\Gamma \phi\nabla\phi\cdot {\bf n} d\Gamma
      $$

      ---

      $$
      E_P = \rho g \iiint_\Omega (z - z_0) d\Omega
      = \rho g \iiint_\Omega \frac{1}{2} \nabla \cdot (0,0,(z - z_0)^2) d\Omega
      = \rho g \iint_\Gamma \frac{1}{2} (0,0,(z - z_0)^2) \cdot {\bf n} d\Gamma
      = \frac{1}{2}\rho g \iint_\Gamma (z - z_0)^2 n_z d\Gamma
      $$

      ---

*/

double KinematicEnergy(const auto &faces) {
  double EK = 0;
  for (const auto &f : faces) {
    auto [p0, p1, p2] = f->getPoints();
    EK += Dot((p0->phiphin[0] + p1->phiphin[0] + p2->phiphin[0]) / 3. * gradPhi(f), f->normal) * f->area;
  }
  return _WATER_DENSITY_ * EK / 2.;
};

double PotentialEnergy(const auto &faces) {
  double EP = 0;
  for (const auto &f : faces) {
    auto [p0, p1, p2] = f->getPoints();
    auto intpX = interpolationTriangleLinear0101(T3Tddd{p0->X, p1->X, p2->X});
    for (const auto &[x0, x1, w0w1] : __GWGW10__Tuple)
      EP += std::pow(intpX(x0, x1)[2], 2) * f->normal[2] * w0w1 * intpX.J(x0, x1);
  }
  return _WATER_DENSITY_ * EP / 2.;
};

double TotalEnergy(const auto &faces) { return KinematicEnergy(faces) + PotentialEnergy(faces); };

/*DOC_EXTRACT 0_7_OTHERS

### 内部流速の計算方法（使わなくてもいい）

[Fochesato2005](https://onlinelibrary.wiley.com/doi/10.1002/fld.838)にあるように，
流体内部の流速$\nabla \phi$は，BIEを微分して求めることができる．

$$
u({\bf a}) = \nabla\phi({\bf a}) = \int_{\partial \Omega} \frac{\partial Q}{\partial n} ({\bf x})Q({\bf x}, {\bf a}) - \phi({\bf x}) \frac{\partial Q}{\partial n} ({\bf x}, {\bf a}) d\Gamma
$$

$$
Q({\bf x},{\bf a}) = \frac{{\bf r}}{4\pi r^3}, \quad \frac{\partial Q}{\partial n} ({\bf x},{\bf a}) = \frac{1}{4\pi r^3} (3 \mathbf{n} - (\mathbf{r} \cdot \mathbf{n}) \frac{\mathbf{r}}{r^2})
$$

*/

#endif