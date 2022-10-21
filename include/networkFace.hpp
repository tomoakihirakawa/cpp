#ifndef networkFace_H
#define networkFace_H

#include "Network.hpp"
#include "NetworkUtility.hpp"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

inline Tddd networkFace::normalVelocityRigidBody(const Tddd &X) const {
   return this->normal * Dot(this->normal, this->network->velocityRigidBody(X));
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// コンストラクタ
inline networkFace::networkFace(Network *network_IN, Network *storage_IN, const T_LLL &Lines_IN, T_3P points)
    : Triangle(ToX(points)),
      network(network_IN),
      status(false),
      Lines(Lines_IN),
      Tetras({nullptr, nullptr}),
      XPoints(0) {
#ifdef DEM
   this->contactP = {};
#endif
   try {
      this->storage = storage_IN;  //なぜか初期化リストに入れれない
      this->storage->add(this);
      std::get<0>(this->Lines)->Add(this);
      std::get<1>(this->Lines)->Add(this);
      std::get<2>(this->Lines)->Add(this);
      // setBounds();
      this->Points = this->getPointsFromLines();
      T3Tddd p0p1p2_X = ToX(this->Points);
      CoordinateBounds::setBounds(p0p1p2_X);
      this->area = TriangleArea(p0p1p2_X);
      this->normal = TriangleNormal(p0p1p2_X);
      this->angles = TriangleAngles(p0p1p2_X);
   } catch (std::exception &e) {
      std::cerr << e.what() << colorOff << std::endl;
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
   };
};
// コンストラクタ
inline networkFace::networkFace(const netFp f)
    : Triangle(extractXtuple(f)),
      network(f->network),
      status(f->status),
      //   Lines(f->Lines),
      Lines(f->Lines),
      XPoints(0),
      /* ------------------------------------------------------ */
      interp_from_p0({nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}),
      interp_from_p1({nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}),
      interp_from_p2({nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}),
      /* ------------------------------------------------------ */
      //   force({0., 0., 0., 0., 0., 0.}),
      //   inertia(6, 0.),
      //   acceleration(3, 0.),
      //   velocity(6, 0.),
      //   mass(0.),
      //   center_of_mass(3, 0.),
      map_Net_ContactPoints({{nullptr, {}}}) {
   // this->grid_pull_factor = f->grid_pull_factor;

   this->Dirichlet = f->Dirichlet;
   this->Neumann = f->Neumann;
   /* ------------------------------------------------------ */
   this->xyzInverse = f->xyzInverse;
   this->normal = f->normal;
   this->angles = f->angles;
   this->area = f->area;
   // this->Lines = f->Lines;
   this->network = f->getNetwork();
   this->storage = f->getStorage();
   this->storage->add(this);
   // this->Points = f->getPoints();
   this->Points = f->getPoints();
   auto [l0, l1, l2] = f->getLinesTuple();
   // this->Lines = {l0, l1, l2};
#ifdef DEM
   this->contactP = f->contactP;
#endif
};
/* ------------------------------------------------------ */
// b% ------------------------------------------------------ */
// b% particlizeDetailsは普通のparticlizeに詳しい情報を加えて返す．
// b% 深さ毎に，面の頂点をシフトしてm線形補間に利用する．2021/11/17
using TPPP = std::tuple<networkPoint *, networkPoint *, networkPoint *>;
inline std::unordered_set<networkPoint *> networkFace::particlize(const double dx, const V_d &depth_list) {
   // depth_list: 法線方向にdx*depthだけ動かす{-1,0,1,2,3,..}など

   std::unordered_set<networkPoint *> ret;
   // double alpha;
   // T3Tddd X0X1X2;
   int count = 0;
   networkPoint *p0, *p1, *p2;
   for (const auto &d : depth_list /*double実際の長さ*/) {
      auto [p0_, p1_, p2_] = this->getPoints();
      if (count % 3 == 1) {
         p0 = p2_;
         p1 = p0_;
         p2 = p1_;
      } else if (count % 3 == 2) {
         p0 = p1_;
         p1 = p2_;
         p2 = p0_;
      } else {
         p0 = p0_;
         p1 = p1_;
         p2 = p2_;
      }
      T3Tddd X0X1X2 = {p0->getXtuple() + d / Dot(p0->getNormalTuple(), this->normal) * p0->getNormalTuple(),
                       p1->getXtuple() + d / Dot(p1->getNormalTuple(), this->normal) * p1->getNormalTuple(),
                       p2->getXtuple() + d / Dot(p2->getNormalTuple(), this->normal) * p2->getNormalTuple()};
      for (const auto &[xyz, t0t1] : triangleIntoPoints(X0X1X2, dx)) {
         auto p = new networkPoint(this->getNetwork(), xyz);
         p->particlize_info = {this, {p0, p1, p2}, t0t1, d, dx};
         ret.emplace(p);
         this->addParametricPoints(p);
      }
      count++;
   }
   return ret;
};

/*getPointsOnLines_detail
`Intersection`（Mathematica like）で面の頂点を順に取得しながら，
もし線上に`xpoint`があればxpointsを`network::sortByDistance`でソートして同時に取得いく．
結果として，this->面の隣接する線に関係する全ての点をccw周りにソートして返す．
getPointsOnLines_detail*/
/*getPointsOnLines_code*/
// V_netPp networkFace::getPointsOnLines() const
// {
// 	V_netPp ret, Ps;
// 	int s = this->Lines.size();
// 	if (s < 3)
// 		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "このFaceはのLines.size()は" + std::to_string(s) + "でfaceを形成していません");

// 	for (auto i = 0; i < s; i++)
// 	{
// 		if ((Ps = Intersection(this->Lines[i]->getPoints(), this->Lines[(i + 1) % s]->getPoints())).size() != 1)
// 			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "このFaceを構成するLineは一部繋がっていません");
// 		else
// 		{
// 			ret.emplace_back(Ps[0]);
// 			if (Lines[(i + 1) % s]->penetrateQ())
// 			{
// 				auto sortedXPoints = network::sortByDistance(Lines[(i + 1) % s]->getXPoints(), Ps[0]);
// 				ret.insert(ret.end(), sortedXPoints.begin(), sortedXPoints.end());
// 			}
// 		}
// 	}
// 	return ret;
// };
/*getPointsOnLines_code*/
/* ------------------------------------------------------ */
// VV_netPp networkFace::getPointsOnLinesDivided() const
// {
// 	auto P = getPointsOnLines();
// 	VV_i index;
// 	int s = P.size(), find1 = -1;
// 	// changing condition value when initialize process
// 	for (size_t i = 0; i < P.size(); i++)
// 	{
// 		if (P[i]->isXPoint() && find1 == -1 /*only when initialized*/)
// 			find1 = i;
// 		else if (P[i]->isXPoint() && find1 != -1)
// 			index.push_back({find1, find1 = i});
// 	}
// 	if (index.empty())
// 	{
// 		return {};
// 	}
// 	//
// 	index.push_back(V_i{(*index.rbegin())[1], index[0][0] + s /*add s for the periodicity used later*/}); // add first and last
// 	//
// 	VV_netPp ret;
// 	for (const auto &ind : index)
// 	{
// 		V_netPp tmp;
// 		for (auto i = ind[0]; i <= ind[1]; i++)
// 			tmp.emplace_back(P[i % s]);
// 		ret.emplace_back(tmp);
// 	}
// 	return ret;

// }; /*getPointsOnLinesDivided_code*/
//////////////////////////////////////////////////////////////
// #define debug_getPointsCutLines

bool isOnSameLine(const netPp p, const netPp q) {
   if (p && q)
      if (p->getXLine() == q->getXLine())
         return true;
   return false;
};
//////////////////////////////////////////////////////////////
// VV_netPp networkFace::getPointsCutLines() const
// {
// 	V_netPp x_ing = this->getPointsPenetrate();
// 	V_netPp x_ed = this->getPointsPenetrated();
// 	VV_netPp ret({});

// 	/// x_ingが2点以上ない場合：ある１辺が唯一他の面と干渉している場合，カットラインはなし，
// 	if (x_ing.size() < 2)
// 		return {};

// 	int c = 0;
// 	while (x_ing.size() >= 2)
// 	{
// #ifdef debug_getPointsCutLines
// 		std::cout << red << "x_ing " << x_ing << std::endl;
// 		std::cout << blue << "x_ed " << x_ed << colorOff << std::endl;
// #endif

// 		auto s = x_ing[0];
// 		network::erase(x_ing, s);

// 		int c_ = 0;
// 		V_netPp a_line({s});

// 		bool hit_edge = false;
// 		bool cannot_cut = false;
// 		do
// 		{
// 			auto neighbors = s->getNeighbors();
// #ifdef debug_getPointsCutLines
// 			std::cout << s << "->getNeighbors() = " << neighbors << ",  a_line=" << a_line << std::endl;
// #endif
// 			for (auto i = 0; i < neighbors.size(); i++)
// 			{
// 				// neighbors[0]<--s-->neighbors[1]

// 				//この面を貫く点を優先的に取得，アペンドしていかなければならない
// 				if (network::erase(x_ed, neighbors[i]))
// 				{
// 					a_line.emplace_back(neighbors[i]);
// 					s = neighbors[i];
// 					break;
// 				}

// 				//貫く点ではなく，この面が貫いている＝エッジ上の点と繋がっている場合，それは1つのcutlineの終点となる．
// 				else if (!isOnSameLine(s, neighbors[i]) && network::erase(x_ing, neighbors[i]))
// 				{
// 					a_line.emplace_back(neighbors[i]);
// 					hit_edge = true;
// 					break;
// 				}

// 				//エッジからエッジまで，カットラインが横断していない場合，たどり着けない場合も確かにありえる．
// 				//そのような面は，おそらく計算には用いない．ここでは，ちゅと半端なこのカットラインを無視し，残りのpenetrate点の数珠つなぎに進む．
// 				if (i == neighbors.size() - 1)
// 				{
// 					cannot_cut = true;
// 					a_line = {};
// 					break;
// 					// throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "s->getNeighbors()は含まれていない");
// 				}
// 			}
// 		} while (!hit_edge && !cannot_cut);

// 		if (!a_line.empty())
// 			ret.emplace_back(a_line);

// 		if (c++ > 100)
// 			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
// 	}
// 	return ret;
// };
// /*getPointsCutLines_code*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int chainAndDelete(V_netPp &base, VV_netPp &PC) {
   V_netPp ret;
   for (size_t i = 0; i < PC.size(); i++) {
      Chain<networkPoint> chn;
      chn.join_front(base, PC[i]);
      if (chn.isComplete /*1 or 2*/) {
         PC.erase(PC.begin() + i);  //使ったPCを消去して，使えなくする
         base = chn.chained;
         return chn.isComplete;
      }
   }
   return 0;
};
// chainAndDeleteFixingOrder ２番目の引数の反転を許さない
int chainAndDeleteFixingOrder(V_netPp &base, VV_netPp &PC) {
   V_netPp ret;
   for (size_t i = 0; i < PC.size(); i++) {
      Chain<networkPoint> chn;
      chn.join_front_fix_order(base, PC[i]);
      if (chn.isComplete /*1 or 2*/) {
         PC.erase(PC.begin() + i);
         base = chn.chained;
         return chn.isComplete;
      }
   }
   return 0;
};
   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   /**
    * これが実行できると言うことは，
    * 完全なループが全てのgetPointsOnLinesDividedで完成しなければこのループは終わらない．
    * よって，getPointsCutLinesが過不足なく接続でき，多角形を生成できたということ
    */

#define debug_getPointsCutFaces
   // VV_netPp
   // networkFace::getPointsCutFaces() const
   // {
   // 	try
   // 	{
   // 		VV_netPp PsLDiv = this->getPointsOnLinesDivided(); /*順番fixed*/
   // 		const VV_netPp PsCutL = this->getPointsCutLines();
   // 		// std::cout << "-------------------------------------------" << std::endl;
   // 		// std::cout << "PsCutL : " << PsCutL << std::endl;
   // 		// std::cout << "PsLDiv : " << PsLDiv << std::endl;
   // 		VV_netPp ret({});
   // 		// std::cout << Red << "this = " << this << colorOff << std::endl;
   // 		VV_netPp available = PsLDiv;
   // 		while (!PsLDiv.empty())
   // 		{
   // 			auto tmp = PsLDiv[0] /*順番fixed*/;
   // 			bool complete = false;
   // 			int c = 0;
   // 			auto tmp_PdCutL = PsCutL;
   // 			PsLDiv.erase(PsLDiv.begin());

   // 			do
   // 			{
   // 				// std::cin.ignore();
   // 				// std::cout << Red << "面の周囲の点群 = " << tmp << std::endl;
   // 				// std::cout << "面を横断する点群 = " << tmp_PdCutL << std::endl;
   // 				// std::cout << "面の周囲の点群 = " << PsLDiv << std::endl;
   // 				// std::cout << "        ret = " << ret << std::endl;

   // 				auto comp = chainAndDelete(tmp, tmp_PdCutL);

   // 				// std::cout << Blue << "comp = " << comp << std::endl;
   // 				// std::cout << "面の周囲の点群 = " << tmp << std::endl;
   // 				// std::cout << "面を横断する点群 = " << tmp_PdCutL << std::endl;
   // 				// std::cout << "面の周囲の点群 = " << PsLDiv << std::endl;
   // 				// std::cout << "        ret = " << ret << std::endl;

   // 				if (comp == 2)
   // 				{
   // 					tmp.pop_back();
   // 					ret.emplace_back(tmp);
   // 					break;
   // 				}
   // 				else if (comp == 1)
   // 				{
   // 					// auto dummy = PsLDiv;
   // 					int COMP = chainAndDeleteFixingOrder(tmp, PsLDiv /*自分以外では？*/ /*順番fixed*/);

   // 					// std::cout << Green << "COMP = " << COMP << std::endl;
   // 					// std::cout << "面の周囲の点群 = " << tmp << std::endl;
   // 					// std::cout << "面を横断する点群 = " << tmp_PdCutL << std::endl;
   // 					// std::cout << "面の周囲の点群 = " << PsLDiv << std::endl;
   // 					// std::cout << "        ret = " << ret << std::endl;

   // 					if (COMP == 2)
   // 					{
   // 						// PsLDiv = dummy;
   // 						tmp.pop_back();
   // 						ret.emplace_back(tmp);
   // 						break;
   // 					}
   // 					else if (COMP != 1)
   // 					{
   // 						break;

   // 						// {
   // 						//     auto tmp = this->getPointsOnLinesDivided();
   // 						//     for(auto k=0; k<tmp.size() ;k++)
   // 						//     {
   // 						//         std::cout << std::setprecision(15) << obj3D::extractX(tmp[k]) << std::endl;
   // 						//         std::cout << std::setprecision(15) << tmp[k] << std::endl;
   // 						//         mk_vtu("./vtu/getPointsOnLinesDivided" + std::to_string(k) + ".vtu", {tmp[k]});
   // 						//     }
   // 						// }

   // 						// {
   // 						// auto tmp = this->getPointsCutLines();
   // 						// for(auto k=0; k<tmp.size() ;k++)
   // 						//     {
   // 						//         std::cout << std::setprecision(15) << obj3D::extractX(tmp[k]) << std::endl;
   // 						//         std::cout << std::setprecision(15) << tmp[k] << std::endl;
   // 						//         mk_vtu("./vtu/getPointsCutLines" + std::to_string(k) + ".vtu", {tmp[k]});
   // 						//     }
   // 						// }

   // 						// std::cout << "面を横断する点群 = " << this->getPointsCutLines() << std::endl;
   // 						// std::cout << "面の周囲の点群 = " << this->getPointsOnLinesDivided() << std::endl;

   // 						// throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "接続できるものを探せなかった．これはエラー");
   // 					}
   // 				}
   // 				if (comp == 0)
   // 				{
   // 					break;

   // 					// {
   // 					//     auto tmp = this->getPointsOnLinesDivided();
   // 					//     for(auto k=0; k<tmp.size() ;k++)
   // 					//     {
   // 					//         std::cout << std::setprecision(15) << obj3D::extractX(tmp[k]) << std::endl;
   // 					//         std::cout << std::setprecision(15) << tmp[k] << std::endl;
   // 					//         mk_vtu("./vtu/getPointsOnLinesDivided" + std::to_string(k) + ".vtu", {tmp[k]});
   // 					//     }
   // 					// }

   // 					// {
   // 					// auto tmp = this->getPointsCutLines();
   // 					// for(auto k=0; k<tmp.size() ;k++)
   // 					//     {
   // 					//         std::cout << std::setprecision(15) << obj3D::extractX(tmp[k]) << std::endl;
   // 					//         std::cout << std::setprecision(15) << tmp[k] << std::endl;
   // 					//         mk_vtu("./vtu/getPointsCutLines" + std::to_string(k) + ".vtu", {tmp[k]});
   // 					//     }
   // 					// }

   // 					// std::cout << "面を横断する点群 = " << this->getPointsCutLines() << std::endl;
   // 					// std::cout << "面の周囲の点群 = " << this->getPointsOnLinesDivided() << std::endl;

   // 					// throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "接続できるものを探せなかった．これはエラー");
   // 				}
   // 			} while (c++ < 1000 /*この場合おそらく中途半端なカットラインが存在する＝エラー*/);
   // 		};
   // 		return ret;
   // 	}
   // 	catch (std::exception &e)
   // 	{
   // 		std::cerr << e.what() << colorOff << std::endl;
   // 		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
   // 	};
   // };
   /*getPointsCutFaces_code*/
   /////////////////////////////////////////////////////////////////////////////////////////////////////////
   /////////////////////////////////////////////////////////////////////////////////////////////////////////
   /////////////////////////////////////////////////////////////////////////////////////////////////////////

   // class XLoops
   // {
   // public:
   // 	VV_netPp loops;
   // 	int failure;
   // 	int success;
   // 	XLoops() : loops({}), failure(0), success(0){};
   // 	XLoops(const V_netFp f);
   // 	XLoops(const netFp f) : loops({}), failure(0), success(0)
   // 	{
   // 		VV_netPp PsLDiv = f->getPointsOnLinesDivided(); /*順番fixed*/
   // 		const VV_netPp PsCutL = f->getPointsCutLines();
   // 		VV_netPp available = PsLDiv;
   // 		while (!PsLDiv.empty())
   // 		{
   // 			auto tmp = PsLDiv[0] /*順番fixed*/;
   // 			bool complete = false;
   // 			int c = 0;
   // 			auto PC = PsCutL;
   // 			PsLDiv.erase(PsLDiv.begin());

   // 			do
   // 			{
   // 				auto comp = chainAndDelete(tmp, PC);

   // 				if (comp == 2)
   // 				{
   // 					tmp.pop_back();
   // 					this->loops.emplace_back(tmp);
   // 					success++;
   // 					break;
   // 				}
   // 				else if (comp == 1)
   // 				{
   // 					// auto dummy = PsLDiv;
   // 					int COMP = chainAndDeleteFixingOrder(tmp, PsLDiv /*自分以外では？*/ /*順番fixed*/);
   // 					if (COMP == 2)
   // 					{
   // 						// PsLDiv = dummy;
   // 						tmp.pop_back();
   // 						this->loops.emplace_back(tmp);
   // 						success++;
   // 						break;
   // 					}
   // 					else if (COMP != 1)
   // 					{
   // 						failure++;
   // 						break;
   // 					}
   // 				}
   // 				if (comp == 0)
   // 				{
   // 					failure++;
   // 					break;
   // 				}
   // 			} while (c++ < 1000 /*この場合おそらく中途半端なカットラインが存在する＝エラー*/);
   // 		};
   // 	};
   // };

   // inline XLoops::XLoops(const V_netFp fs) : loops({}), failure(0), success(0)
   // {
   // 	for (const auto &f : fs)
   // 	{
   // 		XLoops tmp(f);
   // 		this->loops.insert(this->loops.end(), tmp.loops.begin(), tmp.loops.end());
   // 		this->failure += tmp.failure;
   // 		this->success += tmp.success;
   // 	}
   // };

#endif