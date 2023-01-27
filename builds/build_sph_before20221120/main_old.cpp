int time_step;
#define simulation
#include "GNUPLOT.hpp"
#define DEM
#include <unordered_set>
#include "Network.hpp"
#include "SPH_weightingFunctions.hpp"
#include "integrationOfODE.hpp"
#include "minMaxOfFunctions.hpp"

std::string home_dir = std::getenv("HOME");
std::string output_dir = home_dir + "/SPH";

/* ------------------------------------------------------ */
/*
 * クラスの内部変数を外部関すが取り出す場合
 * 外部関数はext???と名付け，という抽出関数を通して行うことにする
 * クラスの自身が自身の内部変数を読み出す場合は，get???メソッドを使う
 * ext???はクラスのget???を使って変数を抽出することもあるだろう．
 */
/* ------------------------------------------------------ */
V_netLp getLinesAround(netPp p)
{
    V_netLp ret = {};
    for (const auto &f : p->getFaces())
        for (const auto &line : f->getLines())
            ret.emplace_back(line);
    return DeleteDuplicates(ret);
};
V_netLp getLinesAround(netLp l)
{
    V_netLp ret = {l};
    for (const auto &p : l->getPoints())
        for (const auto &f : p->getFaces())
            for (const auto &line : f->getLines())
                ret.emplace_back(line);
    return DeleteDuplicates(ret);
};
bool refine(netLp l, double len)
{
    if (l->length() > len)
    {
        l->divide();
        // Print("新しいlineが生成されることは，このループに問題を引き起こさないか？");
        auto tmp = getLinesAround(l);
        for (auto &line : tmp)
        {
            if (isFlat(line))
            {
                line->flipIfIllegal();
                auto fs = line->getFaces();
                for (auto i = 0; i < 2; i++)
                    for (auto &f : fs)
                        LaplacianSmoothingIfFlat(f->getPoints() /*内部でシャッフルする*/); //ちゃんとsetXされているかチェック
            }                                                                              // LaplacianSmoothingIfFlat(line->getPoints() /*内部でシャッフルする*/);
        }
        return true;
    }
    return false;
};
/////////////////////////////////
void remesh(const V_Netp &all_obj, const double lim_len)
{
    for (auto i = 0; i < all_obj.size(); i++)
    {
        auto obj0 = all_obj[i];
        display(obj0);
        bool found;

        //均等なメッシュを作りたいなら，limに該当する線の数を先にチェックし，全体における，
        //該当する線の割合をlimで調整しながら進めるべき．
        //長さが大きく異なる線を，適当に分割すると，見た目，構造が大きく変わってしまう．
        //構造を維持し，徐々に変化させる必要がある．急な変化はよくない
        // 該当する点の割合と調整し徐々に分割していく

        int count = 0;
        try
        {
            Print(extLength(obj0->getLines()));
            for (const auto &lim : Subdivide(2. * lim_len, lim_len, 50))
            {
                int c = 0;
                do
                {
                    found = false;
                    auto line = obj0->getLines();
                    std::shuffle(std::begin(line), std::end(line), std::default_random_engine());
                    //ループ中にlの配置が変わるので，全体をフリップできない場合がある
                    // LaplacianSmoothingIfFlat(obj0->getPoints());
                    int cc = 0;
                    for (auto j = 0; j < line.size(); j++)
                        found = refine(line[j], lim);
                    // LaplacianSmoothingIfFlat(obj0->getPoints());
                    // std::cout << Red << "c = " << c++ << reset << std::endl;
                    c++;
                } while (c <= 5 /*最低回数*/ || (found && c < 100));
                mk_vtu(output_dir + "/remesh" + std::to_string(count++) + ".vtu", obj0->getFaces());
            }
        }
        catch (std::exception &e)
        {
            std::cerr << e.what() << reset << std::endl;
            throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
        };

        // try
        // {
        //   for (const auto &lim : Reverse(Subdivide(2. * lim_len, lim_len/2., 50)))
        //   {
        //     int c = 0;
        //     do
        //     {
        //       found = false;
        //       auto line = obj0->getLines();
        //       // std::shuffle(std::begin(line), std::end(line), std::default_random_engine());
        //       std::shuffle(std::begin(line), std::end(line), std::default_random_engine());
        //       //ループ中にlの配置が変わるので，全体をフリップできない場合がある
        //       // LaplacianSmoothingIfFlat(obj0->getPoints());
        //       for (auto j = 0; j < line.size(); j++)
        //       {
        //         // std::cout << "line[j] = " << line[j] << ", lim = " << lim << std::endl;
        //         // std::cout << "line[j]->length() = " << line[j]->length() << std::endl;
        //         found = coarsen(line[j], lim);
        //         if (found)
        //           break;
        //       }
        //       // LaplacianSmoothingIfFlat(obj0->getPoints());
        //       // std::cout << Red << "c = " << c++ << reset << std::endl;
        //       c++;
        //     } while (c < 200);
        //     mk_vtu(home_dir+"/vtu/remesh" + std::to_string(count++) + ".vtu", obj0->getFaces());
        //   }
        // }
        // catch (error_message &e)
        // {
        //   std::cerr << e.what() << reset << std::endl;
        //   throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
        // }
    }
    Print("remesh done", Red);
};
////////////////////////////////////////////////////////////
bool linkIfCloser(const netPp p, const netPp q, const V_i &lowup = {5, 8})
{
    if (p == q)
    {
        //!近傍探査を続けてくださいという意味でtrue
        return true;
    }
    //!近傍の点の再リンクも同時に行う
    auto neighbors = p->getNeighbors();
    // auto neighbors = Flatten(BFS(p,1,{p->getNetwork()}));
    if (neighbors.size() < lowup[0])
    {
        // Print("link");
        link(p, q, p->getNetwork());
        p->sortLinesByLength();
        return true;
    }
    else if (!MemberQ(neighbors, q))
    {
        /*current farest dist*/
        netLp longestL = (*(p->getLines().rbegin()));
        if (longestL->length() > Norm3d(p->getX() - q->getX()))
        {
            auto replaced_P = (*longestL)(p);
            longestL->Replace(replaced_P, q);
            p->sortLinesByLength();
            q->sortLinesByLength();

            // delete if too many points
            {
                longestL = *(p->getLines().rbegin());
                if (p->getLines().size() > lowup[1] && (*longestL)(p)->getLines().size() > lowup[1])
                    delete longestL;
            }

            // delete if too many points
            {
                longestL = *(q->getLines().rbegin());
                if (q->getLines().size() > lowup[1] && (*longestL)(q)->getLines().size() > lowup[1])
                    delete longestL;
            }

            if (replaced_P->getLines().size() < lowup[0])
                for (const auto &n : neighbors)
                    linkIfCloser(replaced_P, n);

            return true;
        }
        else
        {
            //メンバーではなかったpointだったが，メンバーの中で最遠の点よりも遠かった．
            return false;
        }
    }
    else
        return true;
};

// #define linkIfCloser_findContancP_debug

bool linkIfCloser_findContancP(const netPp p, const netPp q, const V_i &lowup, const double limit_length_to_break = 1E+15)
{
    //*lowupのlowは全ての点は，最低この数以上のリンクを持つ
    //*lowupのupは，p,qの両者がこのupよりも多い数のリンクを持つことはない，という限界．ただし，片方が，100このリンクを持つことはありえる．
    //メンバーと比べて距離が遠い場合のみfalseを返す
    // int lower_lim_vectorlen = lowup[0];
    // int upper_lim_vectorlen = lowup[1];
    if (p == q)
    {
        //!近傍探査を続けてくださいという意味でtrue
        return true;
    }

    double p2q_length = Norm3d(p->getX() - q->getX());

    //!力の計算などに用いる近傍の点の場合，保持する
    //例えば，DEMの場合は，互いの粒子（半径）が重なる場合，保持する
    // SPHの場合は，影響半径内部なら保持する
    if (p2q_length - (p->radius + q->radius) < 0.)
    {
        network::add(p->contactP, q);
        network::add(q->contactP, p);
    };

    //! SPH用の近傍点を取得
    double gamma = 5.; //!各点の近傍の定義は違って良い
    if ((p2q_length - gamma * p->radius) < 0.)
        network::add(p->neighborP, q);
    if ((p2q_length - gamma * q->radius) < 0.)
        network::add(q->neighborP, p);

    //!近傍の点の再リンクも同時に行う
    auto neighbors = p->getNeighbors();
    // auto neighbors = Flatten(BFS(p,1,{p->getNetwork()}));
    if (neighbors.size() < lowup[0])
    {
        //!この点pが保有する点の数が少ない場合はかならずリンクする
#ifdef linkIfCloser_findContancP_debug
        Print("この点pが保有する点の数が少ない場合はかならずリンクする", Green);
#endif
        link(p, q, p->getNetwork());
        p->sortLinesByLength();
        return true;
    }
    else if (!MemberQ(neighbors, q))
    {
        //!この点pが保有する点の数が十分の場合，距離をチェックする
        // ここが終わらない理由は，遠くの点から握られているため．
        // その手は，振り払えない．．．しかし，回数が決まっているはずだが．．
#ifdef linkIfCloser_findContancP_debug
        Print("この点pが保有する点の数が十分の場合，距離をチェックする", Blue);
#endif
        netLp longestL = (*(p->getLines().rbegin()));
        // double p2q_length = Norm3d(p->getX() - q->getX());
        if (longestL->length() > p2q_length)
        {
            auto replaced_P = (*longestL)(p);
            longestL->Replace(replaced_P, q);
            p->sortLinesByLength();
            q->sortLinesByLength();

#ifdef linkIfCloser_findContancP_debug
            Print("相手が十分な数の線を保有している場合，この線は不要と考える");
#endif
            {
                longestL = *(p->getLines().rbegin());
                if (p->getLines().size() > lowup[1] /*十分な線を保有*/ && (*longestL)(p)->getLines().size() > lowup[1] /*十分な線を保有*/)
                    delete longestL;
            }
#ifdef linkIfCloser_findContancP_debug
            Print("相手が十分な数の線を保有している場合，この線は不要と考える");
#endif
            {
                longestL = *(q->getLines().rbegin());
                if (q->getLines().size() > lowup[1] /*十分な線を保有*/ && (*longestL)(q)->getLines().size() > lowup[1] /*十分な線を保有*/)
                    delete longestL;
            }

            if (replaced_P->getLines().size() < lowup[0])
            {
#ifdef linkIfCloser_findContancP_debug
                Print("無限ループの可能性はありえない，なぜならこの点はリンクが欠如しているので，必ずつながる事になる．");
                Print("この点pのneigborsの中から最寄りの物を選んでもらう");
#endif
                double dist = 1E+15;
                netPp closest_n = nullptr;
                for (const auto &n : neighbors)
                {
                    if (n != replaced_P)
                    {
                        auto tmp = Norm3d(n->getX() - replaced_P->getX());
                        if (dist > tmp)
                        {
                            closest_n = n;
                            dist = tmp;
                        }
                    }
                }
                link(replaced_P, closest_n, p->getNetwork());
                replaced_P->sortLinesByLength();
            }
            return true;
        }
        else
        {
            //メンバーではなかったpointだったが，メンバーの中で最遠の点よりも遠かった．
            // std::cout << red << "メンバーではなかったpointだったが，メンバーの中で最遠の点よりも遠かった．" << reset << std::endl;
            if (p2q_length > limit_length_to_break)
            {
#ifdef linkIfCloser_findContancP_debug
                Print("近傍点ではなく，ある範囲に含まれなくなった，ー＞これ以降も含まれることはない．");
#endif
                // std::cout << red << "近傍点ではなく，ある範囲に含まれなくなった，ー＞これ以降も含まれることはない．" << reset << std::endl;
                return false;
            }
            else if (p->getLines().size() >= lowup[1])
            {
#ifdef linkIfCloser_findContancP_debug
                Print("満杯状態の近傍点群に対して，近い点から比較して行ったが，既存の点よりも遠かったー＞近傍点に含まれるものはこれ以降ありえない．");
#endif
                // std::cout << blue << "満杯状態の近傍点群に対して，近い点から比較して行ったが，既存の点よりも遠かったー＞近傍点に含まれるものはこれ以降ありえない．" << reset << std::endl;
                return false;
            }
            else
                return true;
        }
    }
    else
        return true;
};

//計算を終了させる条件を，最長のものを超えた場合にしていた．
//しかし，この条件では，とても長い線を謝って保持している場合，全く終了しない．
//しかし，最低限必要な線は，最寄りの点までの線，たかだか20本(これは，引数として与えられる)なので，そこまでで強制的に終了させることにする．
//ただし，引数として与えられる線は，近い方からソートされておかなければならない
void linkIfCloser_findContancP(const netPp p, const V_netPp closer_sorted_ps, const V_i &lowup, const double limit_length_to_break = 1E+15)
{
    int max_i = closer_sorted_ps.size() < lowup[1] ? closer_sorted_ps.size() : lowup[1];
    for (auto i = 0; i < max_i; i++)
    {
        // Print(i, Red);
        if (!linkIfCloser_findContancP(p, closer_sorted_ps[i], lowup, limit_length_to_break))
            return;
    }
};

//////////////////////////////////////////////////////////////////

V_d outside_direction(const netPp p)
{
    double sum = 0., weight;
    V_d Vtothis, ret(3, 0.);
    for (const auto &q : TakeExcept(Flatten(BFS(p, 2, {p->getNetwork()})), p))
    {
        Vtothis = p->getX() - q->getX();
        weight = 1. / Norm3d(Vtothis);
        sum += weight;
        ret += weight * Vtothis;
    }
    return ret / sum;
};
//////////////////////////////////////////////////////////////////
V_d force(const netPp p)
{
    // // double radius = p->radius;
    // double k = 30.;
    // double dump = 1.;
    // V_d ret(3, 0.);
    // V_d Vtothis;
    // double del;
    // V_d X = p->getX();
    // for (const auto &q : p->getContactPoints())
    // {
    // 	Vtothis = X - q->getX();
    // 	if ((del = Norm3d(Vtothis) - (p->radius + q->radius)) < 0.)
    // 		ret += (k * pow(-del, 1.5) - dump * Dot(p->V - q->V /*dir to this*/, Vtothis)) * Vtothis;
    // }
    // return ret;

    return {0, 0, 0}; //デバッグのために
};

// void calcForceFromContancP(const V_netPp &points)
// {
// 	try
// 	{
// 		double k = 10000.;
// 		double dump = 1.;
// #ifdef _OPENMP
// 		Print("並列化@calcForceFromContancP");
// #pragma omp parallel for
// #endif
// 		for (auto i = 0; i < points.size(); i++)
// 		{ // double radius = p->radius;
// 			auto p = points[i];
// 			V_d Vtothis(3, 0.);
// 			double del;
// 			V_d X = p->getX();
// #ifdef _OPENMP
// #pragma omp critical
// #endif
// 			{
// 				V_netPp ps = p->contactP;
// 				for (const auto &q : ps)
// 				{
// 					network::erase(p->contactP, q);
// 					network::erase(q->contactP, p);
// 					Vtothis = X - q->getX();
// 					del = Norm3d(Vtothis) - (p->radius + q->radius);
// 					if (del < 0.)
// 					{
// 						V_d tmp_force = (k * pow(-del, 1.5) - dump * Dot(p->V - q->V /*dir to this*/, Vtothis)) * Vtothis;
// 						points[i]->F += tmp_force;
// 						q->F -= tmp_force;
// 					}
// 				}
// 			}
// 		}
// 	}
// 	catch (error_message &e)
// 	{
// 		e.print();
// 		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
// 	}
// };

// void calcActionAndReaction(const V_netPp &points, const V_netFp &faces)
// {
// 	try
// 	{
// 		double k = 30. * 2.;
// 		double dump = 1.;
// #ifdef _OPENMP
// 		Print("並列化@calcActionAndReaction");
// #pragma omp parallel for
// #endif
// 		for (auto i = 0; i < points.size(); i++)
// 		{
// 			auto p = points[i];
// 			V_d Vtothis;
// 			V_d X = p->getX();
// 			double del;
// 			for (const auto &f : faces)
// 			{
// 				Vtothis = Dot(X - f->getX(), f->getNormal()) * f->getNormal();
// 				// normalがX-Xfと逆の場合，normalをかけると符号がかわり，Xを向く
// 				del = p->radius - Norm3d(Vtothis);
// 				//まだヒットしたかわからない
// 				if (del > 0.)
// 				{
// 					V_d normal = f->getNormal();
// 					VV_d AB = {X + p->radius * normal, X - p->radius * normal};
// 					if (isIntersectingSurface(f->getLocations(), AB))
// 					{
// 						V_d tmp_force = Vtothis * (k * pow(del, 2.) - dump * Dot(p->V /*- q->V*/ /*dir to this*/, Vtothis));
// #ifdef _OPENMP
// #pragma omp critical
// #endif
// 						{
// 							points[i]->F += tmp_force;
// 							f->F -= tmp_force;
// 						}
// 					}
// 				}
// 			}
// 		}
// 	}
// 	catch (error_message &e)
// 	{
// 		e.print();
// 		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
// 	}
// };

// V_d force(const netPp p, const netFp f)
// {
// 	// double radius = p->radius;
// 	double k = 30.;
// 	double dump = 1.;
// 	V_d ret(3, 0.);
// 	V_d Vtothis;
// 	double del;
// 	V_d X = p->getX();
// 	for (const auto &q : p->getNeighbors())
// 	{
// 		Vtothis = Dot(X - f->getX(), f->getNormal()) * f->getNormal();
// 		// normalがX-Xfと逆の場合，normalをかけると符号がかわり，Xを向く
// 		del = p->radius - Norm3d(Vtothis);
// 		//まだヒットしたかわからない
// 		if (del > 0.)
// 		{
// 			auto normal = f->getNormal();
// 			auto AB = {X + p->radius * normal, X - p->radius * normal};
// 			if (isIntersectingSurface(f->getLocations(), AB))
// 			{
// 				ret += (k * pow(del, 1.5) - dump * Dot(p->V /*- q->V*/ /*dir to this*/, Vtothis)) * Vtothis;
// 			}
// 		}
// 	}
// 	return ret;
// };

// V_d force(const netFp f, const netPp p)
// {
// 	return -force(p, f);
// };

///////////////////////////////////////////////////

// void updatePosition_considering_refrectionDEM(const V_netPp &points, const V_netFp &faces)
// {
// 	double dt = 0.002;
// 	// Print("baloon facesとの反射を考慮したpointsの時間発展");
// 	VVV_d p0p1p2s;
// 	for (const auto &f : faces)
// 		p0p1p2s.emplace_back(f->getLocations());

// // 自分の周り以外との接触回避
// #ifdef _OPENMP
// 	Print("並列化");
// #pragma omp parallel for
// #endif
// 	for (auto i = 0; i < points.size(); i++)
// 	{
// 		auto p = points[i];
// 		auto X = p->getX();
// 		bool ishit = false;
// 		V_d dx = p->V * dt;
// 		VV_d v_Vnew = {p->V + (p->F / p->mass) * dt};
// 		VV_d v_dx = {};
// 		// int hitIndex = -1;
// 		V_i hitIndcies = {};
// 		do
// 		{
// 			// std::cout << "p = " << p << ", hitIndex = " << hitIndex << std::endl;
// 			geometry::intersectionTriangleLine LT(p0p1p2s, X, X + dx, hitIndcies);
// 			// std::cout << "p = " << p << ", LT.isIntersect = " << LT.isIntersect << std::endl;
// 			ishit = LT.isIntersect;
// 			if (ishit)
// 			{
// 				// std::cin.ignore();
// 				// hitIndex = LT.indexOfTriangle;
// 				hitIndcies.emplace_back(LT.indexOfTriangle);
// 				// std::cout << "hitIndcies = " << hitIndcies << std::endl;
// 				X = LT.X; //on wall
// 				// std::cout << "dx = " << dx << std::endl;
// 				dx = LT.vecX2B_; //残りdx
// 				// std::cout << "dx = " << dx << std::endl;
// 				v_dx.emplace_back(LT.vecA2X);
// 				auto Vnew = LT.reflectIfPossible(*v_Vnew.rbegin());
// 				v_Vnew.emplace_back(Vnew);
// 			}
// 			else
// 			{
// 				v_dx.push_back(dx);
// 			}
// 		} while (ishit);

// #ifdef _OPENMP
// #pragma omp critical
// #endif
// 		{
// 			p->setX(p->getX() + Sum(v_dx));
// 			p->V = *v_Vnew.rbegin();
// 		}
// 	}
// };

// void updatePosition_considering_refrectionSPH(const V_netPp &points, const V_netFp &faces, const double dt /*RKに指定された*/)
// {
// 	VVV_d p0p1p2s;
// 	for (const auto &f : faces)
// 		p0p1p2s.emplace_back(f->getLocations());
// // 自分の周り以外との接触回避
// #ifdef _OPENMP
// 	Print("並列化");
// #pragma omp parallel for
// #endif
// 	for (auto i = 0; i < points.size(); i++)
// 	{
// 		auto p = points[i];
// 		auto X = p->getX();
// 		bool ishit = false;
// 		V_d dx = p->U_SPH * dt; //!反射を考慮して計算してあげなければならない．
// 		/* ------------------------ dVdt ------------------------ */
// 		double nu = 0.001005 / p->density;
// 		V_d dVdt(3, 0.);
// 		/* ------------------------ 圧力勾配 ------------------------ */
// 		dVdt -= p->gradP_SPH / p->density;
// 		dVdt += nu * p->lap_U;
// 		dVdt += {0, 0, -9.81};
// 		/* ------------------------ V ------------------------ */
// 		VV_d v_Vnew = {p->U_SPH + dVdt * dt};
// 		VV_d v_dx = {};
// 		// int hitIndex = -1;
// 		V_i hitIndcies = {};
// 		do
// 		{
// 			// std::cout << "p = " << p << ", hitIndex = " << hitIndex << std::endl;
// 			geometry::intersectionTriangleLine LT(p0p1p2s, X, X + dx, hitIndcies);
// 			// std::cout << "p = " << p << ", LT.isIntersect = " << LT.isIntersect << std::endl;
// 			ishit = LT.isIntersect;
// 			if (ishit)
// 			{
// 				// std::cin.ignore();
// 				// hitIndex = LT.indexOfTriangle;
// 				hitIndcies.emplace_back(LT.indexOfTriangle);
// 				// std::cout << "hitIndcies = " << hitIndcies << std::endl;
// 				X = LT.X; //on wall
// 				// std::cout << "dx = " << dx << std::endl;
// 				double vis = 0.5;
// 				dx = LT.vecX2B_ * vis; //残りdx
// 				// std::cout << "dx = " << dx << std::endl;
// 				v_dx.emplace_back(LT.vecA2X);
// 				auto Vnew = LT.reflectIfPossible(*v_Vnew.rbegin());
// 				v_Vnew.emplace_back(Vnew);
// 			}
// 			else
// 			{
// 				v_dx.push_back(dx);
// 			}
// 		} while (ishit);

// #ifdef _OPENMP
// #pragma omp critical
// #endif
// 		{
// 			/* ------------------------ 流速の更新 ----------------------- */
// 			p->setX(p->getX() + Sum(v_dx));

// 			p->U_SPH = *v_Vnew.rbegin();
// 			/* ----------------------- 密度の更新 ----------------------- */
// 			p->density -= p->div_U / p->density * dt;
// 		}
// 	}
// };

//////////////////////////

// void updatePosition_considering_refrectionSPH_mod(const double dt, const V_netPp &points, const V_netFp &faces)
// {
// 	// double dt = 0.02 / 100.;
// 	VVV_d p0p1p2s;
// 	for (const auto &f : faces)
// 		p0p1p2s.emplace_back(f->getLocations());
// // 自分の周り以外との接触回避
// #ifdef _OPENMP
// 	Print("並列化");
// #pragma omp parallel for
// #endif
// 	for (auto i = 0; i < points.size(); i++)
// 	{
// 		auto p = points[i];
// 		auto X = p->getX();
// 		bool ishit = false;
// 		V_d dx = p->U_SPH * dt;
// 		VV_d v_Vnew = {p->U_SPH};
// 		VV_d v_dx = {};
// 		// int hitIndex = -1;
// 		V_i hitIndcies = {};
// 		do
// 		{
// 			// std::cout << "p = " << p << ", hitIndex = " << hitIndex << std::endl;
// 			geometry::intersectionTriangleLine LT(p0p1p2s, X, X + dx, hitIndcies);
// 			// std::cout << "p = " << p << ", LT.isIntersect = " << LT.isIntersect << std::endl;
// 			ishit = LT.isIntersect;
// 			if (ishit)
// 			{
// 				// std::cin.ignore();
// 				// hitIndex = LT.indexOfTriangle;
// 				hitIndcies.emplace_back(LT.indexOfTriangle);
// 				// std::cout << "hitIndcies = " << hitIndcies << std::endl;
// 				X = LT.X; //on wall
// 				// std::cout << "dx = " << dx << std::endl;
// 				double vis = 0.5;
// 				dx = LT.vecX2B_ * vis; //残りdx
// 				// std::cout << "dx = " << dx << std::endl;
// 				v_dx.emplace_back(LT.vecA2X);
// 				auto Vnew = LT.reflectIfPossible(*v_Vnew.rbegin());
// 				v_Vnew.emplace_back(Vnew);
// 			}
// 			else
// 			{
// 				v_dx.push_back(dx);
// 			}
// 		} while (ishit);

// #ifdef _OPENMP
// #pragma omp critical
// #endif
// 		{
// 			/* ------------------------ 流速の更新 ----------------------- */
// 			p->setX(p->getX() + Sum(v_dx));
// 			p->U_SPH = *v_Vnew.rbegin();
// 		}
// 	}
// };

//////////////////////////

// void save_DVdt_DrhoDt_SPH(const netPp p)
// {
// 	double nu = 0.001005 / p->density;
// 	Tddd DVdt = {0., 0., -9.81};
// 	DVdt -= p->gradP_SPH / p->density;
// 	DVdt += nu * p->lap_U;
// 	//////////////////
// 	p->DUDt = DVdt;
// 	p->DrhoDt_SPH = -p->div_U / p->density;
// };

//////////////////////////////////////////////////////////////////
double g = 9.81;
using map_P_Tddd = std::map<networkPoint *, Tddd>;
using map_P_d = std::map<networkPoint *, double>;
using VV_SorIorMap = std::vector<std::vector<std::variant<std::string,
                                                          int,
                                                          map_P_Vd,
                                                          map_P_Tddd,
                                                          map_P_d>>>;

double kernel_Bspline_3(const V_d &x, const V_d &origin, const double h)
{
    //!影響半径hは，平均的粒子間隔の定数倍にすることが多い．3とか
    double s = Norm3d(x - origin) / h;
    if (s <= 1.0)
    {
        if (s < 0.5)
            return (1. - 6. * s * s + 6. * s * s * s) / (M_PI * h * h * h);
        else
            return 2. * pow(1. - s, 3) / (M_PI * h * h * h);
    }
    else
        return 0.;
};

V_netPp extPointsInSmoothingLength(const V_netPp &ps, const netPp origin, const double h)
{
    V_netPp ret(0);
    auto X = origin->getX();
    for (const auto &p : ps)
        if (Norm3d(p->getX() - X) / h <= 2.)
            ret.emplace_back(p);
    return ret;
};

V_d extDensity(const V_netPp &points)
{
    V_d ret(points.size(), 0.);
    int i = 0;
    for (const auto &p : points)
        ret[i++] = p->density;
    return ret;
};

V_d extVolume(const V_netPp &points)
{
    V_d ret(points.size(), 0.);
    int i = 0;
    for (const auto &p : points)
        ret[i++] = p->volume;
    return ret;
};

V_d extMass(const V_netPp &points)
{
    V_d ret(points.size(), 0.);
    int i = 0;
    for (const auto &p : points)
        ret[i++] = p->mass;
    return ret;
};

V_d extPressure_SPH(const V_netPp &points)
{
    V_d ret(points.size(), 0.);
    int i = 0;
    for (const auto &p : points)
        ret[i++] = p->pressure_SPH;
    return ret;
};

V_d extStaticPressure_SPH(const V_netPp &points)
{
    V_d ret(points.size(), 0.);
    int i = 0;
    for (const auto &p : points)
        ret[i++] = p->getStaticPressure();
    return ret;
};

V_d extDrhoDt(const V_netPp &points)
{
    V_d ret(points.size(), 0.);
    int i = 0;
    for (const auto &p : points)
        ret[i++] = p->DrhoDt_SPH;
    return ret;
};

std::vector<Tddd> extU_SPH(const std::unordered_set<networkPoint *> &points)
{
    std::vector<Tddd> ret(points.size());
    int i = 0;
    for (const auto &p : points)
        ret[i++] = p->U_SPH;
    return ret;
};
V_d Norm_extU_SPH(const std::unordered_set<networkPoint *> &points)
{
    V_d ret(points.size());
    int i = 0;
    for (const auto &p : points)
        ret[i++] = Norm(p->U_SPH);
    return ret;
};
V_d Norm_extDUDt_SPH(const std::unordered_set<networkPoint *> &points)
{
    V_d ret(points.size());
    int i = 0;
    for (const auto &p : points)
        ret[i++] = Norm(p->DUDt_SPH);
    return ret;
};
double Max_Norm_extU_SPH(const std::unordered_set<networkPoint *> &points)
{
    double ret = 0, tmp;
    for (const auto &p : points)
        if (ret < (tmp = Norm(p->U_SPH)))
            ret = tmp;
    return ret;
};
double Max_Norm_extDUDt_SPH(const std::unordered_set<networkPoint *> &points)
{
    double ret = 0, tmp;
    for (const auto &p : points)
        if (ret < (tmp = Norm(p->DUDt_SPH)))
            ret = tmp;
    return ret;
};
double expected_Max_U_SPH(const std::unordered_set<networkPoint *> &points, const double dt)
{
    double ret = 0., tmp;
    for (const auto &p : points)
    {
        tmp = Norm(p->U_SPH + p->DUDt_SPH * dt);
        if (ret < tmp)
            ret = tmp;
        tmp = Norm(p->U_SPH);
        if (ret < tmp)
            ret = tmp;
    }
    return ret;
};
std::vector<Tddd> extTmpU_SPH(const V_netPp &points)
{
    std::vector<Tddd> ret(points.size());
    int i = 0;
    for (const auto &p : points)
        ret[i++] = p->tmp_U_SPH;
    return ret;
};

std::vector<Tddd> extlap_U(const V_netPp &points)
{
    std::vector<Tddd> ret(points.size());
    int i = 0;
    for (const auto &p : points)
        ret[i++] = p->lap_U;
    return ret;
};
V_d extDensityInterpolated_SPH(const V_netPp &points)
{
    V_d ret(points.size(), 0.);
    int i = 0;
    for (const auto &p : points)
        ret[i++] = p->density_interpolated_SPH;
    return ret;
};

V_d extRadius_SPH(const std::unordered_set<networkPoint *> &points)
{
    V_d ret(points.size(), 0.);
    int i = 0;
    for (const auto &p : points)
        ret[i++] = p->radius_SPH;
    return ret;
};

V_i extContactPointsSize_SPH(const std::unordered_set<networkPoint *> &points)
{
    V_i ret(points.size());
    int i = 0;
    for (const auto &p : points)
        ret[i++] = p->getContactPoints().size();
    return ret;
};

std::vector<Tddd> extMuLapGX_SPH(const V_netPp &points)
{
    std::vector<Tddd> ret(points.size());
    int i = 0;
    for (const auto &p : points)
        ret[i++] = p->mu_lap_rho_g_SPH;
    return ret;
};

std::vector<Tddd> tmp_V_SPH(const V_netPp &points)
{
    std::vector<Tddd> ret(points.size());
    // double mu = 0.001005;
    // double nu = mu / 1000.;
    V_d gravity = {0., 0., -9.81};
    int i = 0;
    for (const auto &p : points)
        ret[i++] = p->tmp_U_SPH;
    return ret;
};

V_i neighborsNumLim = {7, 7};

V_d Skewness(const VV_d &X, const V_d &a)
{
    V_d ret(3, 0.);
    V_d s(3, 0.);
    for (const auto &x : X)
    {
        for (auto i = 0; i < 3; i++)
        {
            ret += std::pow(x[i] - a[i], 3);
            s += std::pow(x[i] - a[i], 2);
        }
    }
    ret = ret / ((double)X.size());
    s = s / (double)X.size();
    return ret / s;
};

template <typename T>
int Count(const std::vector<std::vector<T *>> &vec, const std::vector<T *> &w)
{
    int ret = 0;
    if (w.empty())
    {
        for (const auto &v : vec)
            if (v.empty())
                ret++;
    }
    else
    {
        for (const auto &v : vec)
            if (v == w)
                ret++;
    }
    return ret;
};

V_netPp ext(const V_netPp &ps, const V_Netp &nets)
{
    V_netPp ret(0);
    ret.reserve(ps.size());
    for (const auto &p : ps)
        if (MemberQ(nets, p->getNetwork()))
            ret.emplace_back(p);
    return ret;
};

VV_d getWiderBounds(VV_d bounds, double scale)
{
    auto xbounds = bounds[0];
    auto ybounds = bounds[1];
    auto zbounds = bounds[2];
    return VV_d{xbounds + scale * (xbounds - Mean(xbounds)),
                ybounds + scale * (ybounds - Mean(ybounds)),
                zbounds + scale * (zbounds - Mean(zbounds))};
};

// V_netPp(const std::vector<PointsBuckets> &pBs, const V_d X, const double smoothing_len)
// {
// 	V_netPp ret = Flatten(pBs[0].getObjects(X, smoothing_len));
// 	for (auto i = 1; i < pBs.size(); i++)
// 		ret += Flatten(pBs[i].getObjects(X, smoothing_len));
// 	return ret;
// }

V_netPp getObjects(const std::vector<Buckets<networkPoint>> &pBs, const Tddd &X, const int limit_num)
{
    V_netPp ret;
    int lim_depth = 10;
    for (const auto &B : pBs)
    {
        auto tmp = Flatten(B.getObjects(X, lim_depth, limit_num));
        ret.insert(ret.end(), tmp.begin(), tmp.end());
    }
    return ret;
}

netPp getFarestPoint(const V_netPp &ps, const netPp a)
{
    netPp ret;
    double far = -0., nr;
    for (const auto &q : ps)
    {
        nr = Norm(q->getX() - a->getX());
        if (nr > far)
        {
            far = nr;
            ret = q;
        }
    }
    return ret;
}

double Median(std::vector<double> scores)
{
    size_t size = scores.size();
    if (size == 0)
    {
        return 0; // Undefined, really.
    }
    else
    {
        sort(scores.begin(), scores.end());
        if (size % 2 == 0)
        {
            return (scores[size / 2 - 1] + scores[size / 2]) / 2.;
        }
        else
        {
            return scores[size / 2];
        }
    }
}

//! ------------------------------------------------------ */
//! ------------------------- 出力 ------------------------ */
//! ------------------------------------------------------ */
uomap_P_d P_density, P_radius_SPH, P_pressure_SPH, P_div_U, P_isSurface, P_ContactMirroedPointsSize;
uomap_P_d P_ContactPointSize, P_ContactFaceSize, P_density_interpolated_SPH, P_DrhoDt_SPH, P_ContactDummyPointsSize;
uomap_P_Tddd P_grad_P, P_lap_U, P_U_SPH, P_DUDt, P_interpolated_normal_SPH, P_tmp_U_SPH, P_X, P_normal_SPH, P_cg_SPH, P_repulsive;
void output(const std::string &name, const std::unordered_set<networkPoint *> &points)
{
    auto s = points.size();
    // P_density.reserve(s);
    // P_radius_SPH.reserve(s);
    // P_pressure_SPH.reserve(s);
    // // P_div_U.reserve(s);
    // P_isSurface.reserve(s);
    // P_ContactPointsSize.reserve(s);
    // // P_density_interpolated_SPH.reserve(s);
    // P_DrhoDt_SPH.reserve(s);
    // P_ContactDummyPointsSize.reserve(s);
    // P_grad_P.reserve(s);
    // P_lap_U.reserve(s);
    // P_U_SPH.reserve(s);
    // P_DUDt.reserve(s);
    // P_interpolated_normal_SPH.reserve(s);
    // P_tmp_U_SPH.reserve(s);
    try
    {
        /* ------------------------- DEM ------------------------ */
        // P_V[p] = {p->volume};
        // P_F[p] = force(p);
        // P_outside[p] = {Norm(outside_direction(p))};
        // P_r[p] = {p->radius};
        /* ------------------------- SPH ------------------------ */
        // P_density[p] = p->density;
        // P_radius_SPH[p] = p->radius_SPH;
        // P_pressure_SPH[p] = p->pressure_SPH;
        // // P_div_U[p] = p->div_U;
        // P_grad_P[p] = p->gradP_SPH;
        // P_lap_U[p] = p->lap_U;
        // // P_mu_lap_g_rho[p] = p->mu_lap_rho_g_SPH;
        // P_DUDt[p] = p->DUDt_SPH;
        // P_DrhoDt_SPH[p] = p->DrhoDt_SPH;
        // P_U_SPH[p] = p->U_SPH;
        // // P_tmp_U_SPH[p] = p->tmp_U_SPH;
        // P_isSurface[p] = p->isSurface ? 1. : 0.;
        // P_ContactPointsSize[p] = (double)p->getContactPoints().size();
        // P_interpolated_normal_SPH[p] = p->interpolated_normal_SPH;
        // P_density_interpolated_SPH[p] = p->density_interpolated_SPH;

        // P_rho_interp[p] = {p->density_interpolated_SPH};
        // P_interpolated_normal[p] = {ToVector(p->interpolated_normal_SPH)};		}

        for (const auto &p : points)
        {
            P_density[p] = p->density;
            P_radius_SPH[p] = p->radius_SPH;
            P_isSurface[p] = p->isSurface ? 1. : 0.;
            P_X[p] = p->getXtuple();
            P_pressure_SPH[p] = p->pressure_SPH;
            P_U_SPH[p] = p->U_SPH;
            // P_div_U[p] = p->div_U;
            P_grad_P[p] = -p->gradP_SPH;
            P_lap_U[p] = p->lap_U;
            // P_mu_lap_g_rho[p] = p->mu_lap_rho_g_SPH;
            P_DUDt[p] = p->DUDt_SPH;
            P_DrhoDt_SPH[p] = p->DrhoDt_SPH;
            P_normal_SPH[p] = p->interpolated_normal_SPH;
            P_cg_SPH[p] = p->cg_neighboring_particles_SPH - p->getXtuple();
            // P_tmp_U_SPH[p] = p->tmp_U_SPH;

            auto ps = p->getContactPoints();
            P_ContactPointSize[p] = (double)ps.size();
            P_ContactFaceSize[p] = (double)(p->getContactFaces().size());

            auto INTXN = IntersectionsSphereTrianglesLines(p->getContactFaces());
            P_ContactMirroedPointsSize[p] = (double)(INTXN.get(p, p->radius_SPH).size());
            P_ContactDummyPointsSize[p] = (double)std::count_if(ps.begin(), ps.end(), [p](const auto q)
                                                                { return p->getNetwork() != q->getNetwork(); });

            P_repulsive[p] = p->repulsive_force_SPH;
        }
        std::cout << "data is made" << std::endl;
    }
    catch (std::exception &e)
    {
        std::cerr << e.what() << reset << std::endl;
        throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
    };
    /* ------------------------------------------------------ */
    try
    {
        mk_vtu(name, {points},
               {{"ρ", P_density},
                {"radius_SPH", P_radius_SPH},
                {"pressure_SPH", P_pressure_SPH},
                {"X", P_X},
                {"-∇P", P_grad_P},
                {"∇・∇U", P_lap_U},
                {"isSurface", P_isSurface},
                {"U", P_U_SPH},
                {"normal_SPH", P_normal_SPH},
                {"center of mass", P_cg_SPH},
                {"DUDt", P_DUDt},
                {"repulsive force", P_repulsive},
                {"contact face size", P_ContactFaceSize},
                {"contact point size", P_ContactPointSize},
                {"contact mirrored points size", P_ContactMirroedPointsSize},
                {"contact dummy points size", P_ContactDummyPointsSize},
                // {"density_interpolated_SPH", P_density_interpolated_SPH},
                {"DρDt_SPH", P_DrhoDt_SPH}});
    }
    catch (std::exception &e)
    {
        std::cerr << e.what() << reset << std::endl;
        throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
    };
};
void output_dummy(const std::string &name, const std::unordered_set<networkPoint *> &points)
{
    uomap_P_Tddd P_oppositeX;
    VV_VarForOutput data;
    for (const auto &p : points)
    {
        if (std::get<0>(p->particlize_info))
            P_oppositeX[p] = oppositeX(p);
    }
    data = {{"oppositeX", P_oppositeX}};
    std::cout << "data is made" << std::endl;
    mk_vtu(name, {points}, data);
};
//!十分な精度で計算できる壁面粒子の圧力を計算する
//!十分な数を取得できなかった点は，その後補完する
/* ------------------------------------------------------ */
// Fluid	   Dynamic Viscosity (Ns/m2)    Kinematic Viscosity (m2/s)
//             mu                           nu
/* ------------------------------------------------------ */
// Water	   1.00 x 10-3                  1.00 x 10-6
// Sea Water   1.07 x 10-3                  1.04 x 10-6
// Mercury     1.56 x 10-3                  1.15 x 10-7
// Kerosene    1.92 x 10-3                  2.39 x 10-4
/* ------------------------------------------------------ */
V_d zeros(3, 0.);
Tddd gravity = {0., 0., -9.81};
//注意：平滑化距離を設定する際に大事なこと
//バケットで平滑化距離範囲内の全ての点を抜き出せている
//核関数の定義は，平滑化距離内の点だけが値を持つ形になっているか．　
int order = 5;
/* ------------------------------------------------------ */
JSON settingSPH("./settingSPH.json");
// JSON water_json(std::ifstream("./water.json"));
//@ ----------------------- 粒子の初期配置条件 ---------------------- */
const double particle_spacing = stod(settingSPH["particle_spacing"])[0]; //粒子間隔        //*=20[ptcl/meter]
V_d xbounds = stod(settingSPH["xbounds"]);
V_d ybounds = stod(settingSPH["ybounds"]);
V_d zbounds = stod(settingSPH["zbounds"]);
//
Tdd buckets_xbounds = ToTdd(stod(settingSPH["buckets_xbounds"]));
Tdd buckets_ybounds = ToTdd(stod(settingSPH["buckets_ybounds"]));
Tdd buckets_zbounds = ToTdd(stod(settingSPH["buckets_zbounds"]));
geometry::CoordinateBounds bucketsBounds({buckets_xbounds, buckets_ybounds, buckets_zbounds});
//

const double C_SML_h = stod(settingSPH["C_SML_h"])[0];
const double C_SML_sigma = stod(settingSPH["C_SML_sigma"])[0];
const double C_SML = C_SML_h * C_SML_sigma;
//
const double C_CFL_velocity = stod(settingSPH["C_CFL_velocity"])[0];
const double C_CFL_accel = stod(settingSPH["C_CFL_accel"])[0];
const double max_dt = stod(settingSPH["max_dt"])[0];
//
const double mu = stod(settingSPH["mu"])[0];
//
const double C_Tait = stod(settingSPH["C_Tait"])[0];
//
const double C_artificial_viscousity_alpha = stod(settingSPH["C_artificial_viscousity_alpha"])[0];
const double C_artificial_viscousity_beta = stod(settingSPH["C_artificial_viscousity_beta"])[0];
//
const int kNS_SML = stoi(settingSPH["kNS_SML"])[0]; // dxを決めるための近傍粒子数
//
const double initial_surface_z_position = stod(settingSPH["initial_surface_z_position"])[0]; // dxを決めるための近傍粒子数
// 準備時間
const double preparation_max_dt = stod(settingSPH["preparation_max_dt"])[0];
const double preparation_time = stod(settingSPH["preparation_time"])[0];
const double preparation_time_step = stoi(settingSPH["preparation_time_step"])[0];
const double preparation_C_artificial_viscousity_alpha = stod(settingSPH["preparation_C_artificial_viscousity_alpha"])[0];
const double preparation_C_artificial_viscousity_beta = stod(settingSPH["preparation_C_artificial_viscousity_beta"])[0];
//
//@ ------------------------------------------------------------- */
PVDWriter waterPVD(output_dir + "/water.pvd");
PVDWriter oppositeXPVD(output_dir + "/oppositeX.pvd");
PVDWriter contact_faces_PVD(output_dir + "/contact_faces.pvd");
PVDWriter reflectPVD(output_dir + "/reflectX.pvd");
PVDWriter SPP_PVD(output_dir + "/SPP_X.pvd");
PVDWriter wave_makerPVD(output_dir + "/wave_maker.pvd");
//* ------------------------------------------------------ */
//*                           メイン                        */
//* ------------------------------------------------------ */
int main()
{
    for (const auto [a, b] : settingSPH())
        std::cout << std::setw(30) << a << std::setw(20) << b << std::endl;
    std::cout << Green << "Enterを押してください" << std::endl;
    std::cin.ignore();

    auto generate_network_from_file = [](const std::ifstream &ifs)
    {
        Network *tank = nullptr;
        // JSON object_JSON(std::ifstream("./tank.json"));
        JSON object_JSON(ifs);
        auto file_directory = object_JSON["objfile"][0];
        auto name = object_JSON["name"][0];
        tank = new Network(file_directory, name);
        std::cout << name << std::endl;
        if (object_JSON.find("rotate"))
        {
            std::cout << "rotate" << std::endl;
            auto theta_angle = stod(object_JSON["rotate"]);
            tank->rotate(theta_angle[0], {theta_angle[1], theta_angle[2], theta_angle[3]});
        }
        if (object_JSON.find("scale"))
        {
            std::cout << "scale" << std::endl;
            auto scale = stod(object_JSON["scale"]);
            tank->scale(scale[0], {scale[1], scale[2], scale[3]});
        }
        if (object_JSON.find("reverseNormal"))
        {
            std::cout << "reverseNormal" << std::endl;
            if (stob(object_JSON["reverseNormal"])[0])
                tank->reverseNormal();
        }
        if (object_JSON.find("translate"))
        {
            std::cout << "translate" << std::endl;
            V_d tmp = stod(object_JSON["translate"]);
            tank->translate({tmp[0], tmp[1], tmp[2]});
        }
        mk_vtu(output_dir + "/" + name + ".vtu", tank->getFaces());
        const V_d depth_list = stod(object_JSON["depth_list"]);
        double volume_of_a_particle = stod(object_JSON["volume_of_a_particle"])[0];
        double density = stod(object_JSON["density"])[0];
        for (const auto &f : tank->getFaces())
        {
            f->clearParametricPoints();
            for (const auto &p : f->particlize(particle_spacing, depth_list))
            {
                // for (const auto [xyz, t0, t1, depth, h_] : particlizeInfo(f, particle_spacing /*粒子間隔*/, depth_list))
                // auto p = new networkPoint(wall_dummy, wall_dummy, xyz);
                // p->initialX = xyz;
                p->setDensityVolume(density, volume_of_a_particle);
                p->radius_SPH = C_SML * particle_spacing;
                // p->radius_SPH = C_SML * std::pow(volume_of_a_particle, 1 / 3.);
                //以下はparticlize_infoから入手できる．
                p->normal_SPH = f->getNormalTuple();
                p->face_org = f;
                // p->particlize_info = {f, t0, t1, depth, h_};
            }
        }
        tank->setBounds();
        tank->resetInitialX();
        tank->setBounds();
        return tank;
    };

    try
    {
        /* ----------------- networkPointの生成，配置 ----------------- */
        auto net = new Network;

        auto X = Subdivide(xbounds[0], xbounds[1], (int)std::round((xbounds[1] - xbounds[0]) / particle_spacing)); //粒子のX座標
        auto Y = Subdivide(ybounds[0], ybounds[1], (int)std::round((ybounds[1] - ybounds[0]) / particle_spacing)); //粒子のY座標
        auto Z = Subdivide(zbounds[0], zbounds[1], (int)std::round((zbounds[1] - zbounds[0]) / particle_spacing)); //粒子のZ座標
        for (const auto &x : X)
            for (const auto &y : Y)
                for (const auto &z : Z)
                    new networkPoint(net, net, {x, y, z});
        //! --------------------- 質量・体積・密度の決定 -------------------- */
        int n = (X.size() - 1) * (Y.size() - 1) * (Z.size() - 1);                                        /*1少なく*/
        double Vtot = (xbounds[1] - xbounds[0]) * (ybounds[1] - ybounds[0]) * (zbounds[1] - zbounds[0]); //不変であるべき全体積
        double density = 1000.;
        double g = 9.81;
        double Mtot = density * Vtot;
        double volume = Vtot / n; //各粒子の初期体積
        for (const auto &p : net->getPoints())
        {
            p->initialX = p->getXtuple();
            p->setDensityVolume(density, volume); //質量(mass)は関数内部で自動で決められる
            p->radius_SPH = C_SML * particle_spacing;
            // p->radius_SPH = C_SML * std::pow(volume, 1 / 3.); //平滑化半径（"距離"は曖昧さがあるので使わない）
            p->pressure_SPH = density * g * (initial_surface_z_position - std::get<2>(p->getXtuple()));
        }
        {
            auto p = (*net->getPoints().begin());
            std::cout << Grid({"total volume", "total mass", "volume", "mass", "radius_SPH", "radius", "distance"}) << std::endl;
            std::cout << Grid({Vtot, Mtot, p->volume, p->mass, p->radius_SPH, p->radius, particle_spacing}) << std::endl;
        }
        //! --------------------------------------------- */
        Network *wave_maker = generate_network_from_file(std::ifstream("./wave_maker.json"));
        Network *tank_init = generate_network_from_file(std::ifstream("./tank_init.json"));
        Network *tank = generate_network_from_file(std::ifstream("./tank.json"));
        //! --------------------------------------------- */
        mk_vtu(output_dir + "/wave_maker_points.vtu", {wave_maker->getPoints()});
        mk_vtu(output_dir + "/tank_points.vtu", {tank->getPoints()});
        V_Netp rigid_bodies = {tank, wave_maker};

/*
WCSPH:
元々，圧縮性流体に対する解析手法だったSPHを，非圧縮性に適用できるように改良したものをWeakly Compressible SPH(WCSPH)と呼ぶ．
これは．Monaghan(1994)から始まったもの．
WCSPHでは，密度をTaitの式に代入してから，圧力は陽に計算する．不自然な圧力振動が生じることが知られている．

EISP:

*/
//@ WCSPH/EISPH
// #define WCSPH
#define EISPH
#define apply_polygon_boundary
        /* ------------------------------------------------------ */
        double dt = max_dt;
        int count = 0;
        std::unordered_set<networkPoint *> related_points;
        TimeWatch watch;
        double real_time = 0;
        std::unordered_set<networkPoint *> dummy_points, water_points;
        /* ------------------------------------------------------ */
        for (const auto &p : net->getPoints())
            p->mu_SPH = mu;
        /* ------------------------------------------------------ */
        //@ ------------------------------------------------------------------- */
        //@ ------------------------------------------------------------------- */
        //@                               メインループ                            */
        //@ ------------------------------------------------------------------- */
        //@ ------------------------------------------------------------------- */
        for (auto step = 0; step < 100000; step++)
        {
            if (step == 0)
            {
                wave_maker->clearBucketParametricPoints();
                wave_maker->makeBucketParametricPoints(particle_spacing / 2.5);
                wave_maker->makeBucketFaces(particle_spacing / 2.5);
                mk_vtu(output_dir + "/wave_maker_parametric_points.vtu", {wave_maker->getParametricPoints()});
            }
            if (step == 0 || step == preparation_time_step)
            {
                std::swap(tank, tank_init);
                tank->clearBucketParametricPoints();
                tank_init->clearBucketParametricPoints();
                tank->makeBucketParametricPoints(particle_spacing / 2.5);
                tank->makeBucketFaces(particle_spacing / 2.5);
                mk_vtu(output_dir + "/tank_parametric_points.vtu", {tank->getParametricPoints()});
                rigid_bodies = {tank, wave_maker};
            }
            /* ------------------------------------------------------ */
            water_points = net->getPoints();

#if defined(EISPH)
            //準備時間は，圧力を静水圧に固定する
            if (step <= 2)
                for (const auto &p : water_points)
                {
                    p->pressure_SPH = density * g * (initial_surface_z_position - std::get<2>(p->getXtuple()));
                    p->pressure_SPH_ = p->pressure_SPH;
                }
#endif
            //! いらない粒子は消して計算を始める
            std::cout << Magenta << "step :" << step << std::endl;
            std::cout << Blue << "----------- 粒子数の情報 ---------" << std::endl;
            std::cout << Grid({"fluid particles", "last dummy_points"}) << std::endl;
            std::cout << Grid({water_points.size(), dummy_points.size()}) << std::endl;
            {
                int count = 0;
                for (const auto &p : water_points)
                    if (!bucketsBounds.isInside(p->getXtuple()))
                    {
                        delete p;
                        count++;
                    }
                if (count > 0)
                    std::cout << Red << "消去粒子数: " << count << reset << std::endl;
                else
                    std::cout << Blue << "消去粒子数: " << count << reset << std::endl;
                water_points = net->getPoints();
            }
            //% ------------------------------------------------------ */
            //%                         バケットの生成                   */
            //% ------------------------------------------------------ */
            Print("バケットの生成", Green);
            std::cout << "net->getLocationsTuple() = " << net->getLocationsTuple() << std::endl;
            net->makeBucketPoints(particle_spacing / 2.);
            std::cout << Green << "Elapsed time: " << Red << watch() << reset << " s\n";
            // b% ------------------------------------------------------ */
            // b%                      近傍粒子の取得                      */
            // b% ------------------------------------------------------ */
            for (const auto &n : Join(rigid_bodies, {net}))
            {
                for (const auto &p : n->getPoints())
                {
                    p->clearContactPoints();
                    p->clearContactFaces();
                }
                for (const auto &f : n->getFaces())
                    f->clearContactPoints();
            }
            // 各流体粒子の近傍点を取得し保存
            std::vector<Tddd> oppositeX = {}, reflectX = {}, SPP_X = {};
            // b$ ------------------------------------------------------ */
            // b$                        面との接触を確認                   */
            // b$ ------------------------------------------------------ */
#ifdef _OPENMP
#pragma omp parallel
#endif
            for (const auto &p : net->getPoints())
#ifdef _OPENMP
#pragma omp single nowait
#endif
            {
                p->clearContactFaces();
                p->radius = p->radius_SPH; // Mean(extLength(p->getLines()));/
                for (const auto &n : rigid_bodies)
                    p->addContactFaces(n->getBucketFaces(), false); /**shadowあり*/
            }

            //
            for (const auto &p : net->getPoints())
            {
                auto INTXN = IntersectionsSphereTrianglesLines(p->getContactFaces());
                for (const auto &[F0, F1, X, Y, N] : INTXN.getFFXYN(p, p->radius_SPH))
                {
                    if (isFinite(X))
                        reflectX.emplace_back(X);
                    if (isFinite(Y))
                    {
                        oppositeX.emplace_back(Y);
                        if (p->isSurface)
                            SPP_X.push_back(Y + Reflect(getVectorToSPP(p), N));
                    }
                };
                if (p->isSurface)
                    SPP_X.push_back(p->getXtuple() + getVectorToSPP(p));
            }
            //% ---------------------- 平滑化距離の計算 ---------------------- */
            {
#ifdef use_variable_smoothing_length
                std::cout << Green << "可変の平滑化距離の計算" << reset << std::endl;
#ifdef _OPENMP
#pragma omp parallel
#endif
                for (const auto &p : water_points)
#ifdef _OPENMP
#pragma omp single nowait
#endif
                {
                    std::vector<networkPoint *> VP = TakeExcept(Flatten(net->getBucketPoints().getObjects(p->getXtuple(), 5 /*深さ上限*/, kNS_SML /*粒子数上限*/)), p);
                    if (!VP.empty())
                    {
                        auto mean = Mean(Distance(p, VP));

                        auto INTXN = IntersectionsSphereTrianglesLines(p->getContactFaces());
                        V_d distances;
                        for (const auto &[X, Y] : INTXN.get(p, mean))
                            distances.emplace_back(Norm(Y - p->getXtuple()));
                        for (const auto &q : VP)
                            distances.emplace_back(Norm(q->getXtuple() - p->getXtuple()));
                        std::sort(distances.begin(), distances.end());

                        if (distances.size() > 5)
                            distances = V_d(distances.begin(), distances.begin() + 5);
                        mean = Mean(distances);

                        if (mean > 1.1 * particle_spacing)
                            p->radius_SPH = p->radius_SPH / 2. + C_SML * 1.1 * particle_spacing / 2.;
                        else if (mean < 0.9 * particle_spacing)
                            p->radius_SPH = p->radius_SPH / 2. + C_SML * 0.9 * particle_spacing / 2.;
                        else
                            p->radius_SPH = p->radius_SPH / 2. + C_SML * mean / 2.;
                        p->isFreeFalling = false;
                    }
                    else
                    {
                        p->radius_SPH = p->radius_SPH / 2. + C_SML * particle_spacing / 2.;
                        p->isFreeFalling = true;
                    }
                }
#else
                std::cout << Green << "固定の平滑化距離の計算: C_SML * particle_spacing = "
                          << C_SML << " * " << particle_spacing
                          << " = " << C_SML * particle_spacing << reset << std::endl;
                for (const auto &p : water_points)
                {
                    p->radius_SPH = C_SML * particle_spacing;
                    p->isFreeFalling = false;
                }
#endif
                //% --------------- p->radius_SPHの範囲だけ点を取得 --------------- */
                std::cout << Green << "p->radius_SPHの範囲だけ点を取得" << reset << std::endl;
#ifdef _OPENMP
#pragma omp parallel
#endif
                for (const auto &p : water_points)
#ifdef _OPENMP
#pragma omp single nowait
#endif
                {
                    //@ unlimiteを使う2021/11/10
                    p->addContactPoints(net->getBucketPoints(), p->getXtuple(), p->radius_SPH * 1.2);
#ifndef apply_polygon_boundary
                    for (const auto &n : rigid_bodies)
                        p->addContactPoints(n->getBucketParametricPoints(), p->getXtuple(), p->radius_SPH * 1.2);
#endif

                    // 物体だけで法線方向は決められない．
                    // 流体と物体が接する面をまず検知し，次に
                }

                // チェック
                // std::cout << "平滑化半径の平均値：" << Mean(extRadius_SPH(net->getPoints())) << std::endl;
                // std::cout << "接触点数平均：" << Mean(extContactPointsSize_SPH(net->getPoints())) << std::endl;
            }
            std::cout << green << "Elapsed time: " << Red << watch() << reset << " s\n";
            //! 近傍粒子探査が終わったら時間ステップを決めることができる
            /* ------------------------------------------------------ */
            Print("近傍粒子探査が終わったら時間ステップを決めることができる", Green);
            //! ------------------------------------------------------ */

            auto p_V = (*std::min_element(water_points.begin(), water_points.end(),
                                          [](auto a, auto b)
                                          {
                                              auto A = a->radius_SPH / Norm(a->U_SPH);
                                              auto B = b->radius_SPH / Norm(b->U_SPH);
                                              if (!isFinite(A))
                                                  return false; //つまり，B,Aの順にする
                                              return A < B;
                                          }));

            auto min_relative_velocity = (p_V->radius_SPH / C_SML_sigma) / Norm(p_V->U_SPH);

            for (const auto &p : water_points)
            {
                Tddd velocity, n, U;
                auto INTXN = IntersectionsSphereTrianglesLines(p->getContactFaces());
                // b$ ------------------------ ポリゴン ------------------------ */
                for (const auto &[F0, F1, X, Y, N] : INTXN.getFFXYN(p, p->radius_SPH))
                {
                    if (F1)
                    {
                        velocity = F0->getNormalTuple() * Dot(Tddd{std::get<0>(F0->getNetwork()->velocity),
                                                                   std::get<1>(F0->getNetwork()->velocity),
                                                                   std::get<2>(F0->getNetwork()->velocity)},
                                                              F0->getNormalTuple()) +
                                   F1->getNormalTuple() * Dot(Tddd{std::get<0>(F1->getNetwork()->velocity),
                                                                   std::get<1>(F1->getNetwork()->velocity),
                                                                   std::get<2>(F1->getNetwork()->velocity)},
                                                              F1->getNormalTuple());
                        n = N;
                    }
                    else
                    {
                        velocity = F0->getNormalTuple() * Dot(Tddd{std::get<0>(F0->getNetwork()->velocity),
                                                                   std::get<1>(F0->getNetwork()->velocity),
                                                                   std::get<2>(F0->getNetwork()->velocity)},
                                                              F0->getNormalTuple());
                        n = N;
                    }
#if defined(free_slip_boundary_condition)
                    U = p->U_SPH - 2 * n * Dot(p->U_SPH, n) + Dot(velocity, n) * n;
#elif defined(no_slip_boundary_condition)
                    U = -p->U_SPH + Dot(velocity, n) * n;
#elif defined(zero_boundary_condition)
                    U = {0, 0, 0};
#else
                    U = -n * Dot(p->U_SPH, n) + Dot(velocity, n) * n;
#endif
#if defined(boundary_gurd)
                    U += Dot(U, n) * n * boundary_gurd_value(p, Norm(Y - p->getXtuple()), p->radius_SPH);
#endif
                    if (isFinite(U))
                        if (min_relative_velocity < (p->radius_SPH / C_SML_sigma) / Norm(U))
                            min_relative_velocity = (p->radius_SPH / C_SML_sigma) / Norm(U);
                }
                // b%$------------------------------------------------------ */
            }
            auto p_A = (*std::min_element(water_points.begin(), water_points.end(),
                                          [](auto a, auto b)
                                          {
                                              auto A = (a->radius_SPH / Norm(a->DUDt_SPH));
                                              auto B = (b->radius_SPH / Norm(b->DUDt_SPH));
                                              if (!isFinite(A))
                                                  return false;
                                              return A < B;
                                          }));
            auto min_relative_accel = std::sqrt((p_A->radius_SPH / C_SML_sigma) / Norm(p_A->DUDt_SPH));
            std::cout << Red << "----------- 時間ステップdt ---------" << std::endl;

            Tddd dt_candidates;
            if (step <= preparation_time_step)
                dt_candidates = {(step == 0 ? 1E-10 : preparation_max_dt),
                                 C_CFL_velocity * min_relative_velocity,
                                 C_CFL_accel * min_relative_accel};
            else
                dt_candidates = {(step == 0 ? 1E-10 : max_dt),
                                 C_CFL_velocity * min_relative_velocity,
                                 C_CFL_accel * min_relative_accel};

            std::cout << Red << "実際の平滑化距離hを使って評価している" << std::endl;
            std::cout << "dtの候補 = " << dt_candidates << " -> " << (dt = FiniteMin(dt_candidates)) << std::endl;
            std::cout << Red << "----------------------------------" << reset << std::endl;
            //! ------------------------------------------------------ */
            if (step % 2 == 0)
            {
                // if (!dummy_points.empty())
                // 	output_dummy(output_dir + "/dummy_points" + std::to_string(count) + ".vtu", dummy_points);
                // mk_vtu(output_dir + "/contact_dummy_point.vtu", {net->getContactPointsOfPoints(rigid_bodies)});
                /* ------------------------------------------------------ */
                std::string filename = "points" + std::to_string(count) + ".vtu";
                output(output_dir + "/" + filename, water_points);
                waterPVD.push(filename, real_time);
                waterPVD.output();
                /* ------------------------------------------------------ */
                //! メモリーリークに注意
                {
                    std::string filename = "wave_maker" + std::to_string(count) + ".vtu";
                    mk_vtu(output_dir + "/" + filename, wave_maker->getFaces());
                    wave_makerPVD.push(filename, real_time);
                    wave_makerPVD.output();
                }

                {
                    std::string filename = "contact_boundary_faces" + std::to_string(count) + ".vtu";
                    mk_vtu(output_dir + "/" + filename, net->getContactFacesOfPoints());
                    contact_faces_PVD.push(filename, real_time);
                    contact_faces_PVD.output();
                }

                {
                    std::string filename = "oppositeX" + std::to_string(count) + ".vtu";
                    mk_vtu(output_dir + "/oppositeX.vtu", {oppositeX});
                    mk_vtu(output_dir + "/" + filename, {oppositeX});
                    oppositeXPVD.push(filename, real_time);
                    oppositeXPVD.output();
                }
                {
                    std::string filename = "reflectX" + std::to_string(count) + ".vtu";
                    mk_vtu(output_dir + "/reflectX.vtu", {reflectX});
                    mk_vtu(output_dir + "/" + filename, {reflectX});
                    reflectPVD.push(filename, real_time);
                    reflectPVD.output();
                }
                {
                    std::string filename = "SPP_X" + std::to_string(count) + ".vtu";
                    mk_vtu(output_dir + "/SPP_X.vtu", {SPP_X});
                    mk_vtu(output_dir + "/" + filename, {SPP_X});
                    SPP_PVD.push(filename, real_time);
                    SPP_PVD.output();
                }

                /* ------------------------------------------------------ */
                Print("出力");
                count++;
            }
            /* ------------------------ チェック ------------------------ */
            // std::unordered_set<networkPoint *> check;
            // for (const auto &p : net->getPoints())
            // 	check.insert(p->getContactPoints().begin(), p->getContactPoints().end());
            // mk_vtu(home_dir + "/vtu/check_all_contact_points.vtu", {check});
            /* ------------------------------------------------------ */

            //@ ------------------------------------------------------ */
            //@                  ルンゲクッタを使った時間積分               */
            //@ ------------------------------------------------------ */
#ifdef WCSPH
            std::map<netPp, RungeKutta_<double> *> P_RK_rho;
#elif defined(EISPH)
            std::map<netPp, RungeKutta_<double> *> P_RK_pressure;
#endif
            std::map<netPp, RungeKutta_<Tddd> *> P_RK_U, P_RK_X;
            for (const auto &p : water_points)
            {
                P_RK_U[p] = (new RungeKutta_(dt, real_time, p->U_SPH, 4));
                P_RK_X[p] = (new RungeKutta_(dt, real_time, p->getXtuple(), 4));
#ifdef WCSPH
                P_RK_rho[p] = (new RungeKutta_(dt, real_time, p->density, 4));
#elif defined(EISPH)
                P_RK_pressure[p] = (new RungeKutta_(dt, real_time, p->pressure_SPH, 4));
#endif
            }
            /* ------------------------------------------------------ */
            for (const auto &n : Join(rigid_bodies, {net}))
                for (const auto &p : n->getPoints())
                {
                    p->lap_U = {0, 0, 0};
                    p->gradP_SPH = {0, 0, 0};
                    p->isFreeFalling = false;
                }
            do
            {
                dt = P_RK_X[*water_points.begin()]->getdt();
                std::cout << "dt = " << dt << std::endl;
                /* ------------------------------------------------------ */
                double start_move = 0.1;
                if (real_time >= start_move)
                {
                    double t = real_time - start_move;
                    double h = 0.1;
                    double L = 0.25;
                    double w = std::sqrt(M_PI * 9.8 / L * tanh(M_PI * h / L));
                    double A = 0.005;
                    wave_maker->acceleration = {-A * w * w * cos(w * t), 0, 0, 0, 0, 0};
                    wave_maker->velocity = {-A * w * sin(w * t), 0, 0, 0, 0, 0};
                    wave_maker->translateFromInitialX({A * (cos(w * t) - 1), 0, 0});
                }
                /* ------------------------------------------------------ */
                //* ------------------------------------------------------------ */
                // substep          壁面粒子またはダミー粒子の流速を設定               */
                //* ------------------------------------------------------------ */
#ifndef apply_polygon_boundary
                Print("接触するダミー粒子を抜き出す", Green);
                dummy_points.clear();
                dummy_points = net->getContactPointsOfPoints(rigid_bodies);
                std::cout << "ダミー粒子の数：" << dummy_points.size() << std::endl;
#else
                std::cout << "ポリゴン壁を使う．ダミー粒子は使わない" << std::endl;
#endif
                //@ ------------------------------------------------------ */
                //@                 ダミー粒子を更新する関数                   */
                //@ ------------------------------------------------------ */
                auto updateDummyPressure = [&net, dummy_points, &rigid_bodies]()
                {
#ifndef apply_polygon_boundary
                    Print("ダミー粒子（鏡映関係を使う）の圧力を計算．WCSPHなら圧力の評価の\"後\"．EISPHなら圧力の評価の\"前\"", Green);
                    auto faces = net->getContactFacesOfPoints();
#ifdef _OPENMP
#pragma omp parallel
#endif
                    for (const auto &p : dummy_points)
#ifdef _OPENMP
#pragma omp single nowait
#endif
                    {
                        //ダミー粒子は，反対側にある粒子と同じ物性値（圧力，質量，密度）を持つ．（see Nomeritae(2016)）
                        // p->pressure_SPH = dummy_pressure_Asai(net->getBucketPoints(), p, 2);
                        p->pressure_SPH = dummy_pressure_Asai_Modified(net->getBucketPoints(),
                                                                       p,
                                                                       faces,
                                                                       4. * 0.08,
                                                                       2);
                    }
#endif
                };
                auto updateDummy_U_tmpU_mass_rho = [&net, &dummy_points, &rigid_bodies]()
                {
#ifndef apply_polygon_boundary
                    auto faces = net->getContactFacesOfPoints();
#ifdef _OPENMP
#pragma omp parallel
#endif
                    for (const auto &cp : dummy_points)
#ifdef _OPENMP
#pragma omp single nowait
#endif
                    {
                        auto opX = getClosestFacePoint(cp, faces, 4. * 0.08);
                        // auto p_at_Y = B_water.getObjects_unorderedset(Y, 8 /*深さ最大*/, 1 /*自身を含めるので２点*/);
                        auto ps = net->getBucketPoints().getObjects_unorderedset(opX, 8 /*深さ最大*/, 1 /*これだけとったら終わり*/);
                        if (!ps.empty())
                        {
                            auto [U, tmp_U, mass, density] = U_tmpU_mass_density_SPH_IDW(ps, opX, 2);
                            //流速をゼロとする場合
                            cp->U_SPH = {0, 0, 0};
                            cp->tmp_U_SPH = {0, 0, 0};
                            //フリースリップ条件
                            // cp->U_SPH = cp->reflect(U);
                            // cp->tmp_U_SPH = cp->reflect(tmp_U);
                            //ノースリップ条件
                            // cp->U_SPH = -U;
                            // cp->tmp_U_SPH = -tmp_U;
                            cp->setDensity_ConstantVolume(density);
                        }
                        else
                        {
                            cp->U_SPH = {0, 0, 0};
                            cp->tmp_U_SPH = {0, 0, 0};
                            cp->setDensity_ConstantVolume(1000.);
                            // cp->mass = 1E-10;
                            // cp->volume = 1E-10;
                            // cp->density = 1E-10;
                            // cp->setDensity(density);
                        }
                    }
#endif
                };

                //! ------------------------------------------------------ */
                //! ------------------------------------------------------ */

                Print("ラプラシアンUの計算の前にダミー粒子の流速を計算", Green);
                // std::cout << "tank->getParametricPoints().size() = " << tank_param_points.size() << std::endl;
                updateDummy_U_tmpU_mass_rho(); //! ダミー

                std::cout << green << "Elapsed time: " << Red << watch() << reset << " s\n";
                // b$ ------------------------------------------------------ */
#ifdef apply_polygon_boundary
                // b$ ------------------------------------------------------ */
                // b$                        面との接触を確認                   */
                // b$ ------------------------------------------------------ */
#ifdef _OPENMP
#pragma omp parallel
#endif
                for (const auto &p : water_points)
#ifdef _OPENMP
#pragma omp single nowait
#endif
                {
                    p->clearContactFaces();
                    p->radius = p->radius_SPH; // Mean(extLength(p->getLines()));/
                    for (const auto &n : rigid_bodies)
                        p->addContactFaces(n->getBucketFaces(), false); /**shadowあり*/
                }
#endif
                // b$ ------------------------------------------------------ */
                //* ------------------------------------------------------ */
                //*                  流体粒子の法線方向を計算                  */
                //*  A. Krimi, M. Jandaghian, and A. Shakibaeinia, Water (Switzerland), vol. 12, no. 11, pp. 1–37, 2020.
                //* ------------------------------------------------------ */
                // 法線方向の計算の理解
                // 参考文献を描いてください．
                // 質量の使い方に注意
                // 表示させてチェック
                std::cout << green << "Elapsed time: " << Red << watch() << reset << " s\n";
                Print("法線方向を計算", Green);
                // auto tank_param_points = tank->getParametricPoints();
#ifdef _OPENMP
#pragma omp parallel
#endif
                for (const auto &p : water_points)
#ifdef _OPENMP
#pragma omp single nowait
#endif
                {
                    std::vector<Tddd> Xs;
                    for (const auto &qq : p->getContactPoints())
                    {
                        auto INTXN = IntersectionsSphereTrianglesLines(qq->getContactFaces());
                        for (const auto &[F0, F1, X, Y, N] : INTXN.getFFXYN(qq, qq->radius_SPH))
                            if (Norm(Y - p->getXtuple()) <= p->radius_SPH)
                                Xs.emplace_back(Y);
                        if (Norm(qq->getXtuple() - p->getXtuple()) <= p->radius_SPH)
                            Xs.emplace_back(qq->getXtuple());
                    }
                    p->interpolated_normal_SPH = -normal(Xs, p, p->radius_SPH);
                    p->cg_neighboring_particles_SPH = cg_neighboring_particles(Xs, p) + p->getXtuple();
                    if (!isFinite(p->interpolated_normal_SPH))
                        p->interpolated_normal_SPH = {1, 0, 0};
#ifndef apply_polygon_boundary
                    p->isSurface = std::none_of(p->getContactPoints().begin(), p->getContactPoints().end(),
                                                [p](const auto &q)
                                                { return (q != p &&
                                                          (Norm(q->getXtuple() - p->getXtuple()) < p->radius_SPH / 1.1) &&
                                                          (MyVectorAngle(p->interpolated_normal_SPH, Normalize(q->getXtuple() - p->getXtuple())) < M_PI / 4.)); });
#else
                    p->isSurface = std::none_of(Xs.begin(), Xs.end(),
                                                [p](const auto &X)
                                                {
                                                    return ((Norm(X - p->getXtuple()) > 1E-10) &&
                                                            (Norm(X - p->getXtuple()) < p->radius_SPH / 1.1) &&
                                                            (MyVectorAngle(p->interpolated_normal_SPH, Normalize(X - p->getXtuple())) < M_PI / 4.));
                                                });
#endif
                }
                std::cout << green << "Elapsed time: " << Red << watch() << reset << " s\n";
#ifdef WCSPH
                // b# ------------------------------------------------------ */
                // b#                          WCSPH                         */
                // b# ------------------------------------------------------ */
                /*
                 * 1. DρDtを後で計算するために∇・Uを計算
                 * 1.2 ∇・Uを計算
                 * 1.3 密度DρDtは，DρDt=-ρ∇・Uにより計算し，密度ρは時間発展させる <-- EISPHと違うところ
                 */
                Print("流体粒子のラプラシアンをまず計算", Green);
                Print("WCSPH: 1. DUDtを後で計算するために∇Pを計算", Blue);
                Print("WCSPH: 1.1 流体粒子（テイトの式を使う）とダミー粒子（鏡映関係を使う）の圧力を計算 & DrhoDtを計算", Blue);
#ifdef _OPENMP
#pragma omp parallel
#endif
                for (const auto &p : water_points)
#ifdef _OPENMP
#pragma omp single nowait
#endif
                {
#ifdef apply_polygon_boundary
                    p->lap_U = laplacian_U_Monaghan1992_polygon_boundary(p->getContactPoints(net), p, p->radius_SPH, C_SML_sigma, dt,
                                                                         (step <= preparation_time_step) ? preparation_C_artificial_viscousity_alpha : C_artificial_viscousity_alpha,
                                                                         (step <= preparation_time_step) ? preparation_C_artificial_viscousity_beta : C_artificial_viscousity_beta);
#else
                    p->lap_U = laplacian_U_Monaghan1992(p->getContactPoints(), p, p->radius_SPH, C_SML_sigma,
                                                        (step <= preparation_time_step) ? preparation_C_artificial_viscousity_alpha : C_artificial_viscousity_alpha,
                                                        (step <= preparation_time_step) ? preparation_C_artificial_viscousity_beta : C_artificial_viscousity_beta);
#endif
                    if (!isFinite(p->lap_U))
                    {
                        std::stringstream ss;
                        ss << "p->getContactPoints(net),size()=" << p->getContactPoints(net).size() << std::endl
                           << "cg_neighboring_particles_SPH=" << p->cg_neighboring_particles_SPH << std::endl
                           << "p->interpolated_normal_SPH = " << p->interpolated_normal_SPH << std::endl
                           << "p->density = " << p->density << std::endl
                           << "p->U_SPH = " << p->U_SPH << std::endl
                           << "p->getXtuple() = " << p->getXtuple() << std::endl
                           << "p->isSurface = " << p->isSurface << std::endl
                           << "p->lap_U = " << p->lap_U << std::endl;
                        throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
                    }

#ifdef apply_polygon_boundary
                    p->div_U = div_U_polygon_boundary(p->getContactPoints(net), p, p->radius_SPH, dt);
#else
                    p->div_U = div_U(p->getContactPoints(), p, p->radius_SPH);
#endif
                    p->DrhoDt_SPH = -p->density * p->div_U;
                    p->pressure_SPH = p->pressure_Tait(p->density, C_Tait); //仮の密度で計算してみる．
                    p->pressure_SPH_ = p->pressure_SPH;
                }
                // b# ------------------------------------------------------ */
                // b#                         圧力勾配                        */
                // b# ------------------------------------------------------ */
                Print("流体粒子の圧力勾配∇Pを計算", Green);
#ifdef _OPENMP
#pragma omp parallel
#endif
                for (const auto &p : water_points)
#ifdef _OPENMP
#pragma omp single nowait
#endif
                {
#ifdef apply_polygon_boundary
                    p->gradP_SPH = grad_P_Monaghan1992_polygon_boundary_(p->getContactPoints(net), p, p->radius_SPH);
#else
                    p->gradP_SPH = grad_P_Monaghan1992(p->getContactPoints(), p, p->radius_SPH);
#endif
                    if (!isFinite(p->gradP_SPH))
                    {
                        std::stringstream ss;
                        ss << "p->getContactPoints(net),size()=" << p->getContactPoints(net).size() << std::endl
                           << "cg_neighboring_particles_SPH=" << p->cg_neighboring_particles_SPH << std::endl
                           << "p->interpolated_normal_SPH = " << p->interpolated_normal_SPH << std::endl
                           << "p->density = " << p->density << std::endl
                           << "p->U_SPH = " << p->U_SPH << std::endl
                           << "p->getXtuple() = " << p->getXtuple() << std::endl
                           << "p->isSurface = " << p->isSurface << std::endl
                           << "p->lap_U = " << p->lap_U << std::endl
                           << "p->gradP_SPH = " << p->gradP_SPH << std::endl;
                        throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
                    }
                }
                std::cout << Green << "Elapsed time: " << Red << watch() << reset << " s\n";
                // b# ------------------------------------------------------ */
                // b# ------------------------------------------------------ */
                // b# ------------------------------------------------------ */

#elif defined(EISPH)
                //! ====================================================================================== */
                //! フラクショナルステップ
                //! ====================================================================================== */
                //@ -------------------------------------------------------------- */
                // step1          　　_　  流体粒子のラプラシアンをまず計算   　  　_      _      */
                //@ -------------------------------------------------------------- */
                Print("流体粒子のラプラシアンをまず計算", Green);
#ifdef _OPENMP
#pragma omp parallel
#endif
                for (const auto &p : water_points)
#ifdef _OPENMP
#pragma omp single nowait
#endif
                {
                    {
#ifdef apply_polygon_boundary
                        p->lap_U = laplacian_U_ShaoAndLo2003_polygon_boundary(p->getContactPoints(net),
                                                                              p,
                                                                              p->radius_SPH,
                                                                              dt);
#else
                        p->lap_U = laplacian_U_ShaoAndLo2003(p->getContactPoints(),
                                                             p,
                                                             p->radius_SPH);
#endif
                    }
                    if (!isFinite(p->lap_U))
                    {
                        std::stringstream ss;
                        ss << "p->getContactPoints(net),size()=" << p->getContactPoints(net).size() << std::endl
                           << "cg_neighboring_particles_SPH=" << p->cg_neighboring_particles_SPH << std::endl
                           << "p->interpolated_normal_SPH = " << p->interpolated_normal_SPH << std::endl
                           << "p->density = " << p->density << std::endl
                           << "p->U_SPH = " << p->U_SPH << std::endl
                           << "p->getXtuple() = " << p->getXtuple() << std::endl
                           << "p->isSurface = " << p->isSurface << std::endl
                           << "p->lap_U = " << p->lap_U << std::endl;
                        throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
                    }
                    p->tmp_U_SPH = p->U_SPH + ((/*動粘性係数*/ p->mu_SPH / p->density) * p->lap_U + gravity) * dt;
                }
                std::cout << green << "Elapsed time: " << Red << watch() << reset << " s\n";
                //@ ------------------------------------------------------- */
                //  step2                   仮想的な更新                       */
                //@ ------------------------------------------------------- */
                Print("仮想的な更新", Green);
                for (const auto &p : water_points)
                    p->setX(p->getXtuple() + p->tmp_U_SPH * dt);

                Print("ダミー粒子の流速を計算U_SPH, tmp_U_SPH！！！！", Magenta);
                updateDummy_U_tmpU_mass_rho(); //! ダミー

                //! ------------------------------------------------------ */
                //!                   　EISPH ∇Pの計算                      */
                //! ------------------------------------------------------ */
#ifdef _OPENMP
#pragma omp parallel
#endif
                for (const auto &p : water_points)
#ifdef _OPENMP
#pragma omp single nowait
#endif
                {
#ifdef apply_polygon_boundary
                    p->div_U = div_tmp_U_polygon_boundary(p->getContactPoints(net), p, p->radius_SPH, dt);
#else
                    p->div_U = div_tmp_U(p->getContactPoints(), p, p->radius_SPH);
#endif
                    // 発散ありの場なので，密度は変化している．これの状態で求めた発散より密度を変化させて，
                    // 圧力の計算に用いれば，単純な方法よりもよりソレノイダルに戻す方向に働くはず．
                    // p->tmp_density = p->density + (-p->density * p->div_U) * dt;
                    p->setDensity(p->density + (-p->density * p->div_U) * dt);
                }
                /* ------------------------------------------------------ */
                Print("ダミー粒子の流速を計算U_SPH, tmp_U_SPH！！！！", Magenta);
                updateDummy_U_tmpU_mass_rho(); //! ダミー
                updateDummyPressure();         //! ダミー
                std::cout << Green << "Elapsed time: " << Red << watch() << reset << " s\n";
                /*
                    U(n+1) = (U*) - 1/(rho*) * ∇P(n+1) * dt　とするのが，EISPH
                    この∇P(n+1)はSPHの補間方法で計算する．
                    そのため，まずはP(n+1)の計算が必要．

                    * 次のステップの圧力P(n+1)を非圧縮性条件をうまく考慮した方法で計算
                    * 流体粒子の圧力勾配∇Pを計算
                    * DUDtを計算
                    * 密度は初期の密度に戻す
                */
                Print("EISPH: 1. DUDtを後で計算するために∇Pを計算", Red);
                Print("EISPH: 1.1 ダミー粒子（鏡映関係を使う）の圧力を計算済みの状態で，流体粒子（pressure_EISPH_Hosseini2007）", Red);

#ifdef _OPENMP
#pragma omp parallel
#endif
                for (const auto &p : water_points)
#ifdef _OPENMP
#pragma omp single nowait
#endif
                {
                    // EISPHの圧力計算は初期の圧力を使うため，ルンゲクッタステップ毎に圧力を初期値に戻す必要がある．
                    p->pressure_SPH = P_RK_pressure[p]->getXinit();
                }

                /* ------------------------------------------------------ */

#ifdef _OPENMP
#pragma omp parallel
#endif
                for (const auto &p : water_points)
#ifdef _OPENMP
#pragma omp single nowait
#endif
                {
                    //@ ------------------ 圧力の計算にルンゲクッタを使わない場合 ----------------- */
                    {
#ifdef apply_polygon_boundary
                        p->pressure_SPH_ = pressure_EISPH_Hosseini2007_polygon_boundary(p->getContactPoints(net), p, p->radius_SPH, dt);
#else
                        p->pressure_SPH_ = pressure_EISPH_Hosseini2007(p->getContactPoints(), p, p->radius_SPH, dt);
#endif
                    }

                    if (!isFinite(p->pressure_SPH_))
                    {
                        std::stringstream ss;
                        ss << "p->getContactPoints(net),size()=" << p->getContactPoints(net).size() << std::endl
                           << "cg_neighboring_particles_SPH=" << p->cg_neighboring_particles_SPH << std::endl
                           << "p->interpolated_normal_SPH = " << p->interpolated_normal_SPH << std::endl
                           << "p->density = " << p->density << std::endl
                           << "p->U_SPH = " << p->U_SPH << std::endl
                           << "p->getXtuple() = " << p->getXtuple() << std::endl
                           << "p->isSurface = " << p->isSurface << std::endl
                           << "p->pressure_SPH_ = " << p->pressure_SPH_ << std::endl
                           << "p->pressure_SPH = " << p->pressure_SPH << std::endl
                           << "p->lap_U = " << p->lap_U << std::endl
                           << "p->div_U = " << p->div_U << std::endl
                           << "p->gradP_SPH = " << p->gradP_SPH << std::endl;
                        throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
                    }
                    //@ ------------------------------------------------------------------- */
                }

                // b% ------------------ 圧力の計算にルンゲクッタを使う場合 ----------------- */
#ifdef _OPENMP
#pragma omp parallel
#endif
                for (const auto &p : water_points)
#ifdef _OPENMP
#pragma omp single nowait
#endif
                {
                    auto rk = P_RK_pressure[p];
                    rk->push(p->DPDt_SPH);
                    p->pressure_SPH = rk->getX();
                }
                // % ---------------------------------------------------------------- */

                // DPDtを計算できればRKを支える．
                // 	Pを計算した後に，密度は仮の密度から一定の値に戻す必要がある．
                //! ------------------------------------------------------ */
                //! ------------------------------------------------------ */

                updateDummyPressure(); //! ダミー

                //@ ------------------------------------------------------ */
                //@                         圧力勾配                        */
                //@ ------------------------------------------------------ */
                Print("流体粒子の圧力勾配∇Pを計算", Green);

#ifdef _OPENMP
#pragma omp parallel
#endif
                for (const auto &p : water_points)
#ifdef _OPENMP
#pragma omp single nowait
#endif
                {
#ifdef apply_polygon_boundary
                    p->gradP_SPH = grad_P_Monaghan1992_polygon_boundary_(p->getContactPoints(net), p, p->radius_SPH);
#else
                    p->gradP_SPH = grad_P_Monaghan1992_(p->getContactPoints(), p, p->radius_SPH);
#endif
                    if (!isFinite(p->gradP_SPH))
                    {
                        std::stringstream ss;
                        ss << "p->getContactPoints(net),size()=" << p->getContactPoints(net).size() << std::endl
                           << "cg_neighboring_particles_SPH=" << p->cg_neighboring_particles_SPH << std::endl
                           << "p->interpolated_normal_SPH = " << p->interpolated_normal_SPH << std::endl
                           << "p->density = " << p->density << std::endl
                           << "p->U_SPH = " << p->U_SPH << std::endl
                           << "p->getXtuple() = " << p->getXtuple() << std::endl
                           << "p->isSurface = " << p->isSurface << std::endl
                           << "p->lap_U = " << p->lap_U << std::endl
                           << "p->gradP_SPH = " << p->gradP_SPH << std::endl;
                        throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
                    }

                    //理想では，この圧力の影響によって，粒子はソレノイダル条件を満たし，非圧縮条件を満たしていてほしい．
                    //その理想的な状態が続くなら，∇・Uは常に，ゼロである．仮想位置
                }
                std::cout << Green << "Elapsed time: " << Red << watch() << reset << " s\n";
#endif
                //@ ------------------------------------------------------ */
                //@                      粒子の時間発展                      */
                //@ ------------------------------------------------------ */
                Print("粒子の時間発展", Green);
                real_time = (P_RK_X.begin()->second)->gett();
                for (const auto &p : water_points)
                {
                    //位置
                    {
                        auto rk = P_RK_X[p];
                        // rk->displayStatus();
                        rk->push(p->U_SPH);
                        p->setX(rk->getX());
                    }
                    //速度
                    {
                        auto rk = P_RK_U[p];
                        /* メモ：
                            U* = U(n) + (ν∇^2U + g)*dt
                            U(n+1) = U* - 1/rho*∇P dtなので，
                            つまり，
                            U(n+1) = U(n) + (ν∇^2U + g)*dt - 1/rho*∇P dt
                            -> U(n+1) = U(n) + DUDt * dt
                            -> DUDt = ν∇^2U + g - 1/rho*∇P	または，(U* - U(n))/dt =  ν∇^2U + g を使って，
                            -> DUDt = (U* - U(n))*dt - 1/rho*∇P　とする．

                            next U  = p->U_SPH + p->DUDt_SPH *dt = -p->gradP_SPH / p->density * dt + p->tmp_U_SPH
                            と計算することからもわかるように．
                            確かに，p->DUDt_SPHは，下でいい．
                        */
                        // p->DUDt_SPH = -p->gradP_SPH / p->density + (p->tmp_U_SPH - p->U_SPH) / dt;
                        double nu = p->mu_SPH / p->density;

                        p->DUDt_SPH = -p->gradP_SPH / p->density + nu * p->lap_U + gravity;
                        /* ------------------------- 修正 ------------------------- */
                        p->repulsive_force_SPH *= 0.;
                        auto INTXN = IntersectionsSphereTrianglesLines(p->getContactFaces());
                        /*
                         * p->radius_SPH = C_SML * particle_spacing なので，
                         * p->radius_SPH / C_SML / 2.はparticle_spacingの半分程度となる．
                         */
                        double critical_distance = (p->radius_SPH / C_SML /*粒子間隔*/) / 3.;
                        for (const auto &[F0, F1, X, Y, N] : INTXN.getFFXYN(p, p->radius_SPH))
                        {
                            Tddd velocity = {0, 0, 0};
                            if (F1)
                                velocity = F0->getNormalTuple() * Dot(ToTddd(F0->getNetwork()->velocity), F0->getNormalTuple()) + F1->getNormalTuple() * Dot(ToTddd(F1->getNetwork()->velocity), F1->getNormalTuple());
                            else
                                velocity = F0->getNormalTuple() * Dot(ToTddd(F0->getNetwork()->velocity), F0->getNormalTuple());

                            if (Norm(X - p->getXtuple()) < critical_distance)
                            {
                                auto d = Norm(X - p->getXtuple());
                                auto r = critical_distance / (2. * M_PI);
                                auto w = (1. + std::tanh(-d / r + M_PI)) / 2.;
#ifdef WCSPH
                                p->repulsive_force_SPH += N * 0. * w;
#elif defined(EISPH)
                                p->repulsive_force_SPH += N * 5000. * w;
                                // p->repulsive_force_SPH += N * 5000. * ((std::tanh(2 * (critical_distance - M_PI / 2. * d) / critical_distance) + 1.) / 2.);
                                // p->DUDt_SPH += N * 2 * (critical_distance - Norm(X - p->getXtuple()) / critical_distance);
                                // auto relative_normal_velocity = Dot(p->U_SPH - velocity /*相対速度*/, N) * N;
                                // if (Dot(p->U_SPH - velocity /*相対速度*/, N) < 0 /*相対速度が壁向き*/)
                                // 	p->DUDt_SPH -= 100. * relative_normal_velocity * (1. - Norm(X - p->getXtuple()) / critical_distance);
#endif
                            }
                        }
                        /* ------------------------------------------------------ */
                        rk->push(p->DUDt_SPH + p->repulsive_force_SPH);
                        p->U_SPH = rk->getX();
                    }
                    //密度
#ifdef WCSPH
                    {
                        auto rk = P_RK_rho[p];
                        rk->push(p->DrhoDt_SPH);
                        p->setDensity(rk->getX()); //(*) -> (k+1)
                    }
#elif defined(EISPH)
                    p->setDensity(1000.); //(*) -> (k+1)
#endif
                }
            } while (!(P_RK_X[*water_points.begin()]->finished));

#ifdef EISPH
#ifndef use_RK_for_pressure
            //@ ------------------ 圧力の計算にルンゲクッタを使わない場合 ----------------- */
            for (const auto &p : water_points)
                p->pressure_SPH = p->pressure_SPH_;
                //@ ------------------------------------------------------------------- */
#endif
#endif
            //! ====================================================================================== */
            //! ====================================================================================== */
            //! ====================================================================================== */
            for (const auto &p : water_points)
            {
#ifdef WCSPH
                delete P_RK_rho[p];
#endif
                delete P_RK_U[p];
                delete P_RK_X[p];
#if defined(EISPH)
                delete P_RK_pressure[p];
#endif
            }

            for (const auto &p : water_points)
                if (!isFinite(p->U_SPH))
                    delete p;

            std::cout << Green << "Elapsed time: " << Red << watch() << reset << " s\n";
        }
    }
    catch (std::exception &e)
    {
        std::cerr << e.what() << reset << std::endl;
        throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
    };
    return 0;
};
