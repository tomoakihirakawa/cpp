#ifndef BEM_BVP_H
#define BEM_BVP_H

#include "Network.hpp"

//* ------------------------------------------------------ */
//*                        境界値問題を解く                   */
//* ------------------------------------------------------ */
// #define solve_equations_on_all_points
#define solve_equations_on_all_points_rigid_mode
#define solveBVP_debug

void addIG(std::unordered_map<netP *, Tdd> &P_phiphin, const netPp P_IN, const Tdd &igign)
{
    std::unordered_map<netP *, Tdd>::iterator it;
    if ((it = P_phiphin.find(P_IN)) != P_phiphin.end())
        std::get<0>(it->second) += std::get<0>(igign); // phiは忘れずに計算
    else
        P_phiphin[P_IN] = {std::get<0>(igign), 0.}; // phiは忘れずに計算
}

void addIGn(std::unordered_map<netP *, Tdd> &P_phiphin, const netPp P_IN, const Tdd &igign)
{
    std::unordered_map<netP *, Tdd>::iterator it;
    if ((it = P_phiphin.find(P_IN)) != P_phiphin.end())
        std::get<1>(it->second) += std::get<1>(igign); // phiは忘れずに計算
    else
        P_phiphin[P_IN] = {0., std::get<1>(igign)}; // phiは忘れずに計算
}

void addIG(std::map<netP *, Tdd> &P_phiphin, const netPp P_IN, const Tdd &igign)
{
    std::map<netP *, Tdd>::iterator it;
    if ((it = P_phiphin.find(P_IN)) != P_phiphin.end())
        std::get<0>(it->second) += std::get<0>(igign); // phiは忘れずに計算
    else
        P_phiphin[P_IN] = {std::get<0>(igign), 0.}; // phiは忘れずに計算
}

void addIGn(std::map<netP *, Tdd> &P_phiphin, const netPp P_IN, const Tdd &igign)
{
    std::map<netP *, Tdd>::iterator it;
    if ((it = P_phiphin.find(P_IN)) != P_phiphin.end())
        std::get<1>(it->second) += std::get<1>(igign); // phiは忘れずに計算
    else
        P_phiphin[P_IN] = {0., std::get<1>(igign)}; // phiは忘れずに計算
}

struct BEM_BVP
{
    std::vector<networkPoint *> Points;
    std::vector<networkFace *> Faces;
    const bool Neumann = false;
    const bool Dirichlet = true;

#if defined(use_lapack)
    lapack_lu *lu;
#else
    ludcmp_parallel *lu;
#endif
    VV_d mat_ukn;
    VV_d mat_kn;
    V_d knowns;
    // V_d unknowns;
    // std::vector<networkPoint *> vec_P;

    std::vector<pair_PB> vec_P;
    BEM_BVP(){};
    ~BEM_BVP() { delete this->lu; };

    void solveForPhiPhin_t() const
    {
        // b* ------------------------------------------------------ */
        // b*                         phiphin_t                      */
        // b* ------------------------------------------------------ */
        V_d knowns(vec_P.size());
        V_d phiORphin(vec_P.size());
        for (auto i = 0; i < this->vec_P.size(); i++)
        {
            auto [p, DorN] = vec_P[i];
            if (DorN == Dirichlet)
                knowns[i] = std::get<0>(p->phiphin_t);
            else
                knowns[i] = std::get<1>(p->phiphin_t);
        }
        this->lu->solve(Dot(mat_kn, knowns) /*既知のベクトル（右辺）*/, phiORphin /*解*/);
        /* ------------------------------------------------------ */
        for (auto i = 0; i < this->vec_P.size(); i++)
        {
            // auto p = this->vec_P[i];
            // if (p->Neumann)
            auto [p, DorN] = vec_P[i];
            if (DorN == Dirichlet)
                std::get<1>(p->phiphin_t) = phiORphin[i];
            else
                std::get<0>(p->phiphin_t) = phiORphin[i];
            p->pressure_BEM = -std::get<0>(p->phiphin_t) - _GRAVITY_ * p->height() - Dot(p->U_BEM, p->U_BEM) / 2.;
            p->pressure_BEM *= _WATER_DENSITY_;
            p->pressure = p->pressure_BEM;
        }
    };

    void solve(const Network &water, const Buckets<networkPoint> &FMM_BucketsPoints, const Buckets<networkFace> &FMM_BucketsFaces)
    {
        this->Points = ToVector(water.getPoints());
        this->Faces = ToVector(water.getFaces());
        using map_P_Vd = std::map<netP *, V_d>;
        //* ------------------------------------------------------ */
        //%                     各点で方程式を作る場合                */
        //* ------------------------------------------------------ */
        std::cout << "各点で方程式を作る場合" << std::endl;
        map_P_P_Tdd P_P_IGIGn;

        //@ 各バケツでのモーメントを次数別に保存する．(ユニーク) p->{k,m,Yn,Y}ベクトル
        using uo_P_uoTiiTdd = std::unordered_map<networkPoint *, std::unordered_map<Tii /*k,m*/, Tdd /*YYn*/>>;
        using V_uo_P_uoTiiTdd = std::vector<uo_P_uoTiiTdd>;
        using VV_uo_P_uoTiiTdd = std::vector<V_uo_P_uoTiiTdd>;
        using VVV_uo_P_uoTiiTdd = std::vector<VV_uo_P_uoTiiTdd>;

        VVV_uo_P_uoTiiTdd VVV_P_km_YnY(FMM_BucketsPoints.xsize, VV_uo_P_uoTiiTdd(FMM_BucketsPoints.ysize, V_uo_P_uoTiiTdd(FMM_BucketsPoints.zsize))); //途中でreverseしているので注意
        int ii = 0;
        //
        map_pairPB_pairPB_Tdd PDN_PDN_IGIGn;
        map_pairPB_Tdd init_PDN_IGIGn;
        /* ------------------------------------------------------ */
        /*                     init_PDN_IGIGn                     */
        /* ------------------------------------------------------ */
        for (const auto &org : Points)
        {
            if (org->CORNER)
            {
                init_PDN_IGIGn[{org, Dirichlet}] = {0., 0.};
                init_PDN_IGIGn[{org, Neumann}] = {0., 0.}; // b!
            }
            else if (org->Dirichlet)
                init_PDN_IGIGn[{org, Dirichlet}] = {0., 0.};
            else if (org->Neumann)
                init_PDN_IGIGn[{org, Neumann}] = {0., 0.};
        }
        /* ------------------------------------------------------ */
        /*                        P_P_IGIGn                       */
        /* ------------------------------------------------------ */
        for (const auto &org : Points)
        {
            P_P_IGIGn[org] = {}; //! P_P_IGIGn[row][col] = {IG,IGn}

            //! これで，PDN_PDN_IGIGnは，行の数がCORNERの分だけ倍になった．
            if (org->CORNER)
            {
                PDN_PDN_IGIGn[{org, Dirichlet}] = init_PDN_IGIGn;
                PDN_PDN_IGIGn[{org, Neumann}] = init_PDN_IGIGn; // b!
            }
            else if (org->Dirichlet)
                PDN_PDN_IGIGn[{org, Dirichlet}] = init_PDN_IGIGn;
            else if (org->Neumann)
                PDN_PDN_IGIGn[{org, Neumann}] = init_PDN_IGIGn;
        }
// #define single_layer_FMM
#if defined(single_layer_FMM)
        // double spacing = Mean(extLength(water.getLines())) * 10;
        // Buckets<networkFace> FMM_BucketsFaces(water.bounds(), spacing);
        // Buckets<networkPoint> FMM_BucketsPoints(water.bounds(), spacing);
        // for (const auto &f : water.getFaces())
        // 	FMM_BucketsFaces.add(f->getXtuple(), f);
        // for (const auto &p : water.getPoints())
        // 	FMM_BucketsPoints.add(p->getXtuple(), p);
        std::cout << "FMMのモーメントを計算" << std::endl;
        std::cout << "xsize = " << FMM_BucketsPoints.xsize << std::endl;
        std::cout << "ysize = " << FMM_BucketsPoints.ysize << std::endl;
        std::cout << "zsize = " << FMM_BucketsPoints.zsize << std::endl;

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (auto i = 0; i < FMM_BucketsPoints.buckets.size(); ++i)
#ifdef _OPENMP
#pragma omp single nowait
#endif
        {
            for (auto j = 0; j < FMM_BucketsPoints.buckets[i].size(); ++j)
            {
                for (auto k = 0; k < FMM_BucketsPoints.buckets[i][j].size(); ++k)
                {
                    // auto &P_km_YnY = VVV_P_km_YnY[i][j][k]; //@このバケツでのモーメント．次数別に保存m,k
                    for (const auto &F : FMM_BucketsFaces.buckets[i][j][k])
                    {
                        for (const auto &P : FMM_BucketsPoints.buckets[i][j][k])
                        {
                            Tddd center = FMM_BucketsFaces.indices2X_center(i, j, k);
                            for (const auto &[p, km_YYn] : BEM::calc_P_MomentQuadTuple(F, center)) //@このバケツでのモーメント．次数別に保存m,kに計算する
                            {
                                /*uoTiiTdd*/
                                VVV_P_km_YnY[i][j][k][p] += km_YYn;
                                // P_km_YnY[p] += km_YYn;
                            }
                        }
                    }
                }
            }
        }
        std::cout << "FMM実行" << std::endl;
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (auto i = 0; i < FMM_BucketsPoints.buckets.size(); ++i)
            for (auto j = 0; j < FMM_BucketsPoints.buckets[i].size(); ++j)
                for (auto k = 0; k < FMM_BucketsPoints.buckets[i][j].size(); ++k)
                {
                    for (const auto &P : FMM_BucketsPoints.buckets[i][j][k])
#ifdef _OPENMP
#pragma omp single nowait
#endif
                    {
                        // std::cout << "(i,j,k) = (" << i << "," << j << "," << k << ")" << std::endl;
                        auto &P_igign = P_P_IGIGn[P];
                        for (auto I = 0; I < FMM_BucketsFaces.buckets.size(); ++I)
                            for (auto J = 0; J < FMM_BucketsFaces.buckets[I].size(); ++J)
                                for (auto K = 0; K < FMM_BucketsFaces.buckets[I][J].size(); ++K)
                                {
                                    if ((i - 1 <= I && I <= i + 1) || (j - 1 <= J && J <= j + 1) || (k - 1 <= K && K <= k + 1))
                                    {
                                        //@ この点に，比較的近い面の積分
                                        // std::cout << "原点を節点にとり，方程式を作成．並列化" << std::endl;
                                        for (const auto &F : FMM_BucketsFaces.buckets[I][J][K])
                                        {
                                            for (const auto &[p, igign] : BEM::calc_P_IGIGnQuadTuple_mod(F, P))
                                            {
                                                addIG(P_igign, p, igign);
                                                if (p != P)
                                                {
                                                    addIGn(P_igign, p, igign);
                                                    //@ 対角成分は，リジッドモード法を使って計算できる
                                                    addIGn(P_igign, P, -igign);
                                                }
                                            }
                                        }
                                    }
                                    else
                                    { //@ この点に，比較的遠い面の積分
                                        Tddd O = FMM_BucketsFaces.indices2X_center(i, j, k);
                                        for (const auto &[p, kmYYn] : VVV_P_km_YnY[i][j][k])
                                        {
                                            Tdd igign = {0., 0.};
                                            for (const auto &[km, YYn] : kmYYn)
                                            {
                                                auto [k, m] = km;
                                                auto [nr, theta, psi] = ToSphericalCoodrinates(p->getXtuple() - O);
                                                igign += YYn * real_sph_scale_ommited(k, m, theta, psi) / std::pow(nr, k + 1);
                                            }
                                            addIG(P_igign, p, igign);
                                            if (p != P)
                                            {
                                                addIGn(P_igign, p, igign);
                                                //@ 対角成分は，リジッドモード法を使って計算できる
                                                addIGn(P_igign, P, -igign);
                                            }
                                        }
                                    }
                                }
                    }
                }
#else

        /* ------------------------------------------------------ */

        // #define quad_element
#define linear_element
        // #define liear_and_quad_element

#ifdef _OPENMP
        std::cout << "原点を節点にとり，方程式を作成．並列化" << std::endl;
#ifdef quad_element
        std::cout << "quad" << std::endl;
#elif defined(linear_element)
        std::cout << "linear" << std::endl;
#endif

#pragma omp parallel
#endif
        for (auto &[p0_bool, P_igign] : PDN_PDN_IGIGn)
#ifdef _OPENMP
#pragma omp single nowait
#endif
        {
            auto p0 = std::get<0>(p0_bool);
            std::map<pair_PB, Tdd>::iterator it;
            double p0_ign = 0.;
            pair_PB p_bool;

            for (const auto &integrating_f : Faces)
            {

                /*         *   *       */
                /*        / \ / \      */
                /*       *---*---*     */
                /*      / \ /F\ / \    */
                /*     *---*---*---*   */
                /*        / \ / \      */
                /*       *---*---*     */
                /*
                積分は各節点の変数と重みによって表される．
                その重みを与える．
                */
#if defined(quad_element)
                for (const auto &[p, igign] : BEM::calc_P_IGIGnQuadTuple_mod(integrating_f, p0))
                {
                    if (p->CORNER)
                    {
                        if (integrating_f->Dirichlet)
                            p_bool = {p, Dirichlet};
                        else if (integrating_f->Neumann)
                            p_bool = {p, Neumann};
                    }
                    else
                        p_bool = {p, p->Dirichlet};

                    P_igign[p_bool] += igign;

                    if (p != p0)
                        p0_ign -= std::get<1>(igign);
                }
#elif defined(linear_element)
                auto func = [&](const std::tuple<netP *, Tdd> &p_igign)
                {
                    const auto [p, igign] = p_igign;
                    if (p->CORNER)
                    {
                        if (integrating_f->Dirichlet)
                            P_igign[{p, Dirichlet}] += igign;
                        else // if (integrating_f->Neumann)
                            P_igign[{p, Neumann}] += igign;
                    }
                    else
                        P_igign[{p, p->Dirichlet}] += igign;
                    if (p != p0)
                        p0_ign -= std::get<1>(igign);
                };
                const auto [q0igign, q1igign, q2igign] = BEM::calc_P_IGIGnLinear3Tuples(integrating_f, p0);
                func(q0igign);
                func(q1igign);
                func(q2igign);
#elif defined(linear_element2) for (const auto &[p, igign] : BEM::calc_P_IGIGnLinearTuple(integrating_f, p0))
                {
                    if (p->CORNER)
                    {
                        if (integrating_f->Dirichlet)
                            p_bool = {p, Dirichlet};
                        else if (integrating_f->Neumann)
                            p_bool = {p, Neumann};
                    }
                    else
                        p_bool = {p, p->Dirichlet};

                    P_igign[p_bool] += igign;

                    if (p != p0)
                        p0_ign -= std::get<1>(igign);
                }
#endif
            }
            std::get<1>(P_igign[p0_bool]) = p0_ign;
        }
#endif
        std::cout << "並列化 DONE" << std::endl;
        std::cout << "2つの係数行列の情報を持つ　P_P_IGIGn　を境界条件に応じて入れ替える（移項）:" << std::endl;

        mat_ukn = VV_d(PDN_PDN_IGIGn.size(), V_d(PDN_PDN_IGIGn.size(), 0.)); //! ok
        mat_kn = VV_d(PDN_PDN_IGIGn.size(), V_d(PDN_PDN_IGIGn.size(), 0.));  //! ok

        // 順番を間違えないようにベクトルを作成
        vec_P = std::vector<pair_PB>(PDN_PDN_IGIGn.size());
        auto i = 0;
        for (auto it = PDN_PDN_IGIGn.begin(); it != PDN_PDN_IGIGn.end(); ++it)
            vec_P[i++] = it->first;

        {
            // b@ ------------------------------------------------------ */
            // b@                 系数行列mat_ukn．mat_knの計算             */
            // b@ ------------------------------------------------------ */

#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (auto i = 0; i < vec_P.size(); ++i)
            {
                /*@ mat_ukn, mat_kn, vec_P
                    +--+--+--+
                |	|--+--+--|
                V	|--+--+--|
                    +--+--+--+
                */
                auto PDN_IGIGn = PDN_PDN_IGIGn.at(vec_P[i]);
                auto p = std::get<0>(vec_P[i]);
                auto DorN = std::get<1>(vec_P[i]);
                /* mat_ukn, mat_kn, vec_P
                       ====>
                    +--+--+--+
                    |--+--+--|
                    |--+--+--|
                    +--+--+--+
                */
                if (!(p->CORNER && DorN == Neumann /*変更する対象の行*/)) //! OK
                {
                    for (auto j = 0; j < vec_P.size(); ++j) //! OK
                    {
                        auto igign = PDN_IGIGn.at(vec_P[j]);
                        /* --------------------------------------- */
                        // 未知変数の係数行列は左，既知変数の係数行列は右
                        if (std::get<1>(vec_P[j]) == Neumann)
                            igign = {-std::get<1>(igign), -std::get<0>(igign)};
                        //% IGIGn は 左辺に IG*φn が右辺に IGn*φ が来るように計算しているため，移項する場合，符号を変える必要がある．
                        /*
                        IG*φn - IGn*φ = 0
                        */
                        /* --------------------------------------- */
                        mat_ukn[i][j] = std::get<0>(igign);
                        mat_kn[i][j] = std::get<1>(igign);
                    }
                }
            }

            /* ------------------------------------------------------ */
            double maxpp = 0;
            for (auto i = 0; i < mat_ukn.size(); i++)
            {
                // LUするのはmat_uknだけなので，mat_knの最大値を使う必要はない
                //  if (maxpp < std::abs(mat_kn[i][i]))
                //  	maxpp = std::abs(mat_kn[i][i]);
                if (maxpp < std::abs(mat_ukn[i][i]))
                    maxpp = std::abs(mat_ukn[i][i]);
            }
            /* ------------------------------------------------------ */
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (auto i = 0; i < vec_P.size(); ++i)
            {
                auto PDN_IGIGn = PDN_PDN_IGIGn.at(vec_P[i]);
                auto p = std::get<0>(vec_P[i]);
                auto DorN = std::get<1>(vec_P[i]);
                if (p->CORNER && DorN == Neumann /*変更する対象の行*/) //! OK
                {
                    for (auto j = 0; j < vec_P.size(); ++j)
                    {
                        auto q = std::get<0>(vec_P[j]);
                        if (p == q)
                        {
                            if (std::get<1>(vec_P[j]) == Neumann)
                            {
                                mat_ukn[i][j] = maxpp; //φの系数
                                mat_kn[i][j] = 0;      //φnの系数
                            }
                            else
                            {
                                mat_ukn[i][j] = 0;    //φnの系数
                                mat_kn[i][j] = maxpp; //φの系数移行したからマイナス？　いいえ，移項を考慮した上でこれでいい．
                            }
                        }
                        else
                        {
                            mat_ukn[i][j] = 0;
                            mat_kn[i][j] = 0;
                        }
                    }
                }
            }

            if (!isFinite(mat_ukn))
                throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "mat_ukn is not finite");
            if (!isFinite(mat_kn))
                throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "mat_kn is not finite");

            // b@ ------------------------------------------------------ */
        }

        {
            // b$ ------------------------------------------------------ */
            // b$　　　                   knownsの計算　                    */
            // b$ ------------------------------------------------------ */
            /**
             * Dot(mat_ukn,phiORphin) = Dot(mat_kn,knowns)
             * => phiORphin = Dot(mat_ukn^-1, Dot(mat_kn,knowns))
             */
            knowns = V_d(vec_P.size()); //! ok
                                        //! Pointsの順番と合わせてとるように注意
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (auto i = 0; i < vec_P.size(); i++)
            {
                auto p = std::get<0>(vec_P[i]);
                auto DorN = std::get<1>(vec_P[i]);
                if (DorN == Dirichlet) //! OK
                    knowns[i] = p->phi_Dirichlet = std::get<0>(p->phiphin);
                else
                    knowns[i] = p->phin_Neumann; // std::get<1>(p->phiphin);
            }

            // b$ ------------------------------------------------------ */
        }

        {
            // b% ------------------------------------------------------ */
            // b%                       境界積分方程式を解く                 */
            // b% ------------------------------------------------------ */
            std::cout << "--------------------- 境界積分方程式を解く ---------------------" << std::endl;
            /* ------------------------------------------------------ */
            V_d phiORphin(knowns.size());
            std::cout << "phiORphin.size()= " << phiORphin.size() << std::endl;
            std::cout << "  mat_kn.size() = " << mat_kn.size() << std::endl;
#if defined(solve_equations_on_all_points) || defined(solve_equations_on_all_points_RBF) || defined(solve_equations_on_all_points_rigid_mode)
            //* 未知変数の計算
#if defined(use_lapack)
            std::cout << "lapack lu decomposition" << std::endl;
            this->lu = new lapack_lu(mat_ukn /*未知の行列係数（左辺）*/);
#else
            std::cout << "parallel lu decomposition" << std::endl;
            this->lu = new ludcmp_parallel(mat_ukn /*未知の行列係数（左辺）*/);
#endif
            std::cout << "try to solve" << std::endl;
            this->lu->solve(Dot(mat_kn, knowns) /*既知のベクトル（右辺）*/, phiORphin /*解*/);
#else
            //* 未知変数の計算
            std::cout << "SVD decomposition" << std::endl;
            SVD svd(mat_ukn);
            svd.solve(Dot(mat_kn, knowns) /*既知のベクトル（右辺）*/, phiORphin /*解*/);
#endif
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (auto i = 0; i < vec_P.size(); ++i)
            {
                auto p = std::get<0>(vec_P[i]);
                auto DorN = std::get<1>(vec_P[i]);
                if (DorN == Dirichlet)
                    std::get<1>(p->phiphin) = p->phin_Dirichlet = phiORphin[i];
                else
                    std::get<0>(p->phiphin) = p->phi_Neumann = phiORphin[i];
            }
            if (!isFinite(phiORphin))
                throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "phiORphin is not finite");
        }
    };
};

#endif