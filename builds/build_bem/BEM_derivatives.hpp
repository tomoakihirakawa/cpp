#ifndef BEM_derivatives_H
#define BEM_derivatives_H

#include "Network.hpp"

Tddd fitToNeumannVelocity(Tddd VECTOR, const networkPoint *const p)
{
    if (p->Neumann || p->CORNER)
    {
        Tddd nf, vn0, vn1;
        double max_w = kernel_Bspline3(0., p->radius);
        for (const auto &[f, hit_X] : Reverse(p->getContactFacesXCloser()) /*遠い方から*/)
        {
            nf = f->getNormalTuple();
            vn0 = Dot(VECTOR, nf) * nf;
            vn1 = Dot(f->getNetwork()->velocityRigidBody(hit_X), nf) * nf;
            VECTOR += kernel_Bspline3(Norm(hit_X - p->getXtuple()), p->radius) / max_w * (vn1 - vn0);
        }
    }
    return VECTOR;
};

Tddd condition_Ua(Tddd VECTOR, const networkPoint *const p)
{
    if (p->CORNER)
    {
        auto cross = Normalize(Cross(p->getNormalNeumann_BEM(), p->getNormalDirichlet_BEM()));
        VECTOR = Dot(VECTOR, cross) * cross;
        // for (const auto &l : p->getLines())
        // 	if (l->CORNER)
        // 	{
        // 		Tddd dir = Normalize((*l)(p)->getXBuffer() - p->getXBuffer());
        // 		VECTOR = Dot(VECTOR, dir) * dir;
        // 	}
        for (const auto &f : p->getFaces())
            if (f->Neumann)
                VECTOR -= Dot(VECTOR, f->getNormalTuple()) * f->getNormalTuple();
        return VECTOR;
    }
    else if (p->Dirichlet)
    {
        return VECTOR - Dot(VECTOR, p->getNormal_BEM()) * p->getNormal_BEM();
    }
    else
    {
        for (const auto &f : p->getFaces())
            VECTOR -= Dot(VECTOR, f->getNormalTuple()) * f->getNormalTuple();
        return VECTOR;
    }
};

Tddd vectorsToSurfaceFromBufferX(const networkPoint *p, const std::vector<T3Tddd> &next_Vrtx)
{
    auto closestXFacing = [](const Tddd &p_next_X, const double radius, const std::vector<T3Tddd> &vertices, const Tddd &n)
    {
        Tddd r = {1E+100, 1E+100, 1E+100};
        for (const auto &vertex : vertices)
        {
            if (isFacing(TriangleNormal(vertex), n, M_PI / 180 * 20))
            {
                auto intxn = IntersectionSphereTriangle_(p_next_X, radius, vertex);
                if (intxn.isIntersecting)
                    if (Norm(r) >= Norm(intxn.X - p_next_X))
                        r = intxn.X - p_next_X;
            }
        }
        return r;
    };

    if (next_Vrtx.empty())
        return {0., 0., 0.};
    else if (p->Neumann || p->CORNER)
    {
        std::vector<Tddd> F_clings;
        for (const auto &f : p->getFacesNeumann())
        {
            // auto [p0, p1, p2] = f->getPointsTuple();
            // auto n = TriangleNormal(p0->getXBuffer() + p0->U_BUFFER, p1->getXBuffer() + p1->U_BUFFER, p2->getXBuffer() + p2->U_BUFFER);
            auto n = f->getNormalTuple();
            auto to_closest_X = closestXFacing(p->getXBuffer() + p->U_BUFFER, p->radius, next_Vrtx, n);
            if (isFinite(to_closest_X))
                F_clings.push_back(Dot(to_closest_X, n) * n);
        }
        Tddd r = Mean(F_clings);
        if (isFinite(r))
            return r;
        else
            return {0., 0., 0.};
    }
    else
        return {0., 0., 0.};
};

double minViewRatio(const networkPoint *const p)
{
    double a = p->getSolidAngle();
    return (2 * M_PI - Min(Tdd{std::abs(a - 2 * M_PI), std::abs(2 * M_PI - a)})) / (2 * M_PI);
};

double normalVariance(const networkPoint *const p)
{
    auto n = p->getNormalDirichlet_BEM();
    double m = 0, s = 0;
    for (const auto &f : p->getFacesDirichlet())
    {
        m += (M_PI / 2. - VectorAngle(n, f->getNormalTuple())) / (M_PI / 2.);
        s += 1;
    }
    return m / s;
};

Tddd vectorTangentialShift(const networkPoint *p)
{
    auto nextX_U_Ua = [](const networkPoint *p)
    {
        return p->getXBuffer() + p->U_BUFFER;
    };

    auto next_length = [nextX_U_Ua](const networkLine *const l)
    {
        auto [p0, p1] = l->getPointsTuple();
        return Norm(nextX_U_Ua(p0) - nextX_U_Ua(p1));
    };

    auto getBaseLength = [next_length](const networkLine *line)
    {
        auto [p0, p1] = line->getPointsTuple();
        std::unordered_set<networkLine *> lc = Join(extLinesCORNER_(p0->getFaces()), extLinesCORNER_(p1->getFaces()));
        if (!lc.empty())
        {
            V_d ll;
            for (const auto &l : lc)
                ll.emplace_back(next_length(l));
            return Mean(ll);
        }
        else
        {
            // pを引っ張る力は，ノイマン面とディリクレ面で干渉しない
            auto [p0, p1] = line->getPointsTuple();
            double m = 1, s = 0;
            for (const auto &l : Join(p0->getLinesAround(), p1->getLinesAround()))
                if (line != l)
                    if ((line->Dirichlet && (l->Dirichlet || l->CORNER)) || (line->Neumann && (l->Neumann || l->CORNER)) || (line->CORNER && l->CORNER))
                    {
                        m *= next_length(l);
                        s += 1;
                    }
            return std::pow(m, 1. / s);
        }
    };

    auto vectorToNextNeighborsCenter = [nextX_U_Ua](const networkPoint *const p)
    {
        double s = 0;
        Tddd ret = {0., 0., 0.};
        Tddd pX = nextX_U_Ua(p);
        /* ------------------------------------------------------ */
        if (p->CORNER)
        {
            for (const auto &l : p->getLines())
                if (l->CORNER)
                {
                    ret += (nextX_U_Ua(((*l)(p))) - pX);
                    s += 1;
                }
        }
        else
        {
            for (const auto &l : p->getLines())
            {
                ret += (nextX_U_Ua(((*l)(p))) - pX);
                s += 1;
            }
        }
        return ret / s;
    };

    /* ------------------------------------------------------ */
    Tddd V = {0., 0., 0.};
    auto faces = p->getFaces();
    for (const auto &f : faces)
    {
        auto [p0, p1, p2] = f->getPointsTuple(p);
        double d0, d1, d2;
        Tddd r0, r1, r2;
        {
            auto intersect = IntersectionSphereLine(nextX_U_Ua(p0), 1E+20, T2Tddd{nextX_U_Ua(p1), nextX_U_Ua(p2)});
            d0 = intersect.distance;
            r0 = intersect.X - nextX_U_Ua(p0);
        }
        {
            auto intersect = IntersectionSphereLine(nextX_U_Ua(p2), 1E+20, T2Tddd{nextX_U_Ua(p0), nextX_U_Ua(p1)});
            d1 = intersect.distance;
            // r1 = intersect.X - nextX_U_Ua(p2);
        }
        {
            auto intersect = IntersectionSphereLine(nextX_U_Ua(p1), 1E+20, T2Tddd{nextX_U_Ua(p2), nextX_U_Ua(p0)});
            d2 = intersect.distance;
            // r2 = intersect.X - nextX_U_Ua(p1);
        }
        V += (d0 - (d0 + d1 + d2) / 3.) * r0;
    }
    /* ------------------------------------------------------ */
    Tddd pX = nextX_U_Ua(p);
    // EMTは節点をおしすぎることがあるようだ
    // やはり接線方向でないといけないようだ
    double c_LS = 0.5 /*0.1~0.5*/, c_EMT = 0.1;
    Tddd V_EMT = {0., 0., 0.};
    auto a = minViewRatio(p);
    double tmp = 1, s = 0;
    auto V_LS = vectorToNextNeighborsCenter(p);
    if (p->CORNER)
    {
        if (a > 1. / 6.)
        {
            return condition_Ua(V_LS * c_LS, p);
        }
        else
            return {0., 0., 0.};
    }
    if (p->Dirichlet)
    {
        for (const auto &l : p->getLines())
            V_EMT += (next_length(l) - getBaseLength(l)) * Normalize(nextX_U_Ua((*l)(p)) - pX);
        return condition_Ua(V_EMT * c_EMT + 0.01 * V + V_LS * c_LS, p);
    }
    else
    {
        if (a > 2. / 3.) //比較的滑らかなノイマン面
        {
            for (const auto &l : p->getLines())
                V_EMT += (next_length(l) - getBaseLength(l)) * Normalize(nextX_U_Ua((*l)(p)) - pX);
            return condition_Ua(V_EMT * c_EMT + 0.01 * V + V_LS * c_LS, p);
        }
        else //比較的滑らかでないノイマン面
        {
            V_netFp f;
            for (const auto &l : p->getLines())
            {
                f = l->getFaces();
                a = VectorAngle(f[0]->getNormalTuple(), f[1]->getNormalTuple()) / (2 * M_PI);
                V_EMT += a * (nextX_U_Ua((*l)(p)) - pX);
            }
            return condition_Ua(V_EMT, p);
        }
    } // else
      // {
      // 	// for (const auto &l : p->getLines())
      // 	// 	V_EMT += (next_length(l) - getBaseLength(l)) * Normalize(nextX_U_Ua((*l)(p)) - pX);
      // 	// V_EMT = condition_Ua(V_EMT * c_EMT, p);

    // 	double tmp = 1, s = 0;
    // 	for (const auto &l : p->getLinesOppsoite())
    // 	{
    // 		auto [p0, p1] = l->getPointsTuple();
    // 		auto intersect = IntersectionSphereLine(nextX_U_Ua(p), 1E+20, T2Tddd{nextX_U_Ua(p0), nextX_U_Ua(p1)});
    // 		tmp *= intersect.distance;
    // 		s += 1;
    // 	}
    // 	for (const auto &l : p->getLinesOppsoite())
    // 	{
    // 		auto [p0, p1] = l->getPointsTuple();
    // 		auto intersect = IntersectionSphereLine(nextX_U_Ua(p), 1E+20, T2Tddd{nextX_U_Ua(p0), nextX_U_Ua(p1)});
    // 		V_EMT += (intersect.distance - std::pow(tmp, 1. / s)) * Normalize(intersect.X - pX);
    // 	}

    // 	return condition_Ua(V_EMT * c_EMT + V_LS, p);
    // 	// if (a > 2. / 3.) //比較的滑らかなノイマン面
    // 	// {
    // 	// 	// for (const auto &l : p->getLines())
    // 	// 	// 	V_EMT += (next_length(l) - getBaseLength(l)) * Normalize(nextX_U_Ua((*l)(p)) - pX);
    // 	// 	// V_EMT = condition_Ua(V_EMT * c_EMT, p);

    // 	// 	double tmp = 1, s = 0;
    // 	// 	for (const auto &l : p->getLinesOppsoite())
    // 	// 	{
    // 	// 		auto [p0, p1] = l->getPointsTuple();
    // 	// 		auto intersect = IntersectionSphereLine(nextX_U_Ua(p), 1E+20, T2Tddd{nextX_U_Ua(p0), nextX_U_Ua(p1)});
    // 	// 		tmp *= intersect.distance;
    // 	// 		s += 1;
    // 	// 	}
    // 	// 	for (const auto &l : p->getLinesOppsoite())
    // 	// 	{
    // 	// 		auto [p0, p1] = l->getPointsTuple();
    // 	// 		auto intersect = IntersectionSphereLine(nextX_U_Ua(p), 1E+20, T2Tddd{nextX_U_Ua(p0), nextX_U_Ua(p1)});
    // 	// 		V_EMT += (intersect.distance - std::pow(tmp, 1. / s)) * Normalize(intersect.X - pX);
    // 	// 	}

    // 	// 	return condition_Ua(V_EMT * c_EMT + V_LS, p);
    // 	// }
    // 	// else //比較的滑らかでないノイマン面
    // 	// {
    // 	// 	V_netFp f;
    // 	// 	for (const auto &l : p->getLines())
    // 	// 	{
    // 	// 		f = l->getFaces();
    // 	// 		a = VectorAngle(f[0]->getNormalTuple(), f[1]->getNormalTuple()) / (2 * M_PI);
    // 	// 		V_EMT += a * (nextX_U_Ua((*l)(p)) - pX);
    // 	// 	}
    // 	// 	return condition_Ua(V_EMT, p);
    // 	// }
    // }
};

std::tuple<Tddd, double, double, networkLine *, networkFace *> vectorToNearestAjacentSurface(const networkPoint *p)
{
    /*
    UartificialClingは，完璧にOmega(t+\delta t)に張り付くようにしなければ，
    面からはなれることで計算の破綻を招く可能性がある．
    */
    std::tuple<Tddd, double, double, networkLine *, networkFace *> ret = {{0, 0, 0}, 0., 0., nullptr, nullptr};
    Tddd pX = p->getXBuffer() + p->U_BUFFER;
    if (p->CORNER)
    {
        // auto a = minViewRatio(p);
        // if (a > 1 / 3)
        for (const auto &l : p->getLines())
            if (l->CORNER)
            {
                auto intxn = IntersectionSphereLine(pX, 1E+20, T2Tddd{p->getXBuffer(), (*l)(p)->getXBuffer()});
                if (Norm(std::get<0>(ret) - pX) >= Norm(intxn.X - pX))
                {
                    std::get<0>(ret) = intxn.X;
                    std::get<1>(ret) = intxn.t;
                    std::get<2>(ret) = 1 - intxn.t;
                    std::get<3>(ret) = l;
                    std::get<4>(ret) = nullptr;
                }
            }
    }
    else
    {
        for (const auto &f : p->getFaces())
        {
            if ((p->Dirichlet && f->Dirichlet) || (p->Neumann && f->Neumann))
            {
                auto [p0, p1, p2] = f->getPointsTuple(p);
                auto intxn = IntersectionSphereTriangle_(pX, 1E+20, T3Tddd{p0->getXBuffer(), p1->getXBuffer(), p2->getXBuffer()});
                if (Norm(std::get<0>(ret) - pX) >= Norm(intxn.X - pX))
                {
                    std::get<0>(ret) = intxn.X;
                    std::get<1>(ret) = intxn.t0;
                    std::get<2>(ret) = intxn.t1;
                    std::get<3>(ret) = nullptr;
                    std::get<4>(ret) = f;
                }
            }
        }
    }
    return ret;
};

void calculateVectorToSurfaceInBuffer(const Network &net, const bool adjust_dirichlet = true)
{
    /*
    @ この方法なら，次の時刻における任意の場所でのポテンシャルを見積もることができる．
    @ このことは，任意のノイマン面上に節点を維持する上で便利である．
    @ Ω(t+δt)をまず見積もり，その面上で最適な格子配置となるように流速を修正する．
    */
    auto Points = ToVector(net.getPoints());
    for (const auto &p : Points)
        p->U_BUFFER = p->U_BUFFER_BUFFER = {0., 0., 0.};

    //@ ------------------------------------------------------ */
    //@           次の時刻で最適な格子を目指す修正流速を計算          */
    //@ ------------------------------------------------------ */
    for (auto kk = 0; kk < 1; ++kk)
    {
        //% ------------------------------------------------------ */
        //%           　　　　 vectorTangentialShift   　 　         */
        //%          ラプラス平滑化と引っ張り合わせた接線方向にシフト      */
        //% ------------------------------------------------------ */
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (const auto &p : Points)
#ifdef _OPENMP
#pragma omp single nowait
#endif
            if (!p->Dirichlet || (p->Dirichlet && adjust_dirichlet))
            {
                p->U_BUFFER_BUFFER = vectorTangentialShift(p);
            }

        for (const auto &p : Points)
        {
            if (isFinite(p->U_BUFFER_BUFFER))
                p->U_BUFFER += p->U_BUFFER_BUFFER;
            p->U_BUFFER_BUFFER = {0., 0., 0.};
        }

        //% ------------------------------------------------------ */
        //%              vectorToNearestAjacentSurface             */
        //%           　　　   周辺ディリクレ面に移動        　　　　    */
        //% ------------------------------------------------------ */

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (const auto &p : Points)
#ifdef _OPENMP
#pragma omp single nowait
#endif
            if (p->CORNER || (p->Dirichlet && adjust_dirichlet))
            {
                if (isFinite(p->U_BUFFER_BUFFER))
                {
                    p->clungSurface = vectorToNearestAjacentSurface(p);
                    p->U_BUFFER_BUFFER = (std::get<0>(p->clungSurface) - (p->getXBuffer() + p->U_BUFFER));
                    //! 角のノイマン面から離れるのを防ぐ
                    if (p->CORNER)
                        for (const auto &f : p->getFacesNeumann())
                            p->U_BUFFER_BUFFER -= Dot(p->U_BUFFER_BUFFER, f->getNormalTuple()) * f->getNormalTuple();
                }
            }

        for (const auto &p : Points)
            if (p->CORNER || (p->Dirichlet && adjust_dirichlet))
            {
                if (isFinite(p->U_BUFFER_BUFFER, 1E+10))
                    p->U_BUFFER += p->U_BUFFER_BUFFER;
                p->U_BUFFER_BUFFER = {0., 0., 0.};
            }

            //% ------------------------------------------------------ */
            //%             vectorsToSurfaceFromBufferX                */
            //%              　　　近傍のノイマン面へ移動           　　　   */
            //% ------------------------------------------------------ */

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (const auto &p : Points)
#ifdef _OPENMP
#pragma omp single nowait
#endif
            if (p->CORNER || p->Neumann)
            {
                {
                    //接触面候補の次の時刻の位置を予測
                    std::unordered_set<networkFace *> Fs = p->getContactFaces();
                    for (auto &f : p->getContactFaces())
                        for (auto &q : f->getPoints())
                            for (auto &F : q->getFaces())
                                Fs.emplace(F);

                    std::vector<T3Tddd> nextBodyVertices;
                    for (auto &f : Fs)
                    {
                        auto [p0, p1, p2] = f->getPointsTuple();
                        // auto X0 = p0->getXtuple() + f->getNetwork()->velocityRigidBody(p0->getXtuple()) * p->RK_X.getdt();
                        // auto X1 = p1->getXtuple() + f->getNetwork()->velocityRigidBody(p1->getXtuple()) * p->RK_X.getdt();
                        // auto X2 = p2->getXtuple() + f->getNetwork()->velocityRigidBody(p2->getXtuple()) * p->RK_X.getdt();
                        /* ------------------------------------------------------ */
                        auto net = f->getNetwork();
                        auto COM = net->RK_COM.getX(net->velocityTranslational());
                        Quaternion q;
                        q = q.d_dt(net->velocityRotational());
                        auto Q = net->RK_Q.getX(q());
                        auto X0 = RigidBodyMove(p0, COM, Q);
                        auto X1 = RigidBodyMove(p1, COM, Q);
                        auto X2 = RigidBodyMove(p2, COM, Q);
                        /* ------------------------------------------------------ */
                        nextBodyVertices.emplace_back(T3Tddd{X0, X1, X2});
                    }
                    p->U_BUFFER_BUFFER = vectorsToSurfaceFromBufferX(p, nextBodyVertices);
                    //! 角のディリクレ面へのめり込みを防止
                    if (p->CORNER)
                        p->U_BUFFER_BUFFER -= Dot(p->U_BUFFER_BUFFER, p->getNormalDirichlet_BEM()) * p->getNormalDirichlet_BEM();
                }
            }

        for (const auto &p : Points)
        {
            if (isFinite(p->U_BUFFER_BUFFER, 1E+10))
                p->U_BUFFER += p->U_BUFFER_BUFFER;
            p->U_BUFFER_BUFFER = {0., 0., 0.};
        }
    }

    std::cout << "p->U_BUFFER finished" << std::endl;
};

void calculateVectorFromBufferToContactFaces(const Network &net)
{
    /*
    @ ノイマン面に貼り付けるための必要な調整
    */
    auto Points = ToVector(net.getPoints());
    // for (const auto &p : Points)
    // 	p->U_BUFFER = p->U_BUFFER_BUFFER = {0., 0., 0.};
    for (auto kk = 0; kk < 50; ++kk)
    {
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (const auto &p : Points)
#ifdef _OPENMP
#pragma omp single nowait
#endif
        {
            if (p->CORNER || p->Neumann)
            {
                //接触面候補の次の時刻の位置を予測
                std::unordered_set<networkFace *> Fs = p->getContactFaces();
                for (auto &f : p->getContactFaces())
                    for (auto &q : f->getPoints())
                        for (auto &F : q->getFaces())
                            Fs.emplace(F);

                std::vector<T3Tddd> nextBodyVertices;
                for (auto &f : Fs)
                {
                    auto [p0, p1, p2] = f->getPointsTuple();

                    // auto X0 = p0->getXtuple() + f->getNetwork()->velocityRigidBody(p0->getXtuple()) * p->RK_X.getdt();
                    // auto X1 = p1->getXtuple() + f->getNetwork()->velocityRigidBody(p1->getXtuple()) * p->RK_X.getdt();
                    // auto X2 = p2->getXtuple() + f->getNetwork()->velocityRigidBody(p2->getXtuple()) * p->RK_X.getdt();
                    /* ------------------------------------------------------ */
                    auto net = f->getNetwork();
                    auto COM = net->RK_COM.getX(net->velocityTranslational());
                    Quaternion q;
                    q = q.d_dt(net->velocityRotational());
                    auto Q = net->RK_Q.getX(q());
                    auto X0 = RigidBodyMove(p0, COM, Q);
                    auto X1 = RigidBodyMove(p1, COM, Q);
                    auto X2 = RigidBodyMove(p2, COM, Q);

                    nextBodyVertices.emplace_back(T3Tddd{X0, X1, X2});
                }
                //% ------------------------------------------------------ */
                //%              vectorsToSurfaceFromBufferX              */
                //%              　　　近傍のノイマン面へ移動           　　　   */
                //% ------------------------------------------------------ */
                p->U_BUFFER_BUFFER = vectorsToSurfaceFromBufferX(p, nextBodyVertices);
                //! 角のディリクレ面へのめり込みを防止
                if (p->CORNER)
                    p->U_BUFFER_BUFFER -= Dot(p->U_BUFFER_BUFFER, p->getNormalDirichlet_BEM()) * p->getNormalDirichlet_BEM();
            }
        }

        for (const auto &p : Points)
        {
            if (isFinite(p->U_BUFFER_BUFFER))
                p->U_BUFFER += 0.1 * p->U_BUFFER_BUFFER;
            p->U_BUFFER_BUFFER = {0., 0., 0.};
        }
    }

    for (const auto &p : Points)
    {
        if (!isFinite(p->U_BUFFER))
        {
            std::cout << "p->RK_X.getdt() = " << p->RK_X.getdt() << std::endl;
            std::cout << "p->U_update_BEM = " << p->U_update_BEM << std::endl;
            std::cout << "p->U_BUFFER = " << p->U_BUFFER << std::endl;
            std::cout << "p->Dirichlet = " << p->Dirichlet << std::endl;
            std::cout << "p->Neumann = " << p->Neumann << std::endl;
            std::cout << "p->CORNER = " << p->CORNER << std::endl;
            throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
        }
    }
    std::cout << "p->U_BUFFER finished" << std::endl;
};

#define derivatives_debug
struct derivatives
{
    // public:
    uomap_P_Tddd P_dxdt_correct, P_gradPhi, P_gradPhi_tangential, P_phin_vector, P_dxdt, P_dxdt_mod, P_laplacian, P_U_dot_gradgrad_U;
    uomap_P_d P_DphiDt, P_kappa, P_pressure, P_aphiat, P_aphiant;
    double mean_surface_height_from_zero = 0;
    ~derivatives(){
        // std::cout << "derivatives 破棄" << std::endl;
    };
    derivatives(const Network &net, bool adjust_dirichlet = false)
    {
        std::unordered_set<networkPoint *> Points = net.getPoints();
        P_dxdt_correct.reserve(Points.size());
        P_gradPhi.reserve(Points.size());
        P_gradPhi_tangential.reserve(Points.size());
        P_phin_vector.reserve(Points.size());
        P_dxdt.reserve(Points.size());
        P_dxdt_mod.reserve(Points.size());
        P_laplacian.reserve(Points.size());
        P_U_dot_gradgrad_U.reserve(Points.size());
        P_DphiDt.reserve(Points.size());
        P_kappa.reserve(Points.size());
        P_pressure.reserve(Points.size());
        P_aphiat.reserve(Points.size());
        P_aphiant.reserve(Points.size());
        // 	this->set(Points);
        // };
        // void set(const std::unordered_set<networkPoint*>& Points)
        // {
#ifdef derivatives_debug
        std::cout << Red << "initialize for parallelization" << reset << std::endl;
#endif
        //! initialize for parallelization
        for (const auto &p : Points)
        {
            this->P_gradPhi[p] = {0, 0, 0};
            this->P_DphiDt[p] = 0.;
        }
        this->P_dxdt_correct = this->P_laplacian = this->P_U_dot_gradgrad_U = this->P_dxdt_mod = this->P_dxdt = this->P_phin_vector = this->P_gradPhi_tangential = this->P_gradPhi;
        this->P_kappa = this->P_pressure = this->P_aphiant = this->P_aphiat = this->P_DphiDt;

        int c = 0;
        for (const auto &p : Points)
        {
            if (p->Dirichlet)
            {
                mean_surface_height_from_zero += p->height();
                c++;
            }
        }

        int depthlimit = 9;
        auto Faces = net.getFaces();

        mean_surface_height_from_zero /= c;
        auto pointsbegin = Points.begin();

        /* ------------------------------------------------------ */
        // //@ pullしたい構造物と，接触を感知した流体格子は角に最大の値を与える．
        std::unordered_set<Network *> nets;
        for (const auto &p : Points)
        {
            auto contactFs = p->getContactFaces();
            if (!contactFs.empty())
                for (const auto &cf : contactFs)
                {
                    auto contact_net = cf->getNetwork();
                    if (p->getNetwork() != contact_net)
                        nets.emplace(contact_net);
                }
        }

#ifdef derivatives_debug
        std::cout << "流速の計算" << std::endl;

#endif
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (const auto &p : Points)
#ifdef _OPENMP
#pragma omp single nowait
#endif
        {
            // std::cout << "曲率の計算" << std::endl;
            if (!isFinite(p->phiphin))
            {
                std::cout << "p->phiphinはfiniteではない！！" << std::endl;
                std::cout << "p->phiphin = " << p->phiphin << std::endl;
                if (p->Neumann)
                    std::cout << "p->Neumann" << std::endl;
                if (p->CORNER)
                    std::cout << "p->CORNER" << std::endl;
                if (p->Dirichlet)
                    std::cout << "p->CORNER" << std::endl;
                throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
            }
            V_netPp ps = Flatten(BFS(p, 3));
            auto interpNormals = InterpolationVectorRBF(ToVector(extX(ps)), ToVector(extNormals(ps)), p->getX());
            p->kappa_BEM = interpNormals.div(p->getX()) / 2.; //中心方向法線ベクトルの場合，マイナスをつける．
            //! https://en.wikipedia.org/wiki/Mean_curvature
            auto grad_N_S_full = gradPhi(p);
            p->U_BEM = std::get<2>(grad_N_S_full);
            p->U_tangential_BEM = std::get<1>(grad_N_S_full);
            p->U_normal_BEM = std::get<0>(grad_N_S_full);
            /* -------------------- おおよそのアップデート流速 ------------------- */
            //@ U_update_BEM は first guess
            // 2022/06/17
            // if (p->CORNER)
            // {
            //     Tddd tang = Normalize(Cross(p->getNormalDirichlet_BEM(), p->getNormalNeumann_BEM()));
            //     p->U_update_BEM = p->U_BEM - Dot(p->U_BEM, tang) * tang;
            // }
            // else

            if (p->Neumann)
                p->U_update_BEM = uNeumann(p);
            else
                p->U_update_BEM = p->U_BEM;
        }

        /* ------------------------------------------------------ */
        //@ この後U_update_BEMをclingなどを使って修正する
        //ここの．BUFFER：Ω(t+δt)はルンゲクッタが見積もる時刻の表面と一致しているか？
        // for (const auto &p : Points)
        // 	p->X_BUFFER = p->getXtuple() + p->U_update_BEM * dt;

        // for (const auto &p : Points)
        // 	p->X_BUFFER = p->RK_X.getXinit() + p->U_update_BEM * p->RK_X.getdt();

        // 初期から考えた流速を与える必要があるのでは？
        //ルンゲクッタに従って次のΩ(t+δt)を予測する
        for (const auto &p : Points)
        {
            p->X_BUFFER = p->RK_X.getX(p->U_update_BEM);
            if (!isFinite(p->X_BUFFER))
                throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
        }
#ifdef derivatives_debug
        std::cout << "ラプラシアンを計算" << std::endl;
#endif
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (const auto &p : Points)
#ifdef _OPENMP
#pragma omp single nowait
#endif
        {
            auto ps = Flatten(BFS(p, 2, {p->getNetwork()}));
            auto intp = InterpolationVectorRBF(obj3D::extractX(ps), extVelocities(ps), p->getX());
            auto tmp = intp.laplacian(p->getX());
            p->laplacian_U_BEM = {tmp[0], tmp[1], tmp[2]};
        }
#ifdef derivatives_debug
        std::cout << "DphiDtを計算" << std::endl;
#endif
        // double gamma = 72.75 * 1E-3;	  //[N/m] 水20度
        // double gravity = _GRAVITY_;		  //[m/s2]
        // double density = _WATER_DENSITY_; //[kg/m3]
        // double nu = 0.01005 / density;

        for (auto &[p, v] : this->P_kappa)
            v = p->kappa_BEM;

        for (auto &[p, v] : this->P_gradPhi)
            v = p->U_BEM;

        for (auto &[p, v] : this->P_laplacian)
            v = p->laplacian_U_BEM;

        for (auto &[p, v] : this->P_dxdt_mod)
            v = p->U_update_BEM;

        for (auto &[p, v] : this->P_gradPhi_tangential)
            v = p->U_tangential_BEM;

        for (auto &[p, v] : this->P_phin_vector)
            v = p->U_normal_BEM;

        /*
            この方法は，ノイマン境界条件のclingにおいて，とても自然に無理なく応用できる．
            ディリクレ境界条件に関しては，計算後に補間によってリグリッドしてもいいかもしれない．
        */
        /* ------------------------------------------------------ */
        /*
            X_BUFFERには，U_update_BEMで単純に予測した節点が保存されている．
            予測したディリクレ面は正しいと考える．
            予測したノイマン面よりも，実際に物体を移動させて作った面の方が正しい．
            そこで，ノイマン面と角点に関しては，物体を移動させて作った面に貼り付ける．
        */
        for (const auto &p : Points)
            p->U_BUFFER = p->U_BUFFER_BUFFER = {0., 0., 0.};

        calculateVectorToSurfaceInBuffer(net, adjust_dirichlet);
        calculateVectorFromBufferToContactFaces(net);
        // /* ------------------------------------------------------ */

        for (const auto &p : Points)
        {
            p->U_update_BEM += p->U_BUFFER / p->RK_X.getdt();
            P_dxdt_correct[p] = p->U_BUFFER / p->RK_X.getdt();
            if (!isFinite(p->U_update_BEM, 1E+10) || !isFinite(p->U_BUFFER, 1E+10))
            {
                std::cout << "p->RK_X.getdt() = " << p->RK_X.getdt() << std::endl;
                std::cout << "p->U_update_BEM = " << p->U_update_BEM << std::endl;
                std::cout << "p->U_BUFFER = " << p->U_BUFFER << std::endl;
                std::cout << "p->Dirichlet = " << p->Dirichlet << std::endl;
                std::cout << "p->Neumann = " << p->Neumann << std::endl;
                std::cout << "p->CORNER = " << p->CORNER << std::endl;
                throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
            }
        }

#ifdef derivatives_debug
        std::cout << " DONE" << std::endl;
#endif

        /* ------------------------------------------------------ */

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (const auto &p : Points)
#ifdef _OPENMP
#pragma omp single nowait
#endif
        {
            this->P_dxdt[p] = p->U_update_BEM; //流速
            // //!この場合マイナスでないと，上の部分が半たんする
            auto DphiDt = p->DphiDt(p->U_update_BEM, 0.);
            //!!ノイマンの場合はこれでDphiDtは計算できませんよ！！！
            this->P_DphiDt[p] = DphiDt; // update用
            if (p->Neumann || p->CORNER)
            {
                auto n = p->getNormalNeumann_BEM();
                T6d VW = velocity_from_Neumann_surface(p);
                Tddd angular_velocity = {std::get<3>(VW), std::get<4>(VW), std::get<5>(VW)};
                auto Q = Quaternion();
                auto dQdt = Q.d_dt(angular_velocity);
                std::get<1>(p->phiphin_t) = accel_normal_from_Neumann_surface(p) - Dot(n, Dot(p->U_BEM, grad_U_LinearElement(p)));
                std::get<1>(p->phiphin_t) += Dot(n, Dot(velocity_normal_from_Neumann_surface(p) - p->U_BEM, dQdt.Rv()));
                /*
                ∇U=∇∇f={{fxx, fyx, fzx},
                        {fxy, fyy, fzy},
                        {fxz, fyz, fzz}}
                なので，∇∇f=∇∇f^T
                */
            }
            if (p->Dirichlet || p->CORNER)
                std::get<0>(p->phiphin_t) = p->DphiDt(0.) - Dot(p->U_BEM, p->U_BEM); //!!ノイマンの場合はこれでDphiDtは計算できませんよ！！！

            this->P_aphiat[p] = std::get<0>(p->phiphin_t);  //!!ノイマンの場合はこれでDphiDtは計算できませんよ！！！
            this->P_aphiant[p] = std::get<1>(p->phiphin_t); //!!ノイマンの場合はこれでDphiDtは計算できませんよ！！！

            this->P_U_dot_gradgrad_U[p] = Dot(p->U_BEM, grad_U_LinearElement(p));
            this->P_pressure[p] = p->pressure_BEM;

            // 10000. * (-1 / 2. * Vphi_Vphi - gravity * (std::get<2>(p->getXtuple())) - aphiat);
            //ここで圧力の項が抜けているが，これは全く流速に関係ないことに気づく．
            //なぜなら，表面上のどこでも同じだけ増加に寄与する大気圧は，
            // grad phiの計算によって，定数のため相殺されるからだ．
        }

#ifdef derivatives_debug
        std::cout << "derivatives終了" << std::endl;
#endif
    }
};

#endif