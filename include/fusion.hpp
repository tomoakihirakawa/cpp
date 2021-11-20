#include "rootFinding.hpp"

double norm_f(const T4d &Q, const Tddd &A, const Tddd &M, const Tddd &G0, const Tddd &M0)
{
    //ゼロにしたい関数
    auto [Ax, Ay, Az] = A;     // measured
    auto [Mx, My, Mz] = M;     // measured
    auto [gx0, gy0, gz0] = G0; //
    auto [mx0, my0, mz0] = M0; //
    auto [a, b, c, d] = Q;     // 計算しているクォータニオン
    double a2 = a * a, b2 = b * b, c2 = c * c, d2 = d * d;
    return (pow(Ax - a2 * gx0 - b2 * gx0 + c2 * gx0 + d2 * gx0 - 2 * b * c * gy0 - 2 * a * d * gy0 + 2 * a * c * gz0 - 2 * b * d * gz0, 2) + pow(Az - 2 * a * c * gx0 - 2 * b * d * gx0 + 2 * a * b * gy0 - 2 * c * d * gy0 - a2 * gz0 + b2 * gz0 + c2 * gz0 - d2 * gz0, 2) + pow(Ay + 2 * a * d * gx0 - a2 * gy0 + b2 * gy0 - c2 * gy0 + d2 * gy0 - 2 * c * d * gz0 - 2 * b * (c * gx0 + a * gz0), 2) + pow(Mx - a2 * mx0 - b2 * mx0 + c2 * mx0 + d2 * mx0 - 2 * b * c * my0 - 2 * a * d * my0 + 2 * a * c * mz0 - 2 * b * d * mz0, 2) + pow(-2 * a * c * mx0 - 2 * b * d * mx0 + 2 * a * b * my0 - 2 * c * d * my0 + Mz - a2 * mz0 + b2 * mz0 + c2 * mz0 - d2 * mz0, 2) + pow(2 * a * d * mx0 + My - a2 * my0 + b2 * my0 - c2 * my0 + d2 * my0 - 2 * c * d * mz0 - 2 * b * (c * mx0 + a * mz0), 2)) / 4.;
};

T4d DDq_norm_f(const T4d &Q, const Tddd &A, const Tddd &M, const Tddd &G0, const Tddd &M0)
{
    //ゼロにしたい関数のqに関する微分
    auto [Ax, Ay, Az] = A;     // measured
    auto [Mx, My, Mz] = M;     // measured
    auto [a, b, c, d] = Q;     // 計算しているクォータニオン
    auto [gx0, gy0, gz0] = G0; //
    auto [mx0, my0, mz0] = M0; //
    double a2 = a * a, b2 = b * b, c2 = c * c, d2 = d * d;
    double a2b2c2d2 = a2 + b2 + c2 + d2;
    double normGM2 = gx0 * gx0 + gy0 * gy0 + gz0 * gz0 + mx0 * mx0 + my0 * my0 + mz0 * mz0;
    return {Az * (-(c * gx0) + b * gy0 - a * gz0) + Ay * (d * gx0 - a * gy0 - b * gz0) + Ax * (-(a * gx0) - d * gy0 + c * gz0) + Mz * (-(c * mx0) + b * my0 - a * mz0) + My * (d * mx0 - a * my0 - b * mz0) + Mx * (-(a * mx0) - d * my0 + c * mz0) + a * a2b2c2d2 * normGM2, Ay * (-(c * gx0) + b * gy0 - a * gz0) + Az * (-(d * gx0) + a * gy0 + b * gz0) + Ax * (-(b * gx0) - c * gy0 - d * gz0) + My * (-(c * mx0) + b * my0 - a * mz0) + Mz * (-(d * mx0) + a * my0 + b * mz0) + Mx * (-(b * mx0) - c * my0 - d * mz0) + a2b2c2d2 * b * normGM2, Ax * (c * gx0 - b * gy0 + a * gz0) + Az * (-(a * gx0) - d * gy0 + c * gz0) + Ay * (-(b * gx0) - c * gy0 - d * gz0) + Mx * (c * mx0 - b * my0 + a * mz0) + Mz * (-(a * mx0) - d * my0 + c * mz0) + My * (-(b * mx0) - c * my0 - d * mz0) + a2b2c2d2 * c * normGM2, Ax * (d * gx0 - a * gy0 - b * gz0) + Ay * (a * gx0 + d * gy0 - c * gz0) + Az * (-(b * gx0) - c * gy0 - d * gz0) + Mx * (d * mx0 - a * my0 - b * mz0) + My * (a * mx0 + d * my0 - c * mz0) + Mz * (-(b * mx0) - c * my0 - d * mz0) + a2b2c2d2 * d * normGM2};
};

T4T4d D2D2q_norm_f(const T4d &Q, const Tddd &A, const Tddd &M, const Tddd &G0, const Tddd &M0)
{
    //ゼロにしたい関数の微分の微分
    auto [Ax, Ay, Az] = A;     // measured
    auto [Mx, My, Mz] = M;     // measured
    auto [a, b, c, d] = Q;     // 計算しているクォータニオン
    auto [gx0, gy0, gz0] = G0; //
    auto [mx0, my0, mz0] = M0; //
    double a2 = a * a, b2 = b * b, c2 = c * c, d2 = d * d;
    double a2b2c2d2 = a2 + b2 + c2 + d2;
    double normGM2 = gx0 * gx0 + gy0 * gy0 + gz0 * gz0 + mx0 * mx0 + my0 * my0 + mz0 * mz0;
    return {{-(Ax * gx0) - Ay * gy0 - Az * gz0 - Mx * mx0 - My * my0 - Mz * mz0 + (3 * a2 + b2 + c2 + d2) * normGM2,
             Az * gy0 - Ay * gz0 + my0 * Mz - My * mz0 + 2 * a * b * normGM2,
             -(Az * gx0) + Ax * gz0 - mx0 * Mz + Mx * mz0 + 2 * a * c * normGM2, Ay * gx0 - Ax * gy0 + mx0 * My - Mx * my0 + 2 * a * d * normGM2},
            {Az * gy0 - Ay * gz0 + my0 * Mz - My * mz0 + 2 * a * b * normGM2,
             -(Ax * gx0) + Ay * gy0 + Az * gz0 - Mx * mx0 + My * my0 + Mz * mz0 + (a2 + 3 * b2 + c2 + d2) * normGM2,
             -(Ay * gx0) - Ax * gy0 - mx0 * My - Mx * my0 + 2 * b * c * normGM2, -(Az * gx0) - Ax * gz0 - mx0 * Mz - Mx * mz0 + 2 * b * d * normGM2},
            {-(Az * gx0) + Ax * gz0 - mx0 * Mz + Mx * mz0 + 2 * a * c * normGM2,
             -(Ay * gx0) - Ax * gy0 - mx0 * My - Mx * my0 + 2 * b * c * normGM2,
             Ax * gx0 - Ay * gy0 + Az * gz0 + Mx * mx0 - My * my0 + Mz * mz0 + (a2 + b2 + 3 * c2 + d2) * normGM2, -(Az * gy0) - Ay * gz0 - my0 * Mz - My * mz0 + 2 * c * d * normGM2},
            {Ay * gx0 - Ax * gy0 + mx0 * My - Mx * my0 + 2 * a * d * normGM2,
             -(Az * gx0) - Ax * gz0 - mx0 * Mz - Mx * mz0 + 2 * b * d * normGM2, -(Az * gy0) - Ay * gz0 - my0 * Mz - My * mz0 + 2 * c * d * normGM2, Ax * gx0 + Ay * gy0 - Az * gz0 + Mx * mx0 + My * my0 - Mz * mz0 + (a2 + b2 + c2 + 3 * d2) * normGM2}};
};

struct Fusion
{
    std::tuple<double, Quaternion, Tddd, Tddd, Tddd, Tddd> time_Q0;
    std::tuple<double, Quaternion, Tddd, Tddd, Tddd, Tddd> time_Q1;
    std::tuple<double, Quaternion, Tddd, Tddd, Tddd, Tddd> time_Q2;
    Tddd G0;
    Tddd M0;
    NewtonRaphson<T4d> NR;
    Fusion(const Tddd &G0_IN, const Tddd &M0_IN)
        : G0(G0_IN), M0(M0_IN), NR({1., 0., 0., 0.}),
          time_Q0({0., Quaternion(), {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}}),
          time_Q1({0., Quaternion(), {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}}),
          time_Q2({0., Quaternion(), {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}}){};
    // gxS, gx, u2, u3, u4, u5はセンサーが実際に計測した値 重力と磁場
    void setG(const Tddd &G0_IN) { this->G0 = G0_IN; };
    void setM(const Tddd &M0_IN) { this->M0 = M0_IN; };
    Quaternion solveForQuaternion(const Tddd &A, const Tddd &M, const double time = 0.)
    {
        // Magnetについては，Normalizeした方がいいかもしれない．
        //  NRは前回の結果を引き継いでいる
        for (auto i = 0; i < 50; ++i)
        {
            NR.update(DDq_norm_f(NR.X, A, M, G0, M0), D2D2q_norm_f(NR.X, A, M, G0, M0)); /*initial X is updated*/
            NR.X = Normalize(NR.X);
            if (Norm(NR.dX) < 1E-8)
                break;
        }
        time_Q1 = time_Q0;
        time_Q2 = time_Q1;
        std::get<0>(time_Q2) = time;
        std::get<1>(time_Q2).set(NR.X);
        std::get<2>(time_Q2) = A;
        std::get<3>(time_Q2) = M;
        std::get<4>(time_Q2) = {0., 0., 0.};
        std::get<5>(time_Q2) = A - std::get<1>(time_Q2).Rs(G0);
        time_Q0 = time_Q2;
        return std::get<1>(time_Q2);
    };
    Quaternion solveForQuaternion(const Tddd &A, const Tddd &M, const Tddd &W, const double time = 0.)
    {
        // NRは前回の結果を引き継いでいる
        for (auto i = 0; i < 50; ++i)
        {
            NR.update(DDq_norm_f(NR.X, A, M, G0, M0), D2D2q_norm_f(NR.X, A, M, G0, M0)); /*initial X is updated*/
            NR.X = Normalize(NR.X);
            if (Norm(NR.dX) < 1E-8)
                break;
        }
        time_Q1 = time_Q0;
        time_Q2 = time_Q1;
        std::get<0>(time_Q2) = time;
        std::get<1>(time_Q2).set(NR.X);
        std::get<2>(time_Q2) = A;
        std::get<3>(time_Q2) = M;
        std::get<4>(time_Q2) = W;
        std::get<5>(time_Q2) = A - std::get<1>(time_Q2).Rs(G0);
        time_Q0 = time_Q2;
        return std::get<1>(time_Q2);
    };
};