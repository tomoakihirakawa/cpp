#include "InterpolationRBF.hpp"
#include "Network.hpp"

// int main()
// {
//     NetworkObj obj("../../obj/pumpkin.obj");
//     obj.scale({1., 1., 1.});
//     mk_vtu("./vtu/obj_org.vtu", obj.getFaces());
//     double dx = .01;
//     auto tmp = new Network;
//     for (const auto &X : Flatten(InterpolateFacesLinear(obj.getFaces(), dx, {0})))
//         new networkPoint(tmp, tmp, X);
//     mk_vtu("./vtu/obj.vtu", {tmp->getPoints()});
// }
//! * ------------------------------------------------------ */
/**
 * BFSの例
 */
// int main()
// {
//     NetworkObj obj("../../obj/bunny.obj");
//     mk_vtu("./vtu/obj.vtu", obj.getFaces());

//     auto i = 0;
//     for (const auto &v_points : BFS(obj.getPoints()[0], 27))
//         mk_vtu("./vtu/bsf_points" + std::to_string(i++) + ".vtu", {v_points});

//     auto j = 0;
//     for (const auto &v_points : BFS(obj.getFaces()[0], 27))
//         mk_vtu("./vtu/bsf_faces" + std::to_string(j++) + ".vtu", {v_points});
// }
//! * ------------------------------------------------------ */
// int main()
// {
//     NetworkObj obj("../../obj/pumpkin.obj");
//     mk_vtu("./vtu/obj.vtu", obj.getFaces());

//     using Var = std::variant<std::string, int, std::map<netPp, V_d>>;
//     using VV_Var = std::vector<std::vector<Var>>;

//     std::map<netPp, V_d> m;
//     VV_Var data;
//     auto i = 0;
//     for (const auto &v_points : BFS(obj.getPoints()[0], 27))
//     {
//         for (const auto &p : v_points)
//             m[p] = {(double)i};
//         i++;
//     }
//     mk_vtu("./vtu/bsf_points_colored_by_depth.vtu", {obj.getPoints()}, {{"depth", 1, m}});

//     //! not yes impolemented
//     // std::map<netFp, V_d> n;
//     // auto j = 0;
//     // for (const auto &v_faces : BFS(obj.getFaces()[0], 27))
//     //     for (const auto &p : v_faces)
//     //         n[p] = {(double)j};
//     // mk_vtu("./vtu/bsf_faces_colored_by_depth.vtu", obj.getFaces(), {{"depth", 1, n}});
// }
/* ------------------------------------------------------ */
/*            極座標パラメトリック空間にRBFでマップする例        */
/* ------------------------------------------------------ */
// int main()
// {
//     NetworkObj obj("../../obj/bunny.obj");
//     mk_vtu("./vtu/obj.vtu", obj.getFaces());

//     Network net;
//     for (auto l = 0; l < 50; l++)
//     {
//         VV_d points = {};
//         VV_d xyz = {};
//         VV_d param = {};
//         // for (const auto &tup : obj.getPoints()[l]->getNeighbors_Depth2_OnPolarAsTuple2();
//         // {
//         //     auto [t0, t1, X] = tup;
//         //     param.emplace_back(V_d{t0, t1});
//         //     xyz.emplace_back(X);
//         // }
//         networkPoint *p;
//         for (const auto &tup : obj.getPoints()[l]->getNeighbors_Depth2_OnPolarAsTuple())
//         {
//             auto [t0, t1, X, p] = tup;
//             param.emplace_back(V_d{t0, t1});
//             xyz.emplace_back(X);
//         }

//         param.emplace_back(V_d{0, 0});
//         xyz.emplace_back(obj.getPoints()[l]->getX());
//         auto intp = InterpolationVectorRBF(param, xyz);
//         int n = 50;
//         for (auto i = 0; i < n; i++)
//             for (auto j = 0; j < n; j++)
//             {
//                 double theta = 2. * M_PI * i / (n - 1.);
//                 double r = j / (n - 1.);
//                 points.emplace_back(intp({r * cos(theta), r * sin(theta)}));
//             }
//         mk_vtu("./vtu/polared_points" + std::to_string(l) + ".vtu", {points});
//     }
// }
/* ------------------------------------------------------ */
/*                         粒子化の例                       */
/* ------------------------------------------------------ */
int main()
{
    auto net = new NetworkObj("../../obj/bunny.obj");
    mk_vtu("./vtu/obj.vtu", net->getFaces());
    for (const auto &f : net->getFaces())
    {
        auto h = sqrt(f->getArea()) / 5. /*粒子間隔*/;
        for (const auto &X : particlize(f, 0.001 / 2))
            new networkPoint(net, net, X);
    }
    mk_vtu("./vtu/obj_particlized.vtu", {net->getPoints()});
    delete net;
}
