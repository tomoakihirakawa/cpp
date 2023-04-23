#ifndef networkTetra_H
#define networkTetra_H

#include "Network.hpp"

inline networkTetra::networkTetra(Network *network_IN,
                                  const T_4P &points,
                                  const T_6L &lines,
                                  const T_4F &faces)
    : Tetrahedron(ToX(points)),
      network(network_IN),
      Points(points),
      Lines(lines),
      Faces(faces) {
   this->network->add(this);
   // 面に四面体を格納
   std::ranges::for_each(faces, [&](const auto &f) { f->Tetras = {this, std::get<0>(f->Tetras)}; });
};

#endif