#pragma once

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
   std::ranges::for_each(faces, [&](const auto &f) {
      if (std::get<0>(f->Tetras) == nullptr)
         std::get<0>(f->Tetras) = this;
      else if (std::get<1>(f->Tetras) == nullptr)
         std::get<1>(f->Tetras) = this;
      else
         throw std::runtime_error("networkTetra::networkTetra: faces are full");
   });
   std::ranges::for_each(lines, [&](const auto &l) { l->Tetras.emplace_back(this); });
   std::ranges::for_each(points, [&](const auto &p) { p->Tetras.emplace_back(this); });
};
