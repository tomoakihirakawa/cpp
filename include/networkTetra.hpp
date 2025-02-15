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

inline networkTetra::~networkTetra() {
   // FacesからのTetraの削除
   for (auto &f : this->Faces) {
      if (f != nullptr) {
         if (f->Tetras[0] == this)
            f->Tetras[0] = nullptr;
         else if (f->Tetras[1] == this)
            f->Tetras[1] = nullptr;
      }
      f = nullptr;
   }

   // LinesからのTetraの削除
   for (auto &l : this->Lines) {
      erase_element(l->Tetras, this);  // l->Tetrasがvector/unordered_setかに応じて動作
      l = nullptr;
   }

   // PointsからのTetraの削除
   for (auto &p : this->Points) {
      erase_element(p->Tetras, this);  // p->Tetrasがvector/unordered_setかに応じて動作
      p = nullptr;
   }

   // ネットワークからのTetra削除
   network->erase_element(this);  // unordered_set対応版erase_element
   network = nullptr;
}
