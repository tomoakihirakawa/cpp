#ifndef vtkWriter_H
#define vtkWriter_H
#include "XML.hpp"
#include "basic_arithmetic_vector_operations.hpp"

Tddd ToX(const std::shared_ptr<Tddd> &X) { return *X; };

template <typename T>
struct vtkPolygonWriter : XMLElement {
   /*
   verticesはインデックスと座標ベクトルを含んでいる．
   その組み合わせをidとしている．インデックスは大抵は０．
   ポリゴンやラインを加える際は，このidの列をaddする．
   vtpでは座標ベクトルのインデックスを与えることで，ポリゴンを作るが，idが与えられれば，
   id={インデックス，座標}なので，インデックスもわかるので，自動でインデックスを取得し出力してくれる．

   具体的な使い方については，octreeを見てください．
   */
   /* ------------------------------------------------------ */
   std::shared_ptr<XMLElement> PolyData;
   std::shared_ptr<XMLElement> Piece;
   std::shared_ptr<XMLElement> Points;
   std::shared_ptr<XMLElement> PointData, CellData;
   std::shared_ptr<XMLElement> Polys;
   std::shared_ptr<XMLElement> Lines;
   using id = std::tuple<T, int>;
   std::unordered_map<id, std::tuple<int, Tddd>> verticies;
   std::vector<std::vector<id>> connectivity_lines;
   std::vector<std::tuple<id, id, id>> connectivity3;
   std::vector<std::tuple<id, id, id, id>> connectivity4;
   void clear() {
      verticies.clear();
   };
   void add(const T &v) { std::get<1>(this->verticies[{v, 0}]) = ToX(v); };
   void add(const T &v, const T &u) {
      std::get<1>(this->verticies[{v, 0}]) = ToX(v);
      std::get<1>(this->verticies[{u, 0}]) = ToX(u);
   };
   void add(const T &v, const T &u, const T &w) {
      std::get<1>(this->verticies[{v, 0}]) = ToX(v);
      std::get<1>(this->verticies[{u, 0}]) = ToX(u);
      std::get<1>(this->verticies[{w, 0}]) = ToX(w);
   };
   void add(const T &v, const T &u, const T &w, const T &x) {
      std::get<1>(this->verticies[{v, 0}]) = ToX(v);
      std::get<1>(this->verticies[{u, 0}]) = ToX(u);
      std::get<1>(this->verticies[{w, 0}]) = ToX(w);
      std::get<1>(this->verticies[{x, 0}]) = ToX(x);
   };
   void add(const std::vector<std::vector<T>> &verts) {
      for (const auto &V : verts)
         for (const auto &v : V)
            std::get<1>(this->verticies[{v, 0}]) = ToX(v);
   };
   void add(const std::vector<T> &verts) {
      for (const auto &v : verts)
         std::get<1>(this->verticies[{v, 0}]) = ToX(v);
   };
   void add(const std::unordered_set<T> &verts) {
      for (const auto &v : verts)
         std::get<1>(this->verticies[{v, 0}]) = ToX(v);
   };
   void addLine(const std::vector<T> &conns) {
      std::vector<id> tmp(conns.size());
      for (auto i = 0; i < conns.size(); ++i)
         tmp[i] = {conns[i], 0};
      this->connectivity_lines.emplace_back(tmp);
   };
   void addLine(const T &p0, const T &p1) { this->connectivity_lines.push_back({{p0, 0}, {p1, 0}}); };
   void addLine(const T &p0, const T &p1, const T &p2) { this->connectivity_lines.push_back({{p0, 0}, {p1, 0}, {p2, 0}}); };
   void addLine(const T &p0, const T &p1, const T &p2, const T &p3) { this->connectivity_lines.push_back({{p0, 0}, {p1, 0}, {p2, 0}, {p3, 0}}); };
   void addLine(const std::vector<std::vector<T>> &conns) {
      this->addLine(conns);
   };
   void addPolygon(const std::tuple<T, T, T> &conns) {
      this->connectivity3.push_back({{std::get<0>(conns), 0},
                                     {std::get<1>(conns), 0},
                                     {std::get<2>(conns), 0}});
   };
   void addPolygon(const T &c0, const T &c1, const T &c2) {
      this->connectivity3.push_back({{c0, 0}, {c1, 0}, {c2, 0}});
   };
   void addPolygon(const T &c0, const T &c1, const T &c2, const T &c3) {
      this->connectivity4.push_back({{c0, 0}, {c1, 0}, {c2, 0}, {c3, 0}});
   };
   void addPolygon(const std::tuple<T, T, T, T> &conns) {
      this->connectivity4.push_back({{std::get<0>(conns), 0},
                                     {std::get<1>(conns), 0},
                                     {std::get<2>(conns), 0},
                                     {std::get<3>(conns), 0}});
   };
   void addPolygon(const std::vector<std::tuple<T, T, T>> &conns) {
      for (const auto &c : conns)
         this->addPolygon(c);
   };
   void addPolygon(const std::vector<std::tuple<T, T, T, T>> &conns) {
      for (const auto &c : conns)
         this->addPolygon(c);
   };
   vtkPolygonWriter()
       : XMLElement("VTKFile", {{"type", "PolyData"}, {"version", "0.1"}, {"byte_order", "LittleEndian"}}),  //<- constant setting
         PolyData(new XMLElement("PolyData")),
         Piece(new XMLElement("Piece", {{"NumberOfPoints", "4"}, {"NumberOfVerts", "0"}, {"NumberOfLines", "0"}, {"NumberOfStrips", "0"}, {"NumberOfPolys", "4"}})),
         PointData(new XMLElement("PointData")),
         CellData(new XMLElement("CellData", {{"Scalars", "cell_scalars"}, {"Normals", "cell_normals"}})),
         Points(new XMLElement("Points")),
         Lines(new XMLElement("Lines")),
         Polys(new XMLElement("Polys")) {
      this->XMLElement::add(PolyData);
      PolyData->add(Piece);
      Piece->add(Points);
      Piece->add(PointData);
      Piece->add(CellData);
      Piece->add(Polys);
      Piece->add(Lines);
      /* ------------------------ */
      /* -------------------------------------------------------------------------- */
      {
         std::shared_ptr<XMLElement> DataArray(new XMLElement("DataArray", {{"type", "Float32"}, {"NumberOfComponents", "3"}, {"format", "ascii"}}));
         DataArray->writer = [&](std::stringstream &ofs) {
            int i = 0;
            for (auto &[_, INT_X] : this->verticies) {
               // 書き出すともに番号を保存
               std::get<0>(INT_X) = i++;
               for_each(std::get<1>(INT_X), [&](auto &x) { ofs << std::setprecision(7) << (float)x << " "; });
            }
         };
         this->Points->add(DataArray);
      }
      /* ------------------------------------------------------ */
      /*                        Polygons                        */
      /* ------------------------------------------------------ */
      {  // connectivity
         std::shared_ptr<XMLElement> DataArray(new XMLElement("DataArray", {{"type", "Int32"}, {"Name", "connectivity"}, {"format", "ascii"}}));
         DataArray->writer = [&](std::stringstream &ofs) {
            for (const auto &[ID0, ID1, ID2] : this->connectivity3)
               ofs << std::get<0>(verticies[ID0]) << " " << std::get<0>(verticies[ID1]) << " " << std::get<0>(verticies[ID2]) << " ";
            for (const auto &[ID0, ID1, ID2, ID3] : this->connectivity4)
               ofs << std::get<0>(verticies[ID0]) << " " << std::get<0>(verticies[ID1]) << " " << std::get<0>(verticies[ID2]) << " " << std::get<0>(verticies[ID3]) << " ";
         };
         this->Polys->add(DataArray);
      }
      {  // offsets
         std::shared_ptr<XMLElement> DataArray(new XMLElement("DataArray", {{"type", "Int32"}, {"Name", "offsets"}, {"format", "ascii"}}));
         DataArray->writer = [&](std::stringstream &ofs) {
            int j = 0;
            for (auto i = 0; i < this->connectivity3.size(); ++i)
               ofs << (j += 3) << " ";
            for (auto i = 0; i < this->connectivity4.size(); ++i)
               ofs << (j += 4) << " ";
         };
         this->Polys->add(DataArray);
      }
      /* ------------------------------------------------------ */
      /*                          Lines                         */
      /* ------------------------------------------------------ */
      {  // connectivity
         std::shared_ptr<XMLElement> DataArray(new XMLElement("DataArray", {{"type", "Int32"}, {"Name", "connectivity"}, {"format", "ascii"}}));
         DataArray->writer = [&](std::stringstream &ofs) {
            for (const auto &verts : this->connectivity_lines)
               for (const auto &ID : verts)
                  ofs << std::get<0>(verticies[ID]) << " ";
         };
         this->Lines->add(DataArray);
      }
      {  // offsets
         std::shared_ptr<XMLElement> DataArray(new XMLElement("DataArray", {{"type", "Int32"}, {"Name", "offsets"}, {"format", "ascii"}}));
         DataArray->writer = [&](std::stringstream &ofs) {
            int j = 0;
            for (auto i = 0; i < this->connectivity_lines.size(); ++i)
               ofs << (j += this->connectivity_lines[i].size()) << " ";
         };
         this->Lines->add(DataArray);
      }
   };
   /* ------------------------------------------------------ */
   std::unordered_map<std::string, std::unordered_map<T, Tddd>> data_Tddd;
   std::unordered_map<std::string, std::unordered_map<T, double>> data_double;
   // std::unordered_map<std::string, std::unordered_map<T, int>> data_int;
   void addPointData(const std::string &name, const std::unordered_map<T, Tddd> &V) {
      auto &data = data_Tddd[name];
      data.insert(V.cbegin(), V.cend());
      /* ------------------------------------------------------ */
      for (const auto &elem : this->PointData->elements_shared)
         if (elem->attributes["Name"] == name)
            return;
      /* ------------------------------------------------------ */
      std::shared_ptr<XMLElement> DataArray(new XMLElement("DataArray", {{"type", "Float32"}, {"Name", name}, {"format", "ascii"}, {"NumberOfComponents", "3"}}));
      this->PointData->add(DataArray);
      DataArray->writer = [&](std::stringstream &ofs) {
         auto it = data.begin();
         for (const auto &[ID, value] : this->verticies) {
            it = data.find(std::get<0>(ID));
            if (it != data.end()) {
               for_each(it->second,
                        [&](const auto &x) {
                           if (isFinite(x))
                              ofs << std::setprecision(6) << (Between(x, {-1E-10, 1E-10}) ? 0 : (float)x) << " ";
                           else
                              ofs << std::setprecision(6) << "NaN ";
                        });
            } else
               ofs << "Nan Nan Nan ";
         }
      };
   };
   void addPointData(const std::string &name, const std::unordered_map<T, double> &V) {
      auto &data = data_double[name];
      data.insert(V.cbegin(), V.cend());
      /* ------------------------------------------------------ */
      for (const auto &elem : this->PointData->elements_shared)
         if (elem->attributes["Name"] == name)
            return;
      /* ------------------------------------------------------ */
      std::shared_ptr<XMLElement> DataArray(new XMLElement("DataArray", {{"type", "Float32"}, {"Name", name}, {"format", "ascii"}}));
      this->PointData->add(DataArray);
      DataArray->writer = [&](std::stringstream &ofs) {
         auto it = data.begin();
         for (const auto &[ID, value] : this->verticies) {
            it = data.find(std::get<0>(ID));
            if (it != data.end()) {
               if (isFinite(it->second))
                  ofs << std::setprecision(6) << (Between(it->second, {-1E-10, 1E-10}) ? 0 : (float)it->second) << " ";
               else
                  ofs << std::setprecision(6) << "NaN ";
            } else
               ofs << "Nan ";
         }
      };
   };
   template <typename V>
   V &write(V &ofs) const {
      ofs << "<?xml version=\"1.0\"?>\n";
      Piece->attributes["NumberOfPoints"] = std::to_string(this->verticies.size());
      Piece->attributes["NumberOfPolys"] = std::to_string(this->connectivity3.size() + this->connectivity4.size());
      Piece->attributes["NumberOfLines"] = std::to_string(this->connectivity_lines.size());
      if (this->elements.empty() && this->elements_shared.empty()) {
         if (writer) {
            std::stringstream ss;
            writer(ss);
            ofs << this->getStartTag() << "\n"
                << ss.str() << "\n"
                << this->getEndTag() << "\n";
         } else
            ofs << this->getStartTag() << "\n"
                << this->getEndTag() << "\n";
      } else {
         ofs << this->getStartTag() << "\n";
         for (const auto &elem : this->elements)
            elem->write(ofs);
         for (const auto &elem : this->elements_shared)
            elem->write(ofs);
         ofs << this->getEndTag() << "\n";
      }
      return ofs;
   };
};
/* ------------------------------------------------------ */
/*                         T4Tddd                         */
/* ------------------------------------------------------ */
void vtkPolygonWrite(std::ofstream &ofs, const std::vector<T4Tddd> &V) {
   vtkPolygonWriter<std::shared_ptr<Tddd>> vtp;
   for (const auto &[X0, X1, X2, X3] : V) {
      std::shared_ptr<Tddd> x0(new Tddd(X0));
      vtp.add(x0);
      std::shared_ptr<Tddd> x1(new Tddd(X1));
      vtp.add(x1);
      std::shared_ptr<Tddd> x2(new Tddd(X2));
      vtp.add(x2);
      std::shared_ptr<Tddd> x3(new Tddd(X3));
      vtp.add(x3);
      vtp.addPolygon({x0, x1, x2, x3});
   }
   vtp.write(ofs);
};
// void vtkPolygonWrite(const std::string &name, const std::vector<T4Tddd> &V) {
//    std::ofstream ofs(name);
//    vtkPolygonWrite(ofs, V);
//    ofs.close();
// };
/* ------------------------------------------------------ */
/*                         T3Tddd                         */
/* ------------------------------------------------------ */
void vtkPolygonWrite(std::ofstream &ofs, const std::vector<T3Tddd> &V) {
   vtkPolygonWriter<std::shared_ptr<Tddd>> vtp;
   for (const auto &[X0, X1, X2] : V) {
      std::shared_ptr<Tddd> x0(new Tddd(X0));
      vtp.add(x0);
      std::shared_ptr<Tddd> x1(new Tddd(X1));
      vtp.add(x1);
      std::shared_ptr<Tddd> x2(new Tddd(X2));
      vtp.add(x2);
      vtp.addPolygon({x0, x1, x2});
   }
   vtp.write(ofs);
};
// void vtkPolygonWrite(const std::string &name, const std::vector<T3Tddd> &V) {
//    std::ofstream ofs(name);
//    vtkPolygonWrite(ofs, V);
//    ofs.close();
// };
/* ------------------------------------------------------ */
/*                          Tddd                          */
/* ------------------------------------------------------ */
void vtkPolygonWrite(std::ofstream &ofs, const std::vector<Tddd> &V) {
   vtkPolygonWriter<std::shared_ptr<Tddd>> vtp;
   for (const auto &X : V) {
      std::shared_ptr<Tddd> x(new Tddd(X));
      vtp.add(x);
   }
   vtp.write(ofs);
};
// void vtkPolygonWrite(const std::string &name, const std::vector<Tddd> &V) {
//    std::ofstream ofs(name);
//    vtkPolygonWrite(ofs, V);
//    ofs.close();
// };

/* ------------------------------------------------------ */
/*                         T2Tddd                         */
/* ------------------------------------------------------ */
void vtkPolygonWrite(std::ofstream &ofs, const std::vector<T2Tddd> &V) {
   vtkPolygonWriter<std::shared_ptr<Tddd>> vtp;
   for (const auto &[X0, X1] : V) {
      std::shared_ptr<Tddd> x0(new Tddd(X0));
      vtp.add(x0);
      std::shared_ptr<Tddd> x1(new Tddd(X1));
      vtp.add(x1);
      vtp.addLine({x0, x1});
   }
   vtp.write(ofs);
};
// void vtkPolygonWrite(const std::string &name, const std::vector<T2Tddd> &V) {
//    std::ofstream ofs(name);
//    vtkPolygonWrite(ofs, V);
//    ofs.close();
// };
template <typename T, typename U>
void vtkPolygonWrite(std::ofstream &ofs, const std::vector<std::tuple<T, T>> &V, const std::unordered_map<T, U> &data_double) {
   vtkPolygonWriter<T> vtp;
   for (const auto &[X0, X1] : V) {
      vtp.add(X0, X1);
      vtp.addLine(std::vector<T>{X0, X1});
   }
   vtp.addPointData("data", data_double);
   vtp.write(ofs);
};
/* ------------------------------------------------------ */
void vtkPolygonWrite(std::ofstream &ofs, const Tddd &X) {
   vtkPolygonWriter<std::shared_ptr<Tddd>> vtp;
   std::shared_ptr<Tddd> x(new Tddd(X));
   vtp.add(x);
   vtp.write(ofs);
};

// void vtkPolygonWrite(const std::string &name, const Tddd &V) {
//    std::ofstream ofs(name);
//    vtkPolygonWrite(ofs, V);
//    ofs.close();
// };

template <typename T, typename U>
void vtkPolygonWrite(std::ofstream &ofs, const std::unordered_set<T> &V, const std::unordered_map<T, U> &data_double) {
   vtkPolygonWriter<T> vtp;
   for (const auto &X : V)
      vtp.add(X);
   vtp.addPointData("data", data_double);
   vtp.write(ofs);
};

#if defined(Network_H)
template <typename U>
void vtkPolygonWrite(std::ofstream &ofs, const std::unordered_set<networkFace *> &uo_f, const std::unordered_map<networkPoint *, U> &data_double = {}) {
   vtkPolygonWriter<networkPoint *> vtp;
   for (const auto &f : uo_f) {
      auto abc = f->getPoints();
      vtp.add(std::get<0>(abc), std::get<1>(abc), std::get<2>(abc));
      vtp.addPolygon(abc);
   }
   if (!data_double.empty())
      vtp.addPointData("data", data_double);
   vtp.write(ofs);
};

void vtkPolygonWrite(std::ofstream &ofs, const T_4F &tuple_f) {
   vtkPolygonWriter<networkPoint *> vtp;
   for_each(tuple_f, [&](const auto &f) {
      auto abc = f->getPoints();
      vtp.add(std::get<0>(abc), std::get<1>(abc), std::get<2>(abc));
      vtp.addPolygon(abc);
   });
   vtp.write(ofs);
};

void vtkPolygonWrite(std::ofstream &ofs, const std::unordered_set<networkLine *> &uoL) {
   vtkPolygonWriter<networkPoint *> vtp;
   for (const auto &l : uoL) {
      auto [a, b] = l->getPoints();
      vtp.add(a, b);
      vtp.addLine(a, b);
   }
   vtp.write(ofs);
};

void vtkPolygonWrite(std::ofstream &ofs, const std::unordered_set<networkFace *> &uoF) {
   vtkPolygonWriter<networkPoint *> vtp;
   for (const auto &f : uoF) {
      auto abc = f->getPoints();
      vtp.add(std::get<0>(abc), std::get<1>(abc), std::get<2>(abc));
      vtp.addPolygon(abc);
   };
   vtp.write(ofs);
};

void vtkPolygonWrite(std::ofstream &ofs, const std::unordered_set<networkTetra *> &uoTet) {
   vtkPolygonWriter<networkPoint *> vtp;
   for (const auto &tet : uoTet)
      for_each(tet->Faces, [&](const auto &f) {
         auto abc = f->getPoints();
         vtp.add(std::get<0>(abc), std::get<1>(abc), std::get<2>(abc));
         vtp.addPolygon(abc);
      });
   vtp.write(ofs);
};

#endif
   // template <typename T, typename U>
   // void vtkPolygonWrite(const std::string &name, const std::unordered_set<T> &V, const std::unordered_map<T, U> &data_double) {
   //    std::ofstream ofs(name);
   //    vtkPolygonWrite(ofs, V, data_double);
   //    ofs.close();
   // };

#endif