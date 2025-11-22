#include "basic.hpp"
#include "vtkWriter.hpp"
/* ------------------------------------------------------ */

#define check_tetra_from_Tetrahedron
#if defined(check_tetra_from_Tetrahedron)
int main(int argc, char **argv)
{
    std::shared_ptr<geometry::CoordinateBounds> cube(new geometry::CoordinateBounds(T2Tddd{{0., 0., 0.}, {1., 1., 1.}}));
    std::shared_ptr<Tddd> x0(new Tddd({0, 0, 0}));
    std::shared_ptr<Tddd> x1(new Tddd({1, 0, 0}));
    std::shared_ptr<Tddd> x2(new Tddd({1, 1, 0}));
    std::shared_ptr<Tddd> x3(new Tddd({0, 1, 0}));
    std::shared_ptr<Tddd> x4(new Tddd({0, 0, 1}));
    std::shared_ptr<Tddd> x5(new Tddd({1, 0, 1}));
    std::shared_ptr<Tddd> x6(new Tddd({1, 1, 1}));
    std::shared_ptr<Tddd> x7(new Tddd({0, 1, 1}));
    vtkPolygonWriter<std::shared_ptr<Tddd>> vtp;

    vtp.add({x0, x1, x2, x3, x4, x5, x6, x7});
    // using id = std::tuple<std::shared_ptr<Tddd>, int>;
    vtp.addPolygon({x0, x1, x2, x3});
    vtp.addPolygon({{x0, x1, x2, x3},
                    {x4, x5, x6, x7},
                    {x0, x1, x5, x4},
                    {x2, x3, x7, x6},
                    {x0, x4, x7, x3},
                    {x1, x2, x6, x5}});
    vtp.addPolygon({x0, x4, x7, x3});
    vtp.addLine({x1, x4, x2, x5, x6});
    vtp.addPolygon({x1, x6, x5});
    std::ofstream ofs("./output2.vtp");
    ofs << "<?xml version=\"1.0\"?>\n";
    vtp.write(ofs);
    ofs.close();
    //
    std::stringstream ss;
    vtp.write(ss);
    std::cout << ss.str() << std::endl;
};
#elif defined(check_tetra)
int main(int argc, char **argv)
{
    std::shared_ptr<XMLElement> VTKFile(new XMLElement("VTKFile", {{"type", "PolyData"}, {"version", "0.1"}, {"byte_order", "LittleEndian"}}));
    std::shared_ptr<XMLElement> PolyData(new XMLElement("PolyData"));
    std::shared_ptr<XMLElement> Piece(new XMLElement("Piece", {{"NumberOfPoints", "4"},
                                                               {"NumberOfVerts", "0"},
                                                               {"NumberOfLines", "0"},
                                                               {"NumberOfStrips", "0"},
                                                               {"NumberOfPolys", "4"}}));
    VTKFile->add(PolyData);
    PolyData->add(Piece);
    //
    std::shared_ptr<XMLElement> Points(new XMLElement("Points"));
    std::shared_ptr<XMLElement> PointData(new XMLElement("PointData", {{"Scalars", "my_scalars"}}));
    std::shared_ptr<XMLElement> CellData(new XMLElement("CellData", {{"Scalars", "cell_scalars"}, {"Normals", "cell_normals"}}));
    std::shared_ptr<XMLElement> Polys(new XMLElement("Polys"));
    Piece->add(Points);
    Piece->add(PointData);
    Piece->add(CellData);
    Piece->add(Polys);
    /* ------------------------- Points <- 座標 ------------------------- */
    {
        std::shared_ptr<XMLElement> DataArray(new XMLElement("DataArray", {{"type", "Float32"}, {"NumberOfComponents", "3"}, {"format", "ascii"}}));
        // DataArray->content = "0 0 0 1 0 0 1 1 0 0 1 0 0 0 1 1 0 1 1 1 1 0 1 1 ";
        DataArray->writer = [&](std::stringstream &ofs)
        { ofs << "0 0 0 1 0 0 1 1 0 1 1 1"; };
        Points->add(DataArray);
    }
    /* ----------------- PointData <- 点上のスカラー ---------------------- */
    {
        std::shared_ptr<XMLElement> DataArray(new XMLElement("DataArray", {{"type", "Float32"}, {"Name", "my_scalars"}, {"format", "ascii"}}));
        DataArray->writer = [&](std::stringstream &ofs)
        { ofs << "0 1 2 3"; };
        PointData->add(DataArray);
    }
    /* ------------------ CellData <- セル上のスカラー ---------------------- */
    // {
    //     std::shared_ptr<XMLElement> DataArray(new XMLElement("DataArray", {{"type", "Float32"}, {"Name", "cell_scalars"}, {"format", "ascii"}}));
    //     DataArray->writer = [&](std::stringstream &ofs)
    //     { ofs << "0 1 2 3 4 5"; };
    //     CellData->add(DataArray);
    // }
    // {
    //     std::shared_ptr<XMLElement> DataArray(new XMLElement("DataArray", {{"type", "Float32"}, {"Name", "cell_normals"}, {"NumberOfComponents", "3"}, {"format", "ascii"}}));
    //     DataArray->writer = [&](std::stringstream &ofs)
    //     { ofs << "0 0 -1 0 0 1 0 -1 0 0 1 0 -1 0 0 1 0 0"; };
    //     CellData->add(DataArray);
    // }
    /* --------------------- Polys <- 接続関係とオフセット------------------------ */
    {
        std::shared_ptr<XMLElement> DataArray(new XMLElement("DataArray", {{"type", "Int32"}, {"Name", "connectivity"}, {"format", "ascii"}}));
        DataArray->writer = [&](std::stringstream &ofs)
        { ofs << "0 1 2 0 1 3 2 1 3 0 2 3"; };
        Polys->add(DataArray);
    }
    {
        std::shared_ptr<XMLElement> DataArray(new XMLElement("DataArray", {{"type", "Int32"}, {"Name", "offsets"}, {"format", "ascii"}}));
        DataArray->writer = [&](std::stringstream &ofs)
        { ofs << "3 6 9 12"; };
        Polys->add(DataArray);
    }
    /* ------------------------------------------------------ */
    std::ofstream ofs("./output.vtp");
    ofs << "<?xml version=\"1.0\"?>\n";
    VTKFile->write(ofs);
    ofs.close();
};
#elif defined(check_cube)
int main(int argc, char **argv)
{
    std::shared_ptr<XMLElement> VTKFile(new XMLElement("VTKFile", {{"type", "PolyData"}, {"version", "0.1"}, {"byte_order", "LittleEndian"}}));
    std::shared_ptr<XMLElement> PolyData(new XMLElement("PolyData"));
    std::shared_ptr<XMLElement> Piece(new XMLElement("Piece", {{"NumberOfPoints", "8"},
                                                               {"NumberOfVerts", "0"},
                                                               {"NumberOfLines", "0"},
                                                               {"NumberOfStrips", "0"},
                                                               {"NumberOfPolys", "6"}}));
    VTKFile->add(PolyData);
    PolyData->add(Piece);
    //
    std::shared_ptr<XMLElement> Points(new XMLElement("Points"));
    std::shared_ptr<XMLElement> PointData(new XMLElement("PointData", {{"Scalars", "my_scalars"}}));
    std::shared_ptr<XMLElement> CellData(new XMLElement("CellData", {{"Scalars", "cell_scalars"}, {"Normals", "cell_normals"}}));
    std::shared_ptr<XMLElement> Polys(new XMLElement("Polys"));
    Piece->add(Points);
    Piece->add(PointData);
    Piece->add(CellData);
    Piece->add(Polys);
    /* ------------------------- Points <- 座標 ------------------------- */
    {
        std::shared_ptr<XMLElement> DataArray(new XMLElement("DataArray", {{"type", "Float32"}, {"NumberOfComponents", "3"}, {"format", "ascii"}}));
        // DataArray->content = "0 0 0 1 0 0 1 1 0 0 1 0 0 0 1 1 0 1 1 1 1 0 1 1 ";
        DataArray->writer = [&](std::stringstream &ofs)
        { ofs << "0 0 0 1 0 0 1 1 0 0 1 0 0 0 1 1 0 1 1 1 1 0 1 1 "; };
        Points->add(DataArray);
    }
    /* ----------------- PointData <- 点上のスカラー ---------------------- */
    {
        std::shared_ptr<XMLElement> DataArray(new XMLElement("DataArray", {{"type", "Float32"}, {"Name", "my_scalars"}, {"format", "ascii"}}));
        DataArray->writer = [&](std::stringstream &ofs)
        { ofs << "0 1 2 3 4 5 6 7"; };
        PointData->add(DataArray);
    }
    /* ------------------ CellData <- セル上のスカラー ---------------------- */
    {
        std::shared_ptr<XMLElement> DataArray(new XMLElement("DataArray", {{"type", "Float32"}, {"Name", "cell_scalars"}, {"format", "ascii"}}));
        DataArray->writer = [&](std::stringstream &ofs)
        { ofs << "0 1 2 3 4 5"; };
        CellData->add(DataArray);
    }
    {
        std::shared_ptr<XMLElement> DataArray(new XMLElement("DataArray", {{"type", "Float32"}, {"Name", "cell_normals"}, {"NumberOfComponents", "3"}, {"format", "ascii"}}));
        DataArray->writer = [&](std::stringstream &ofs)
        { ofs << "0 0 -1 0 0 1 0 -1 0 0 1 0 -1 0 0 1 0 0"; };
        CellData->add(DataArray);
    }
    /* --------------------- Polys <- 接続関係とオフセット------------------------ */
    {
        std::shared_ptr<XMLElement> DataArray(new XMLElement("DataArray", {{"type", "Int32"}, {"Name", "connectivity"}, {"format", "ascii"}}));
        DataArray->writer = [&](std::stringstream &ofs)
        { ofs << "0 1 2 3 4 5 6 7 0 1 5 4 2 3 7 6 0 4 7 3 1 2 6 5 "; };
        Polys->add(DataArray);
    }
    {
        std::shared_ptr<XMLElement> DataArray(new XMLElement("DataArray", {{"type", "Int32"}, {"Name", "offsets"}, {"format", "ascii"}}));
        DataArray->writer = [&](std::stringstream &ofs)
        { ofs << "4 8 12 16 20 24"; };
        Polys->add(DataArray);
    }
    /* ------------------------------------------------------ */
    std::ofstream ofs("./output.vtp");
    ofs << "<?xml version=\"1.0\"?>\n";
    VTKFile->write(ofs);
    ofs.close();
};
#else
int main(int argc, char **argv)
{
    auto elemHTML = new XMLElement("html");
    auto elemBody = new XMLElement("body");
    //
    auto elemFamily = new XMLElement("family");
    auto elemFather = new XMLElement("father");
    auto elemMother = new XMLElement("mother");
    auto elemSister = new XMLElement("sister");
    auto elemCar = new XMLElement("car");
    //
    elemHTML->add(elemBody);
    elemBody->add(elemFamily);
    //
    elemFamily->add(elemFather);
    elemFather->add(elemCar);
    elemFamily->add(elemMother);
    elemFamily->add(elemSister);
    //
    elemFather->attributes["name"] = "Taro";
    elemMother->attributes["name"] = "Hanako";
    elemSister->attributes["name"] = "Noriko";
    elemCar->attributes["name"] = "BMW";
    //
    elemFather->content = "sick";
    elemMother->content = "fine";
    elemSister->content = "great";
    std::cout << elemHTML->getChildren() << std::endl;

    auto fp = fopen("./output.xml", "wb");
    fprintf(fp, elemHTML->getChildren().c_str());
    fclose(fp);

    delete elemHTML, elemBody, elemFamily, elemFather, elemMother, elemSister, elemCar;
};

#endif