#pragma once
#include "BEM.hpp"
#include "OutputCommon.hpp"
#include "vtkWriter.hpp"
#include "VPM.hpp"

namespace OutputParaView {

// Element-type-aware VTU output: linear faces → VTK_TRIANGLE (type 5),
// pseudo-quadratic/true-quadratic faces → VTK_QUADRATIC_TRIANGLE (type 22, 6 nodes).
inline void mk_vtu_quadratic(const std::string &filename, const V_netFp &Faces, const VV_VarForOutput &VV_name_comp_mapPVd = {}) {
  try {
#if defined(debug_mk_vtu)
    std::cout << Magenta << filename << colorReset;
    std::cout << "  Faces.size() : " << std::to_string(Faces.size()) << colorReset << " ";
#endif
    struct PointEntry {
      Tddd position;
      networkPoint *ptr;      // non-null for vertex nodes
      networkPoint *mid_pA;   // midpoint endpoint A
      networkPoint *mid_pB;   // midpoint endpoint B
      networkLine *mid_line;  // midpoint edge
      bool is_true_quad_mid;
    };

    std::vector<PointEntry> all_points;
    std::vector<int> cell_sizes;
    std::vector<uint8_t> cell_types;

    for (const auto &f : Faces) {
      auto [p0, p1, p2] = f->getPoints();
      all_points.push_back({p0->getXtuple(), p0, nullptr, nullptr, nullptr, false});
      all_points.push_back({p1->getXtuple(), p1, nullptr, nullptr, nullptr, false});
      all_points.push_back({p2->getXtuple(), p2, nullptr, nullptr, nullptr, false});

      if (f->isPseudoQuadraticElement || f->isTrueQuadraticElement) {
        auto [p0_, l0, p1_, l1, p2_, l2] = f->PLPLPL;
        Tddd mid0, mid1, mid2;

        if (f->isTrueQuadraticElement) {
          mid0 = l0->X_mid;
          mid1 = l1->X_mid;
          mid2 = l2->X_mid;
        } else {
          // Pseudo-quadratic: use DodecaPoints
          // dodecaPoints[0] origin=p0: p0 at (1,0), p1 at (0,1), p2 at (0,0)
          if (f->dodecaPoints[0]) {
            auto &dp = f->dodecaPoints[0];
            mid0 = dp->X(0.5, 0.5) + dp->corner_offset(0.5, 0.5); // l0(p0-p1)
            mid1 = dp->X(0.0, 0.5) + dp->corner_offset(0.0, 0.5); // l1(p1-p2)
            mid2 = dp->X(0.5, 0.0) + dp->corner_offset(0.5, 0.0); // l2(p2-p0)
          } else {
            mid0 = 0.5 * (p0->getXtuple() + p1->getXtuple());
            mid1 = 0.5 * (p1->getXtuple() + p2->getXtuple());
            mid2 = 0.5 * (p2->getXtuple() + p0->getXtuple());
          }
        }

        auto [pa0, pb0] = l0->getPoints();
        auto [pa1, pb1] = l1->getPoints();
        auto [pa2, pb2] = l2->getPoints();
        all_points.push_back({mid0, nullptr, pa0, pb0, l0, f->isTrueQuadraticElement});
        all_points.push_back({mid1, nullptr, pa1, pb1, l1, f->isTrueQuadraticElement});
        all_points.push_back({mid2, nullptr, pa2, pb2, l2, f->isTrueQuadraticElement});

        cell_sizes.push_back(6);
        cell_types.push_back(22); // VTK_QUADRATIC_TRIANGLE
      } else {
        cell_sizes.push_back(3);
        cell_types.push_back(5); // VTK_TRIANGLE
      }
    }

    const int num_points = (int)all_points.size();
    const int num_cells = (int)cell_sizes.size();

    FILE *fp = fopen(filename.c_str(), "wb");
    if (!fp)
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, filename + " can not be opened");

    fprintf(fp, "<?xml version='1.0' encoding='UTF-8'?>\n");
    fprintf(fp, "<VTKFile xmlns='VTK' byte_order='LittleEndian' version='0.1' type='UnstructuredGrid'>\n");
    fprintf(fp, "<UnstructuredGrid>\n");
    fprintf(fp, "<Piece NumberOfCells='%d' NumberOfPoints='%d'>\n", num_cells, num_points);

    // Points
    fprintf(fp, "<Points>\n");
    fprintf(fp, "<DataArray NumberOfComponents='3' type='Float32' Name='Position' format='ascii'>\n");
    for (const auto &pe : all_points)
      fprintf(fp, "%s %s %s ",
              NumtoString(std::get<0>(pe.position)).c_str(),
              NumtoString(std::get<1>(pe.position)).c_str(),
              NumtoString(std::get<2>(pe.position)).c_str());
    fprintf(fp, "\n</DataArray>\n");
    fprintf(fp, "</Points>\n");

    // PointData
    if (!VV_name_comp_mapPVd.empty()) {
      fprintf(fp, "<PointData>\n");
      for (const auto &V_var : VV_name_comp_mapPVd) {
        std::string Name = std::get<std::string>(V_var[0]);

        if (V_var.size() > 1 && std::holds_alternative<uomap_P_d>(V_var[1])) {
          const auto &smap = std::get<uomap_P_d>(V_var[1]);
          fprintf(fp, "<DataArray NumberOfComponents='1' type='Float32' Name='%s' format='ascii'>\n", Name.c_str());
          for (const auto &pe : all_points) {
            if (pe.ptr) {
              auto it = smap.find(pe.ptr);
              if (it != smap.end() && isFinite(it->second))
                fprintf(fp, "%s ", NumtoString(it->second).c_str());
              else
                fprintf(fp, "NaN ");
            } else {
              double val = 0.;
              bool found = false;
              if (pe.mid_line) {
                if (Name == "direction_info_count") { val = pe.mid_line->debug_direction_info_count; found = true; }
                else if (Name == "contact_faces_count") { val = pe.mid_line->debug_contact_faces_count; found = true; }
                else if (Name == "body_vertices_count") { val = pe.mid_line->debug_body_vertices_count; found = true; }
                else if (Name == "isInContact_pass_count") { val = pe.mid_line->debug_isInContact_pass_count; found = true; }
                else if (pe.is_true_quad_mid) {
                  if (Name == "φ") { val = pe.mid_line->phiphin[0]; found = true; }
                  else if (Name == "φn") { val = pe.mid_line->phiphin[1]; found = true; }
                  else if (Name == "φt") { val = pe.mid_line->phiphin_t[0]; found = true; }
                  else if (Name == "φnt") { val = pe.mid_line->phiphin_t[1]; found = true; }
                  else if (Name == "diag") { val = pe.mid_line->diag_coeff_BEM; found = true; }
                }
              }
              if (!found) {
                auto itA = smap.find(pe.mid_pA);
                auto itB = smap.find(pe.mid_pB);
                if (itA != smap.end() && itB != smap.end() &&
                    isFinite(itA->second) && isFinite(itB->second)) {
                  val = 0.5 * (itA->second + itB->second);
                  found = true;
                }
              }
              if (found && isFinite(val))
                fprintf(fp, "%s ", NumtoString(val).c_str());
              else
                fprintf(fp, "NaN ");
            }
          }
          fprintf(fp, "\n</DataArray>\n");
        } else if (V_var.size() > 1 && std::holds_alternative<uomap_P_Tddd>(V_var[1])) {
          const auto &vmap = std::get<uomap_P_Tddd>(V_var[1]);
          fprintf(fp, "<DataArray NumberOfComponents='3' type='Float32' Name='%s' format='ascii'>\n", Name.c_str());
          for (const auto &pe : all_points) {
            if (pe.ptr) {
              auto it = vmap.find(pe.ptr);
              if (it != vmap.end()) {
                const auto &v = it->second;
                fprintf(fp, "%s %s %s ",
                        isFinite(std::get<0>(v)) ? NumtoString(std::get<0>(v)).c_str() : "NaN",
                        isFinite(std::get<1>(v)) ? NumtoString(std::get<1>(v)).c_str() : "NaN",
                        isFinite(std::get<2>(v)) ? NumtoString(std::get<2>(v)).c_str() : "NaN");
              } else {
                fprintf(fp, "NaN NaN NaN ");
              }
            } else {
              auto itA = vmap.find(pe.mid_pA);
              auto itB = vmap.find(pe.mid_pB);
              if (itA != vmap.end() && itB != vmap.end()) {
                auto avg = 0.5 * (itA->second + itB->second);
                fprintf(fp, "%s %s %s ",
                        isFinite(std::get<0>(avg)) ? NumtoString(std::get<0>(avg)).c_str() : "NaN",
                        isFinite(std::get<1>(avg)) ? NumtoString(std::get<1>(avg)).c_str() : "NaN",
                        isFinite(std::get<2>(avg)) ? NumtoString(std::get<2>(avg)).c_str() : "NaN");
              } else {
                fprintf(fp, "NaN NaN NaN ");
              }
            }
          }
          fprintf(fp, "\n</DataArray>\n");
        }
      }
      fprintf(fp, "</PointData>\n");
    }

    // Cells
    fprintf(fp, "<Cells>\n");
    fprintf(fp, "<DataArray type='Int32' Name='connectivity' format='ascii'>\n");
    for (int i = 0; i < num_points; ++i)
      fprintf(fp, "%d ", i);
    fprintf(fp, "\n</DataArray>\n");
    fprintf(fp, "<DataArray type='Int32' Name='offsets' format='ascii'>\n");
    { int sum = 0; for (const auto &sz : cell_sizes) fprintf(fp, "%d ", sum += sz); }
    fprintf(fp, "\n</DataArray>\n");
    fprintf(fp, "<DataArray type='UInt8' Name='types' format='ascii'>\n");
    for (const auto &t : cell_types)
      fprintf(fp, "%d ", t);
    fprintf(fp, "\n</DataArray>\n");
    fprintf(fp, "</Cells>\n");

    fprintf(fp, "</Piece>\n");
    fprintf(fp, "</UnstructuredGrid>\n");
    fprintf(fp, "</VTKFile>\n");
#if defined(debug_mk_vtu)
    std::cout << Red << "|" << colorReset << std::endl;
#endif
    fclose(fp);
  } catch (std::exception &e) {
    std::cerr << e.what() << colorReset << std::endl;
    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
  }
}

void write_step(const OutputContext &ctx, const std::map<std::string, outputInfo> &NetOutputInfo, const std::vector<Network *> &FluidObject, const std::vector<Network *> &RigidBodyObject, const std::vector<Network *> &SoftBodyObject, const std::unordered_set<networkFace *> & /*allFaces*/) {
  // Fluid meshes
  for (auto *net : FluidObject) {
    auto it = NetOutputInfo.find(net->getName());
    if (it == NetOutputInfo.end())
      continue;
    const auto &info = it->second;

    std::filesystem::path filename = info.vtu_file_name + std::to_string(ctx.time_step) + ".vtu";
    mk_vtu_quadratic(ctx.output_directory / filename, net->getBoundaryFaces(), dataForOutput(net, ctx.dt));
    if (info.PVD) {
      info.PVD->push(filename, ctx.simulation_time);
      info.PVD->output();
    }
  }

  // Fluid tetrahedral meshes
  for (auto *net : FluidObject) {
    auto it = NetOutputInfo.find(net->getName() + "_tetra");
    if (it == NetOutputInfo.end())
      continue;
    const auto &info = it->second;

    std::filesystem::path filename = info.vtu_file_name + std::to_string(ctx.time_step) + ".vtu";
    std::ofstream ofs(ctx.output_directory / filename);
    vtkUnstructuredGridWrite(ofs, net->getTetras());
    if (info.PVD) {
      info.PVD->push(filename, ctx.simulation_time);
      info.PVD->output();
    }
  }

  // Rigid meshes
  for (auto *net : RigidBodyObject) {
    auto it = NetOutputInfo.find(net->getName());
    if (it == NetOutputInfo.end())
      continue;
    const auto &info = it->second;

    std::filesystem::path filename = info.vtu_file_name + std::to_string(ctx.time_step) + ".vtu";
    mk_vtu_quadratic(ctx.output_directory / filename, net->getBoundaryFaces(), dataForOutput(net, ctx.dt));
    if (info.PVD) {
      info.PVD->push(filename, ctx.simulation_time);
      info.PVD->output();
    }

    // Mooring lines (VTP)
    for (auto *mooring : net->mooringLines) {
      auto itm = NetOutputInfo.find(mooring->getName());
      if (itm == NetOutputInfo.end())
        continue;
      const auto &minfo = itm->second;

      std::filesystem::path filename_m = minfo.vtu_file_name + std::to_string(ctx.time_step) + ".vtp";
      std::ofstream ofs(ctx.output_directory / filename_m);
      vtkPolygonWrite(ofs, mooring->getLines());
      ofs.close();
      if (minfo.PVD) {
        minfo.PVD->push(filename_m, ctx.simulation_time);
        minfo.PVD->output();
      }
    }
  }

  // Soft meshes
  for (auto *net : SoftBodyObject) {
    auto it = NetOutputInfo.find(net->getName());
    if (it == NetOutputInfo.end())
      continue;
    const auto &info = it->second;

    std::filesystem::path filename = info.vtu_file_name + std::to_string(ctx.time_step) + ".vtu";
    mk_vtu_quadratic(ctx.output_directory / filename, net->getBoundaryFaces(), dataForOutput(net, ctx.dt));
    if (info.PVD) {
      info.PVD->push(filename, ctx.simulation_time);
      info.PVD->output();
    }
  }
}

void write_vpm(const OutputContext &ctx, const VortexMethod &vpm, PVDWriter &pvd) {
  std::string filename = "vpm_" + std::to_string(ctx.time_step) + ".vtp";
  std::filesystem::path path = ctx.output_directory / filename;

  std::ofstream ofs(path);
  if (!ofs) {
    std::cerr << "Error: Cannot open file for writing: " << path << std::endl;
    return;
  }

  const auto &particles = vpm.getParticles();
  const size_t num_particles = particles.size();

  ofs << "<?xml version=\"1.0\"?>\n";
  ofs << "<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n";
  ofs << "  <PolyData>\n";
  ofs << "    <Piece NumberOfPoints=\"" << num_particles << "\" NumberOfVerts=\"" << num_particles << "\">\n";

  // Points coordinates
  ofs << "      <Points>\n";
  ofs << "        <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for (const auto &p : particles) {
    ofs << "          " << p.x[0] << " " << p.x[1] << " " << p.x[2] << "\n";
  }
  ofs << "        </DataArray>\n";
  ofs << "      </Points>\n";

  // Vertices (cells)
  ofs << "      <Verts>\n";
  ofs << "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
  for (size_t i = 0; i < num_particles; ++i) {
    ofs << "          " << i << "\n";
  }
  ofs << "        </DataArray>\n";
  ofs << "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
  for (size_t i = 0; i < num_particles; ++i) {
    ofs << "          " << i + 1 << "\n";
  }
  ofs << "        </DataArray>\n";
  ofs << "      </Verts>\n";

  // PointData
  ofs << "      <PointData>\n";
  ofs << "        <DataArray type=\"Float64\" Name=\"Vorticity\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for (const auto &p : particles) {
    ofs << "          " << p.alpha[0] << " " << p.alpha[1] << " " << p.alpha[2] << "\n";
  }
  ofs << "        </DataArray>\n";
  ofs << "        <DataArray type=\"Float64\" Name=\"u_total\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for (const auto &p : particles) {
    ofs << "          " << p.u_total[0] << " " << p.u_total[1] << " " << p.u_total[2] << "\n";
  }
  ofs << "        </DataArray>\n";
  ofs << "        <DataArray type=\"Float64\" Name=\"u_potential_BEM\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for (const auto &p : particles) {
    ofs << "          " << p.u_potential_BEM[0] << " " << p.u_potential_BEM[1] << " " << p.u_potential_BEM[2] << "\n";
  }
  ofs << "        </DataArray>\n";
  ofs << "        <DataArray type=\"Float64\" Name=\"u_omega_VPM\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for (const auto &p : particles) {
    ofs << "          " << p.u_omega_VPM[0] << " " << p.u_omega_VPM[1] << " " << p.u_omega_VPM[2] << "\n";
  }
  ofs << "        </DataArray>\n";
  ofs << "        <DataArray type=\"Float64\" Name=\"Sigma\" NumberOfComponents=\"1\" format=\"ascii\">\n";
  for (const auto &p : particles) {
    ofs << "          " << p.sigma << "\n";
  }
  ofs << "        </DataArray>\n";
  ofs << "      </PointData>\n";

  ofs << "    </Piece>\n";
  ofs << "  </PolyData>\n";
  ofs << "</VTKFile>\n";

  ofs.close();

  pvd.push(filename, ctx.simulation_time);
  pvd.output();
}

} // namespace OutputParaView
