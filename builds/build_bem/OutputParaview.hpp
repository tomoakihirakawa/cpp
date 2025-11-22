#pragma once
#include "BEM.hpp"
#include "OutputCommon.hpp"
#include "vtkWriter.hpp"

namespace OutputParaView {

void write_step(const OutputContext &ctx, const std::map<std::string, outputInfo> &NetOutputInfo, const std::vector<Network *> &FluidObject, const std::vector<Network *> &RigidBodyObject, const std::vector<Network *> &SoftBodyObject, const std::unordered_set<networkFace *> & /*allFaces*/) {
  // Fluid meshes
  for (auto *net : FluidObject) {
    auto it = NetOutputInfo.find(net->getName());
    if (it == NetOutputInfo.end())
      continue;
    const auto &info = it->second;

    std::filesystem::path filename = info.vtu_file_name + std::to_string(ctx.time_step) + ".vtu";
    mk_vtu(ctx.output_directory / filename, net->getSurfaces(), dataForOutput(net, ctx.dt));
    info.PVD->push(filename, ctx.simulation_time);
    info.PVD->output();
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
    info.PVD->push(filename, ctx.simulation_time);
    info.PVD->output();
  }

  // Rigid meshes
  for (auto *net : RigidBodyObject) {
    auto it = NetOutputInfo.find(net->getName());
    if (it == NetOutputInfo.end())
      continue;
    const auto &info = it->second;

    std::filesystem::path filename = info.vtu_file_name + std::to_string(ctx.time_step) + ".vtu";
    mk_vtu(ctx.output_directory / filename, net->getSurfaces(), dataForOutput(net, ctx.dt));
    info.PVD->push(filename, ctx.simulation_time);
    info.PVD->output();

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
      minfo.PVD->push(filename_m, ctx.simulation_time);
      minfo.PVD->output();
    }
  }

  // Soft meshes
  for (auto *net : SoftBodyObject) {
    auto it = NetOutputInfo.find(net->getName());
    if (it == NetOutputInfo.end())
      continue;
    const auto &info = it->second;

    std::filesystem::path filename = info.vtu_file_name + std::to_string(ctx.time_step) + ".vtu";
    mk_vtu(ctx.output_directory / filename, net->getSurfaces(), dataForOutput(net, ctx.dt));
    info.PVD->push(filename, ctx.simulation_time);
    info.PVD->output();
  }
}

} // namespace OutputParaView