#pragma once
#include "BEM.hpp"
#include "OutputCommon.hpp"

namespace OutputJSON {

void write_step(
    const OutputContext& ctx,
    JSONoutput& jsonout,
    const std::vector<Network*>& FluidObject,
    const std::vector<Network*>& RigidBodyObject,
    const std::vector<Network*>& SoftBodyObject,
    const std::vector<JSON>& MeasurementJSONs,
    const std::unordered_set<networkFace*>& allFaces,
    const std::vector<double>& eq_of_motion,
    double unknownsize,
    double time_setIGIGn,
    double time_solve) {
   // timing
   const double cpu_time_in_seconds = static_cast<double>(std::clock() - ctx.cpu_clock_start) / CLOCKS_PER_SEC;
   const double wall_clock_time_in_seconds = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - ctx.wall_clock_start).count();

   // basic
   jsonout.push("simulation_time", ctx.simulation_time);
   jsonout.push("time_step", ctx.time_step);
   jsonout.push("dt", ctx.dt);
   jsonout.push("eq_of_motion", eq_of_motion);
   jsonout.push("unknownsize", unknownsize);
   jsonout.push("time_setIGIGn", time_setIGIGn);
   jsonout.push("time_solve", time_solve);
   jsonout.push("cpu_time", cpu_time_in_seconds);
   jsonout.push("wall_clock_time", wall_clock_time_in_seconds);

   // fluid mesh size
   for (const auto* water : FluidObject) {
      jsonout.push(water->getName() + "_face_size", static_cast<double>(water->getSurfaces().size()));
      jsonout.push(water->getName() + "_point_size", static_cast<double>(water->getPoints().size()));
   }

   // rigid/soft basic states (safe subset)
   for (const auto* net : RigidBodyObject) {
      jsonout.push(net->getName() + "_COM", net->COM);
      jsonout.push(net->getName() + "_velocity", net->velocity);
      jsonout.push(net->getName() + "_accel", net->acceleration);
   }
   for (const auto* net : SoftBodyObject) {
      jsonout.push(net->getName() + "_point_size", static_cast<double>(net->getPoints().size()));
      jsonout.push(net->getName() + "_face_size", static_cast<double>(net->getSurfaces().size()));
   }

   // measurements: line intersections (optional)
   for (const auto& js : MeasurementJSONs) {
      js.find("position", [&](const std::vector<std::string>& value) {
         if (value.size() != 6)
            throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "position must have 6 elements");
         const Tddd X0 = {std::stod(value[0]), std::stod(value[1]), std::stod(value[2])};
         const Tddd X1 = {std::stod(value[3]), std::stod(value[4]), std::stod(value[5])};
         const auto name_key = js.at("name")[0] + "_intersection";

         Tddd nearest_intersection = {1E+20, 1E+20, 1E+20};
         const Tddd mid = (X0 + X1) * 0.5;

         for (auto* f : allFaces) {
            auto [isintersect, X, _] = IntersectQ_(T2Tddd{X0, X1}, (T3Tddd)(*f));
            if (isintersect && Norm(X - mid) < Norm(nearest_intersection - mid))
               nearest_intersection = X;
         }
         jsonout.push(name_key, nearest_intersection);
      });
   }

   // write result.json
   std::ofstream os(ctx.output_directory / "result.json");
   jsonout.output(os);
   os.close();
}

}  // namespace OutputJSON