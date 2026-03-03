#pragma once

#include <iostream>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// ============================================================================
// Surface collision detection and resolution
//
// Problem A: Non-adjacent surface collision (plunging breaker)
// Problem B: Adjacent face folding (degenerate mesh deformation)
//
// Resolution strategy:
//   1. Detect collision zones (problem faces + 1-ring boundary)
//   2. Tetrahedralize zone points using TetGen (point-cloud Delaunay)
//   3. Classify tetrahedra as water/air using boundary face normals
//   4. Extract water-boundary faces as the new surface
//   5. Replace old problem faces with new surface topology
//   6. Interpolate phi for any new (Steiner) points
// ============================================================================

struct CollisionZone {
   std::unordered_set<networkFace*> problem_faces;
   std::unordered_set<networkFace*> boundary_faces;
   std::unordered_set<networkPoint*> interior_points;
   std::unordered_set<networkPoint*> boundary_points;
   std::unordered_set<networkPoint*> all_points;
   std::unordered_set<networkLine*> all_lines;
   std::unordered_map<networkPoint*, Tdd> saved_phiphin;
   std::unordered_map<networkPoint*, Tdd> saved_phiphin_t;
};

// ---------------------------------------------------------------------------
// Problem B: Adjacent face folding detection
// ---------------------------------------------------------------------------
inline std::unordered_set<networkFace*> detectFoldedFaces(
    const Network& water,
    double normal_reversal_cos) {
   std::unordered_set<networkFace*> folded;

   for (auto* l : water.getBoundaryLines()) {
      auto faces = l->getBoundaryFaces();
      if (faces.size() != 2)
         continue;

      auto* f0 = faces[0];
      auto* f1 = faces[1];
      if (!f0 || !f1)
         continue;

      double dot = Dot(f0->normal, f1->normal);
      if (dot < normal_reversal_cos) {
         folded.insert(f0);
         folded.insert(f1);
      }
   }
   return folded;
}

// ---------------------------------------------------------------------------
// Topological distance check
// ---------------------------------------------------------------------------
inline bool isTopologicallyClose(
    const networkFace* f, const networkFace* g, int max_depth) {
   auto f_points = f->getPoints();
   auto g_points = g->getPoints();

   std::unordered_set<networkPoint*> visited;
   std::vector<networkPoint*> current_ring;
   for (auto* p : f_points) {
      current_ring.push_back(p);
      visited.insert(p);
   }

   for (int depth = 0; depth < max_depth; ++depth) {
      std::vector<networkPoint*> next_ring;
      for (auto* p : current_ring) {
         for (auto* neighbor : p->getNeighbors()) {
            if (visited.count(neighbor))
               continue;
            for (auto* gp : g_points) {
               if (neighbor == gp)
                  return true;
            }
            visited.insert(neighbor);
            next_ring.push_back(neighbor);
         }
      }
      current_ring = std::move(next_ring);
   }
   return false;
}

// ---------------------------------------------------------------------------
// Triangle-to-triangle minimum distance
// ---------------------------------------------------------------------------
inline double triangleTriangleDistance(const networkFace* f, const networkFace* g) {
   auto fp = f->getPoints();
   auto gp = g->getPoints();
   T3Tddd tri_f = {fp[0]->X, fp[1]->X, fp[2]->X};
   T3Tddd tri_g = {gp[0]->X, gp[1]->X, gp[2]->X};

   double min_dist = 1e30;

   for (auto* p : fp) {
      double d = Norm(Nearest(p->X, tri_g) - p->X);
      if (d < min_dist) min_dist = d;
   }
   for (auto* p : gp) {
      double d = Norm(Nearest(p->X, tri_f) - p->X);
      if (d < min_dist) min_dist = d;
   }
   {
      Tddd cf = (fp[0]->X + fp[1]->X + fp[2]->X) / 3.0;
      double d = Norm(Nearest(cf, tri_g) - cf);
      if (d < min_dist) min_dist = d;
   }
   {
      Tddd cg = (gp[0]->X + gp[1]->X + gp[2]->X) / 3.0;
      double d = Norm(Nearest(cg, tri_f) - cg);
      if (d < min_dist) min_dist = d;
   }
   return min_dist;
}

// ---------------------------------------------------------------------------
// Problem A: Non-adjacent surface collision detection
// ---------------------------------------------------------------------------
inline std::unordered_set<networkFace*> detectNonAdjacentCollisions(
    Network& water,
    double proximity_threshold) {
   std::unordered_set<networkFace*> colliding;

   water.makeBucketFaces(proximity_threshold * 2.0);

   auto boundary_faces = water.getBoundaryFaces();

   for (auto* f : boundary_faces) {
      auto fp = f->getPoints();
      Tddd centroid_f = (fp[0]->X + fp[1]->X + fp[2]->X) / 3.0;

      auto nearby = water.BucketFaces.getData(centroid_f, proximity_threshold * 2.0);

      for (auto* g : nearby) {
         if (g == f)
            continue;

         auto gp = g->getPoints();
         bool shares_vertex = false;
         for (auto* p : fp)
            for (auto* q : gp)
               if (p == q) {
                  shares_vertex = true;
                  break;
               }
         if (shares_vertex)
            continue;

         if (isTopologicallyClose(f, g, 2))
            continue;

         double dist = triangleTriangleDistance(f, g);
         if (dist >= proximity_threshold)
            continue;

         Tddd diff = centroid_f - ((gp[0]->X + gp[1]->X + gp[2]->X) / 3.0);
         double d_norm = Norm(diff);
         if (d_norm < 1e-15)
            continue;

         Tddd dir = diff / d_norm;
         if (Dot(f->normal, -dir) > 0 && Dot(g->normal, dir) > 0) {
            colliding.insert(f);
            colliding.insert(g);
         }
      }
   }
   return colliding;
}

// ---------------------------------------------------------------------------
// Build connected collision zones from problem faces
// ---------------------------------------------------------------------------
inline std::vector<CollisionZone> buildCollisionZones(
    const std::unordered_set<networkFace*>& problem_faces,
    int min_zone_faces) {
   std::vector<CollisionZone> zones;
   std::unordered_set<networkFace*> visited;

   for (auto* seed : problem_faces) {
      if (visited.count(seed))
         continue;

      CollisionZone zone;
      std::queue<networkFace*> queue;
      queue.push(seed);
      visited.insert(seed);

      while (!queue.empty()) {
         auto* f = queue.front();
         queue.pop();
         zone.problem_faces.insert(f);

         for (auto* neighbor : f->getNeighbors()) {
            if (neighbor && problem_faces.count(neighbor) && !visited.count(neighbor)) {
               visited.insert(neighbor);
               queue.push(neighbor);
            }
         }
      }

      if (static_cast<int>(zone.problem_faces.size()) < min_zone_faces)
         continue;

      std::unordered_set<networkPoint*> problem_points;
      for (auto* f : zone.problem_faces)
         for (auto* p : f->getPoints())
            problem_points.insert(p);

      for (auto* f : zone.problem_faces) {
         for (auto* neighbor : f->getNeighbors()) {
            if (neighbor && !zone.problem_faces.count(neighbor))
               zone.boundary_faces.insert(neighbor);
         }
      }

      for (auto* f : zone.boundary_faces)
         for (auto* p : f->getPoints())
            zone.boundary_points.insert(p);

      for (auto* p : problem_points) {
         if (!zone.boundary_points.count(p))
            zone.interior_points.insert(p);
      }

      zone.all_points = zone.interior_points;
      zone.all_points.insert(zone.boundary_points.begin(), zone.boundary_points.end());

      for (auto* f : zone.problem_faces)
         for (auto* l : f->getLines())
            zone.all_lines.insert(l);

      for (auto* p : zone.all_points) {
         zone.saved_phiphin[p] = p->phiphin;
         zone.saved_phiphin_t[p] = p->phiphin_t;
      }

      zones.push_back(std::move(zone));
   }
   return zones;
}

// ===========================================================================
// Phase 2: TetGen-based collision resolution
// ===========================================================================

#ifdef tetgenH

// ---------------------------------------------------------------------------
// Resolve a single collision zone using TetGen tetrahedralization.
//
// Algorithm:
//   1. Tetrahedralize zone points (Delaunay point cloud)
//   2. Create tetrahedra in the Network using genTetra()
//   3. Classify each tetrahedron as water or air
//      (centroid on water side of nearest boundary face)
//   4. Identify new surface faces (faces with exactly one water tetra)
//   5. Delete all zone tetrahedra
//   6. Delete old problem faces not reused by new surface
//   7. Clean up orphaned lines and interior points
//   8. Interpolate phi for any Steiner points
// ---------------------------------------------------------------------------
inline bool resolveCollisionZone(
    Network& water,
    CollisionZone& zone,
    int time_step,
    int zone_idx) {

   if (zone.all_points.size() < 4) {
      std::cout << "  zone " << zone_idx << ": too few points (" << zone.all_points.size() << "), skipping" << std::endl;
      return false;
   }

   // === Step 1: Prepare TetGen input from zone points ===
   std::vector<networkPoint*> points_vec(zone.all_points.begin(), zone.all_points.end());
   tetgenio in = generate_tetgenio_input_from_points(points_vec);

   // === Step 2: Run TetGen (point-cloud Delaunay) ===
   tetgenbehavior b;
   b.parse_commandline(const_cast<char*>(""));
   tetgenio out;
   try {
      ::tetrahedralize(&b, &in, &out);
   } catch (...) {
      std::cerr << Red << "[collision] zone " << zone_idx
                << ": TetGen failed, skipping resolution" << colorReset << std::endl;
      return false;
   }

   if (out.numberoftetrahedra == 0) {
      std::cout << "  zone " << zone_idx << ": TetGen produced 0 tetrahedra, skipping" << std::endl;
      return false;
   }

   std::cout << "  zone " << zone_idx << ": TetGen produced "
             << out.numberoftetrahedra << " tetrahedra from "
             << points_vec.size() << " points" << std::endl;

   // === Step 3: Build index-to-point map ===
   // generate_tetgenio_input_from_points sets p->index for each point
   std::vector<networkPoint*> idx_to_point(out.numberofpoints, nullptr);
   for (auto* p : points_vec) {
      idx_to_point[p->index] = p;
   }

   // Handle Steiner points (unlikely for unconstrained Delaunay, but be safe)
   std::vector<networkPoint*> steiner_points;
   if (out.numberofpoints > static_cast<int>(points_vec.size())) {
      water.makeBucketPoints(water.getScale() / 50.);
      for (int i = static_cast<int>(points_vec.size()); i < out.numberofpoints; i++) {
         Tddd X = {out.pointlist[i * 3], out.pointlist[i * 3 + 1], out.pointlist[i * 3 + 2]};
         auto* p = new networkPoint(&water, X);
         water.BucketPoints.add(X, p);
         idx_to_point[i] = p;
         steiner_points.push_back(p);
      }
      std::cout << "  zone " << zone_idx << ": " << steiner_points.size() << " Steiner points added" << std::endl;
   }

   // === Step 4: Classify tetrahedra as water/air ===
   // Use boundary faces (healthy faces surrounding the zone) for classification.
   // Normal convention: normal points from water -> air (outward).
   // Dot(centroid - nearest_surface_point, face->normal) < 0 => water side.
   std::vector<networkFace*> ref_faces(zone.boundary_faces.begin(), zone.boundary_faces.end());

   // Build a quick spatial structure for reference faces
   // (zone boundary faces are usually small in number, linear scan is OK)
   auto classifyPoint = [&](const Tddd& X) -> bool {
      double min_dist = 1e30;
      networkFace* nearest_face = nullptr;
      Tddd nearest_pt;
      for (auto* f : ref_faces) {
         Tddd pt = Nearest(X, f);
         double d = Norm(X - pt);
         if (d < min_dist) {
            min_dist = d;
            nearest_pt = pt;
            nearest_face = f;
         }
      }
      if (!nearest_face) return false;
      return Dot(X - nearest_pt, nearest_face->normal) < 0;
   };

   std::vector<bool> is_water(out.numberoftetrahedra, false);
   for (int i = 0; i < out.numberoftetrahedra; i++) {
      Tddd centroid = {0, 0, 0};
      for (int j = 0; j < 4; j++) {
         int idx = out.tetrahedronlist[i * 4 + j];
         centroid += Tddd{out.pointlist[idx * 3], out.pointlist[idx * 3 + 1], out.pointlist[idx * 3 + 2]};
      }
      centroid /= 4.0;
      is_water[i] = classifyPoint(centroid);
   }

   int n_water = 0;
   for (auto w : is_water)
      if (w) n_water++;
   std::cout << "  zone " << zone_idx << ": " << n_water << " water / "
             << (out.numberoftetrahedra - n_water) << " air tetrahedra" << std::endl;

   if (n_water == 0) {
      std::cout << Red << "  zone " << zone_idx
                << ": no water tetrahedra found, skipping resolution" << colorReset << std::endl;
      // Clean up Steiner points
      for (auto* p : steiner_points) delete p;
      return false;
   }

   // === Step 5: Create tetrahedra in the Network ===
   // genTetra() creates faces and lines automatically (reusing existing ones).
   // At this point, all existing tetrahedra have been deleted by DeleteInteriorTetras(),
   // so all faces have Tetras = {nullptr, nullptr}.
   std::vector<networkTetra*> zone_tetras;
   zone_tetras.reserve(out.numberoftetrahedra);

   for (int i = 0; i < out.numberoftetrahedra; i++) {
      auto* p0 = idx_to_point[out.tetrahedronlist[i * 4 + 0]];
      auto* p1 = idx_to_point[out.tetrahedronlist[i * 4 + 1]];
      auto* p2 = idx_to_point[out.tetrahedronlist[i * 4 + 2]];
      auto* p3 = idx_to_point[out.tetrahedronlist[i * 4 + 3]];

      if (!p0 || !p1 || !p2 || !p3) continue;

      try {
         auto [success, tet] = genTetra(&water, p0, p1, p2, p3);
         if (tet) {
            zone_tetras.push_back(tet);
         }
      } catch (...) {
         // Skip degenerate tetrahedra
      }
   }

   // === Step 6: Identify new surface faces ===
   // A face is on the new water surface if exactly one of its two adjacent
   // tetrahedra is classified as water (the other is air or nullptr).
   std::unordered_set<networkTetra*> water_tetra_set;
   for (size_t i = 0; i < zone_tetras.size(); i++) {
      // Map zone_tetras index back to TetGen index for water classification
      // Since genTetra might skip some, we need to track the mapping.
      // However, genTetra creates tetras in order, so we track via a separate index.
   }
   // Re-classify using the actually-created tetras.
   // We classify each tetra by its centroid position.
   for (auto* tet : zone_tetras) {
      auto [p0, p1, p2, p3] = tet->Points;
      Tddd centroid = (p0->X + p1->X + p2->X + p3->X) / 4.0;
      if (classifyPoint(centroid)) {
         water_tetra_set.insert(tet);
      }
   }

   // Collect all faces involved in zone tetrahedra
   std::unordered_set<networkFace*> zone_tetra_faces;
   for (auto* tet : zone_tetras) {
      for (auto* f : tet->Faces) {
         if (f) zone_tetra_faces.insert(f);
      }
   }

   // Identify new surface faces
   std::unordered_set<networkFace*> new_surface_faces;
   for (auto* f : zone_tetra_faces) {
      auto [t0, t1] = f->Tetras;
      bool w0 = t0 && water_tetra_set.count(t0);
      bool w1 = t1 && water_tetra_set.count(t1);
      // Surface if exactly one side is water
      if (w0 != w1) {
         new_surface_faces.insert(f);
      }
   }

   std::cout << "  zone " << zone_idx << ": " << new_surface_faces.size()
             << " new surface faces identified" << std::endl;

   if (new_surface_faces.empty()) {
      std::cout << Red << "  zone " << zone_idx
                << ": no new surface faces, aborting resolution" << colorReset << std::endl;
      // Delete all zone tetrahedra to restore state
      for (auto* tet : zone_tetras) delete tet;
      for (auto* p : steiner_points) delete p;
      return false;
   }

   // === Step 7: Build set of faces to keep ===
   std::unordered_set<networkFace*> faces_to_keep;
   faces_to_keep.insert(new_surface_faces.begin(), new_surface_faces.end());
   faces_to_keep.insert(zone.boundary_faces.begin(), zone.boundary_faces.end());

   // === Step 8: Delete all zone tetrahedra ===
   // The tetra destructor properly clears face->Tetras pointers.
   for (auto* tet : zone_tetras) {
      delete tet;
   }
   zone_tetras.clear();

   // === Step 9: Delete old problem faces not in keep set ===
   // Also delete faces created by genTetra that are not part of the new surface.
   std::unordered_set<networkFace*> faces_to_delete;
   for (auto* f : zone.problem_faces) {
      if (!faces_to_keep.count(f))
         faces_to_delete.insert(f);
   }
   for (auto* f : zone_tetra_faces) {
      if (!faces_to_keep.count(f))
         faces_to_delete.insert(f);
   }
   // Don't delete boundary faces that are outside the zone
   for (auto* f : zone.boundary_faces) {
      faces_to_delete.erase(f);
   }

   for (auto* f : faces_to_delete) {
      delete f;
   }

   // === Step 10: Delete orphaned lines ===
   // Lines from the original problem zone that now have no faces
   // Also check lines involving interior points
   std::unordered_set<networkLine*> lines_to_check;
   for (auto* l : zone.all_lines) {
      lines_to_check.insert(l);
   }
   // Also check lines of interior points (new lines might have been created)
   for (auto* p : zone.interior_points) {
      for (auto* l : p->getLines())
         lines_to_check.insert(l);
   }
   for (auto* p : steiner_points) {
      for (auto* l : p->getLines())
         lines_to_check.insert(l);
   }

   for (auto* l : lines_to_check) {
      if (l->getFaces().empty()) {
         delete l;
      }
   }

   // === Step 11: Delete orphaned interior points ===
   for (auto* p : zone.interior_points) {
      if (p->getLines().empty()) {
         delete p;
      }
   }
   // Delete orphaned Steiner points
   std::vector<networkPoint*> surviving_steiner;
   for (auto* p : steiner_points) {
      if (p->getLines().empty()) {
         delete p;
      } else {
         surviving_steiner.push_back(p);
      }
   }

   // === Step 12: Interpolate phi for Steiner points ===
   for (auto* p : surviving_steiner) {
      auto neighbors = p->getNeighborPointsOnSurfaces();
      if (neighbors.empty()) continue;

      Tdd avg_phiphin = {0, 0};
      Tdd avg_phiphin_t = {0, 0};
      int count = 0;
      for (auto* q : neighbors) {
         // Use saved values from before resolution (more reliable)
         auto it = zone.saved_phiphin.find(q);
         if (it != zone.saved_phiphin.end()) {
            avg_phiphin += it->second;
            avg_phiphin_t += zone.saved_phiphin_t[q];
         } else {
            avg_phiphin += q->phiphin;
            avg_phiphin_t += q->phiphin_t;
         }
         count++;
      }
      if (count > 0) {
         p->phiphin = avg_phiphin / static_cast<double>(count);
         p->phiphin_t = avg_phiphin_t / static_cast<double>(count);
      }
   }

   std::cout << Green << "  zone " << zone_idx << ": resolution complete ("
             << new_surface_faces.size() << " new surface faces, "
             << surviving_steiner.size() << " Steiner points)" << colorReset << std::endl;

   return true;
}

#endif // tetgenH

// ---------------------------------------------------------------------------
// Top-level collision detection and resolution
// ---------------------------------------------------------------------------
inline void detectAndResolveCollisions(
    Network& water,
    int time_step,
    const SimulationSettings::RemeshingSettings::CollisionSettings& settings) {
   if (!settings.enabled)
      return;

   double mean_len = Mean(extLength(water.getLines()));
   double threshold = mean_len * settings.proximity_factor;

   // Phase 1: Detection
   std::unordered_set<networkFace*> problem_faces;

   if (settings.detect_folding) {
      auto folded = detectFoldedFaces(water, settings.normal_reversal_cos);
      if (!folded.empty()) {
         std::cout << Red << "[collision] time_step " << time_step
                   << ": detected " << folded.size()
                   << " folded faces (adjacent face folding)" << colorReset << std::endl;
         problem_faces.insert(folded.begin(), folded.end());
      }
   }

   if (settings.detect_non_adjacent) {
      auto colliding = detectNonAdjacentCollisions(water, threshold);
      if (!colliding.empty()) {
         std::cout << Red << "[collision] time_step " << time_step
                   << ": detected " << colliding.size()
                   << " colliding faces (non-adjacent collision, threshold=" << threshold << ")"
                   << colorReset << std::endl;
         problem_faces.insert(colliding.begin(), colliding.end());
      }
   }

   if (problem_faces.empty()) {
      return;
   }

   // Phase 2: Build collision zones
   auto zones = buildCollisionZones(problem_faces, settings.min_zone_faces);
   std::cout << Green << "[collision] time_step " << time_step
             << ": " << zones.size() << " collision zone(s) identified" << colorReset << std::endl;

   for (size_t i = 0; i < zones.size(); ++i) {
      const auto& zone = zones[i];
      std::cout << "  zone " << i << ": "
                << zone.problem_faces.size() << " problem faces, "
                << zone.boundary_faces.size() << " boundary faces, "
                << zone.interior_points.size() << " interior points, "
                << zone.boundary_points.size() << " boundary points"
                << std::endl;
   }

   // Phase 3: Resolution (TetGen-based)
#ifdef tetgenH
   int resolved_count = 0;
   for (size_t i = 0; i < zones.size(); ++i) {
      try {
         if (resolveCollisionZone(water, zones[i], time_step, static_cast<int>(i))) {
            resolved_count++;
         }
      } catch (std::exception& e) {
         std::cerr << Red << "[collision] zone " << i
                   << ": resolution failed with exception: " << e.what()
                   << colorReset << std::endl;
      }
   }

   if (resolved_count > 0) {
      water.setGeometricPropertiesForce();
      water.makeBuckets();
      std::cout << Green << "[collision] time_step " << time_step
                << ": resolved " << resolved_count << " / " << zones.size()
                << " collision zone(s)" << colorReset << std::endl;
   }
#else
   std::cout << "[collision] TetGen not available, resolution skipped" << std::endl;
#endif
}
