#pragma once

#include <cmath>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "BEM_collision.hpp"

inline bool flipIfOnce(Network& water, const Tdd& limit_Dirichlet, const Tdd& limit_Neumann, bool force = false, int iteration = 0) {
  try {
    auto [target_of_max_normal_diffD, acceptable_normal_change_by_flipD] = limit_Dirichlet;
    auto [target_of_max_normal_diffN, acceptable_normal_change_by_flipN] = limit_Neumann;
    std::cout << "flipIf" << std::endl;
    water.setGeometricPropertiesForce();
    int count = 0;
    auto V = ToVector(water.getLines());
    for (const auto& l : RandomSample(V)) {
      if (!l->CORNER) {
        if (force && (iteration == 0 || count < iteration)) {
          if (l->Dirichlet) {
            if (l->flipIfTopologicallyBetter(target_of_max_normal_diffD, acceptable_normal_change_by_flipD))
              count++;
          } else {
            if (l->flipIfTopologicallyBetter(target_of_max_normal_diffN, acceptable_normal_change_by_flipN))
              count++;
          }
        } else {
          if (l->Dirichlet) {
            if (l->flipIfBetter(target_of_max_normal_diffD, acceptable_normal_change_by_flipD, 5))
              count++;
          } else {
            if (l->flipIfBetter(target_of_max_normal_diffN, acceptable_normal_change_by_flipN, 4))
              count++;
          }
        }
      }
    }
    if (count > 0) {
      std::cout << Green << "  " << count << " edges flipped." << colorReset << std::endl;
      water.setGeometricPropertiesForce();
      water.checkConnectivity();
      return true;
    }
    std::cout << Green << "  No edges flipped." << colorReset << std::endl;
    return false;
  } catch (std::exception& e) {
    std::cerr << e.what() << colorReset << std::endl;
    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
  };
}

inline bool flipIfBatched(Network& water, const Tdd& limit_Dirichlet, const Tdd& limit_Neumann, const char* tag = nullptr) {
  auto [target_of_max_normal_diffD, acceptable_normal_change_by_flipD] = limit_Dirichlet;
  auto [target_of_max_normal_diffN, acceptable_normal_change_by_flipN] = limit_Neumann;
  std::cout << "flipIf (batched)" << (tag ? std::string(" ") + tag : "") << std::endl;
  water.setGeometricPropertiesForce();
  auto line_alive = [&](const networkLine* l) { return l && (water.Lines.find(const_cast<networkLine*>(l)) != water.Lines.end()); };
  auto collect_non_adjacent = [&](const std::vector<networkLine*>& candidates) {
    std::vector<networkLine*> batch;
    std::unordered_set<networkPoint*> used_points;
    batch.reserve(candidates.size());
    used_points.reserve(candidates.size() * 2);
    for (auto* l : candidates) {
      if (!l)
        continue;
      auto [p0, p1] = l->getPoints();
      if (!p0 || !p1)
        continue;
      if (used_points.contains(p0) || used_points.contains(p1))
        continue;
      used_points.insert(p0);
      used_points.insert(p1);
      batch.emplace_back(l);
    }
    return batch;
  };
  auto gather_neighbor_lines = [&](const networkLine* l, std::unordered_set<networkLine*>& out) {
    if (!l)
      return;
    auto [p0, p1] = l->getPoints();
    auto add_from_point = [&](networkPoint* p) {
      if (!p)
        return;
      for (auto* nl : p->getBoundaryLines())
        if (nl)
          out.insert(nl);
    };
    add_from_point(p0);
    add_from_point(p1);
    for (auto* f : l->getBoundaryFaces())
      if (f)
        for (auto* nl : f->getLines())
          if (nl)
            out.insert(nl);
  };

  int total_flipped = 0;
  std::unordered_set<networkLine*> dirty;
  for (auto* l : water.getBoundaryLines())
    dirty.insert(l);
  for (auto iter = 0; iter < 20; iter++) {
    std::vector<networkLine*> candidates;
    candidates.reserve(dirty.size());
    for (auto* l : dirty) {
      if (!line_alive(l) || l->CORNER)
        continue;
      candidates.emplace_back(l);
    }
    if (candidates.empty())
      break;
    auto shuffled = RandomSample(candidates);
    auto batch = collect_non_adjacent(shuffled);
    if (batch.empty())
      break;
    int flipped_in_batch = 0;
    std::unordered_set<networkLine*> touched;
    for (auto* l : batch) {
      if (!line_alive(l) || l->CORNER)
        continue;
      bool flipped = false;
      if (l->Dirichlet) {
        if (l->flipIfBetter(target_of_max_normal_diffD, acceptable_normal_change_by_flipD, 5))
          flipped = true;
      } else {
        if (l->flipIfBetter(target_of_max_normal_diffN, acceptable_normal_change_by_flipN, 4))
          flipped = true;
      }
      if (flipped) {
        flipped_in_batch++;
        total_flipped++;
        gather_neighbor_lines(l, touched);
      }
    }
    if (flipped_in_batch == 0)
      break;
    std::cout << Green << "  [batched] iter=" << iter << " batch_size=" << batch.size() << " flipped=" << flipped_in_batch << colorReset << std::endl;
    water.setGeometricPropertiesForce();
    water.checkConnectivity();
    dirty = std::move(touched);
    if (dirty.empty())
      break;
  }
  if (total_flipped > 0) {
    std::cout << Green << "  " << total_flipped << " edges flipped (batched)." << colorReset << std::endl;
    return true;
  }
  std::cout << Green << "  No edges flipped (batched)." << colorReset << std::endl;
  return false;
}

inline void remesh_for_main_loop(Network& water, const int time_step, const double min_edge_length, const bool tetrahedralize, const bool surface_flip,
                                 const SimulationSettings::RemeshingSettings::CollisionSettings& collision_settings = {}) {
  const double rad = M_PI / 180.0;
  const double cos_3rad = std::cos(3.0 * rad);
  const double cos_20rad = std::cos(20.0 * rad);
  const double cos_rad = std::cos(rad);
  const double global_mean_len = Mean(extLength(water.getLines()));
  const double limit_len = (min_edge_length > 0.0) ? min_edge_length : global_mean_len * 0.1;

  /* --------------------------------- 四面体の削除 --------------------------------- */
  std::cout << "孤立した四面体の削除を開始" << std::endl;
  water.DeleteIsolatedTetras();
  /* -------------------------------------------------------------------------- */
  std::cout << "内部要素の削除を開始" << std::endl;
  water.DeleteInteriorTetras();

  // Surface collision detection and resolution
  detectAndResolveCollisions(water, time_step, collision_settings);

  water.setGeometricPropertiesForce();
  water.checkConnectivity();

  const int iter_divide_collapse = 3;
  for (auto i = 0; i < iter_divide_collapse; i++) {

    // 適当な周辺の線分長さの平均を計算
    // Bench (BEM_BENCH_STEPS=3, remesh_total avg, baseline=1.00):
    // - collapse local_mean_len cache: 1.14 (slower)
    // - local_line_length vector+sort unique: 1.02 (slower)
    // - skip post-split/post-collapse setGeometricProperties/checkConnectivity when no change: 0.96 (faster)
    auto local_line_length = [](const networkLine* const l) {
      static std::unordered_set<networkPoint*> adjacent_points;
      adjacent_points.clear();
      for (const auto& f : l->getBoundaryFaces()) {
        auto points = f->getPoints();
        adjacent_points.insert(points.begin(), points.end());
      }
      static std::unordered_set<networkLine*> adjacent_lines;
      adjacent_lines.clear();
      for (const auto& p : adjacent_points)
        for (const auto& L : p->getBoundaryLines())
          if (L != l)
            adjacent_lines.insert(L);
      auto ret = 0.;
      for (const auto& L : adjacent_lines)
        ret += L->length();
      return ret / adjacent_lines.size();
    };
    /* --------------------------------- split (batched) --------------------------------- */
    auto line_alive = [&](const networkLine* l) { return l && (water.Lines.find(const_cast<networkLine*>(l)) != water.Lines.end()); };
    auto collect_non_adjacent = [&](const std::vector<networkLine*>& candidates) {
      std::vector<networkLine*> batch;
      std::unordered_set<networkPoint*> used_points;
      batch.reserve(candidates.size());
      used_points.reserve(candidates.size() * 2);
      for (auto* l : candidates) {
        if (!l)
          continue;
        auto [p0, p1] = l->getPoints();
        if (!p0 || !p1)
          continue;
        if (used_points.contains(p0) || used_points.contains(p1))
          continue;
        used_points.insert(p0);
        used_points.insert(p1);
        batch.emplace_back(l);
      }
      return batch;
    };
    auto gather_neighbor_lines = [&](const networkLine* l, std::unordered_set<networkLine*>& out) {
      if (!l)
        return;
      auto [p0, p1] = l->getPoints();
      auto add_from_point = [&](networkPoint* p) {
        if (!p)
          return;
        for (auto* nl : p->getBoundaryLines())
          if (nl)
            out.insert(nl);
      };
      add_from_point(p0);
      add_from_point(p1);
      for (auto* f : l->getBoundaryFaces())
        if (f)
          for (auto* nl : f->getLines())
            if (nl)
              out.insert(nl);
    };
    auto should_split = [&](networkLine* l, double local_mean_len) {
      if (!l)
        return false;
      if (local_mean_len <= 0.)
        return false;
      auto len = l->length();
      if (len <= 2.0 * limit_len || len <= 1.4 * local_mean_len)
        return false;
      auto surfaces = l->getBoundaryFaces();
      if (surfaces.size() != 2)
        return false;

      // Split連鎖を防ぐ: 隣接線との長さ比が大きすぎる場合のみsplit
      // 最短の隣接線の2倍以上でなければsplitしない
      double min_neighbor_len = 1e20;
      for (const auto& f : surfaces) {
        for (const auto& nl : f->getLines()) {
          if (nl != l && nl->length() > 0.)
            min_neighbor_len = std::min(min_neighbor_len, static_cast<double>(nl->length()));
        }
      }
      if (min_neighbor_len < 1e19 && len < 2.0 * min_neighbor_len)
        return false;

      return (l->Dirichlet || l->CORNER || Dot(surfaces[0]->normal, surfaces[1]->normal) < cos_3rad);
    };

    {
      std::unordered_set<networkLine*> dirty;
      for (auto* l : water.getBoundaryLines())
        dirty.insert(l);
      for (auto iter = 0; iter < 20; iter++) {
        std::vector<networkLine*> candidates;
        candidates.reserve(dirty.size());
        std::unordered_map<const networkLine*, double> local_mean_cache;
        local_mean_cache.reserve(dirty.size());
        for (auto* l : dirty) {
          if (!line_alive(l))
            continue;
          auto it = local_mean_cache.find(l);
          auto local_mean_len = (it != local_mean_cache.end()) ? it->second : local_line_length(l);
          if (it == local_mean_cache.end())
            local_mean_cache.emplace(l, local_mean_len);
          if (should_split(l, local_mean_len))
            candidates.emplace_back(l);
        }
        if (candidates.empty())
          break;
        auto batch = collect_non_adjacent(candidates);
        if (batch.empty())
          break;
        bool divided_any = false;
        std::unordered_set<networkLine*> touched;
        for (auto* l : batch) {
          if (!line_alive(l))
            continue;
          auto len = l->length();
          auto local_mean_len = local_line_length(l);
          if (!should_split(l, local_mean_len))
            continue;
          l->Split();
          divided_any = true;
          std::cout << Red << "time_step " << time_step << ": line divided due to large length. length = " << len << ", local mean length = " << local_mean_len << colorReset << std::endl;
          gather_neighbor_lines(l, touched);
        }
        if (!divided_any)
          break;
        water.setGeometricPropertiesForce();
        water.checkConnectivity();
        dirty = std::move(touched);
        if (dirty.empty())
          break;
      }
    }

    water.setGeometricPropertiesForce();
    water.checkConnectivity();
    for (const auto& l : water.getLines())
      if (!l->checkTopology())
        throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "topology error detected in line after division");

    /* ---------------------------------------------------------------------------- */

    if (surface_flip) {
      flipIfBatched(water, {10 * rad /*target n diff*/, 10 * rad /*change n diff*/}, {3 * rad, 3 * rad}, "pre-collapse"); // Stricter for Neumann
    }

    /* ---------------------------------- collapse --------------------------------- */

    for (auto iter = 0; iter < 20; iter++) {
      std::unordered_set<networkLine*> dirty;
      for (auto* l : water.getBoundaryLines())
        dirty.insert(l);
      bool found_small_line = false;
      while (!dirty.empty()) {
        std::vector<networkLine*> candidates;
        candidates.reserve(dirty.size());
        for (auto* l : dirty) {
          if (!line_alive(l))
            continue;
          candidates.emplace_back(l);
        }
        if (candidates.empty())
          break;
        auto batch = collect_non_adjacent(candidates);
        if (batch.empty())
          break;
        bool changed_in_batch = false;
        std::unordered_set<networkLine*> touched;
        for (auto* l : batch) {
          if (!line_alive(l))
            continue;

          // Skip CORNER lines (similar to flip operations)
          if (l->CORNER)
            continue;

          auto surfaces = l->getBoundaryFaces();
          if (surfaces.size() != 2)
            continue;
          auto f0 = surfaces[0];
          auto f1 = surfaces[1];
          if (!f0 || !f1)
            continue;
          auto len = l->length();

          auto [a, b, c] = f0->getPoints(l);
          auto [_, __, d] = f1->getPoints(l);

          auto local_mean_len = local_line_length(l);

          // Neumann lines should only be collapsed if extremely small to prevent large normal changes
          if (l->Neumann && len > limit_len * 0.2)
            continue;

          // 線の長さが平均の0.4倍以下ならマージ (dtの低下を防ぐため閾値を上げる)
          if ((len < local_mean_len * 0.4 || len < limit_len)) {
            gather_neighbor_lines(l, touched);
            l->Collapse();
            found_small_line = true;
            changed_in_batch = true;
            std::cout << "time_step " << time_step << ": line merged due to small length. length = " << len << ", local mean length = " << local_mean_len << std::endl;
            continue;
          }
          //
          double volume = TetrahedronVolume(a->X, b->X, c->X, d->X);
          double ref_volume = std::pow(local_mean_len, 3) * 0.5;
          (void)volume;
          (void)ref_volume;

          auto [p0, p1, p2] = f0->getPoints(l);
          auto ci0 = CircumradiusToInradius(p0->X, p1->X, p2->X);
          auto [q0, q1, q2] = f1->getPoints(l);
          auto ci1 = CircumradiusToInradius(q0->X, q1->X, q2->X);

          // For Neumann lines, check that normals are nearly opposite before allowing collapse
          const double cos_5rad = std::cos(5.0 * rad);
          if (l->Neumann) {
            // Only allow collapse if normals are truly opposite (angle > 175 degrees)
            if (Dot(f0->normal, -f1->normal) < cos_5rad)
              continue;
          }

          // 接する面の法線が反対方向でかつ近ければマージ
          if (l->length() < 0.1 * (Norm(p2->X - p1->X) + Norm(p2->X - p0->X) + Norm(q2->X - q1->X) + Norm(q2->X - q0->X)) / 4.)
            if (l->Collapse()) {
              found_small_line = true;
              changed_in_batch = true;
              std::cout << "time_step " << time_step << ": line merged due to opposite normals. n0 = " << f0->normal << ", n1 = " << f1->normal << std::endl;
              continue;
            }

          if (ci0 > 10. || ci1 > 10.)
            if (l->Collapse()) {
              found_small_line = true;
              changed_in_batch = true;
              std::cout << "time_step " << time_step << ": line merged due to opposite normals. n0 = " << f0->normal << ", n1 = " << f1->normal << std::endl;
              continue;
            }
          if (Norm(p2->X - q2->X) < std::sqrt(TriangleArea(p0->X, p1->X, p2->X) + TriangleArea(q0->X, q1->X, q2->X)) / 10.) {
            l->Flip(true); // できればでいい
            l->Collapse();
            found_small_line = true;
            changed_in_batch = true;
            std::cout << "time_step " << time_step << ": line merged due to opposite normals. n0 = " << f0->normal << ", n1 = " << f1->normal << std::endl;
            continue;
          }
          if (Dot(f0->normal, -f1->normal) > cos_rad)
            if (l->Collapse()) {
              found_small_line = true;
              changed_in_batch = true;
              std::cout << "time_step " << time_step << ": line merged due to opposite normals. n0 = " << f0->normal << ", n1 = " << f1->normal << std::endl;
              continue;
            }
        }
        if (changed_in_batch) {
          water.setGeometricPropertiesForce();
          water.checkConnectivity();
        }
        if (!changed_in_batch)
          break;
        dirty = std::move(touched);
      }
      if (!found_small_line)
        break;
    }

    water.setGeometricPropertiesForce();
    water.checkConnectivity();

    for (const auto& l : water.getLines())
      if (!l->checkTopology())
        throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "topology error detected in line after merge");

    /* --------------------------- エッジフリップによる表面メッシュの改善-------------------------- */

    if (surface_flip) {
      flipIfBatched(water, {20 * rad /*target n diff*/, 20 * rad /*change n diff*/}, {5 * rad, 5 * rad}, "post-collapse"); // Stricter for Neumann
    }

    water.setGeometricPropertiesForce();
    water.checkConnectivity();

    for (const auto& l : water.getLines())
      if (!l->checkTopology())
        throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "topology error after flip");
  }

  if (tetrahedralize)
    water.tetrahedralize();
  water.setGeometricPropertiesForce();

  water.improveTetrahedraDelaunay();
}
