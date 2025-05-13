#ifndef fmm_h
#define fmm_h

#include "lib_multipole_expansion.hpp"
#include "lib_spatial_partitioning.hpp"

struct FMM {
  public:
   FMM(std::unique_ptr<Network> obj, int max_level)
       : obj(std::move(obj)),
         max_level(max_level),
         B_poles(this->obj->scaledBounds(1.1), this->obj->getScale() / 5.0),
         initialized(false) {}

   void initialize() {
      // ... (省略)
   }

   void generateTree() {
      // ... (省略)
   }

   void multipoleExpansion(bool reuse = false) {
      // ... (省略)
   }

   void M2M_M2L_L2L_L2P() {
      TimeWatch tw;
#if defined(_DEBUG_FMM_)
      std::cout << Magenta << "M2M ..." << colorReset << std::endl;
#endif
      M2M(B_poles);
#if defined(_DEBUG_FMM_)
      std::cout << Magenta << "M2M" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;
      std::cout << Magenta << "M2L ..." << colorReset << std::endl;
#endif
      M2L(B_poles);
#if defined(_DEBUG_FMM_)
      std::cout << Magenta << "M2L" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;
      std::cout << Magenta << "L2L ..." << colorReset << std::endl;
#endif
      L2L(B_poles);
#if defined(_DEBUG_FMM_)
      std::cout << Magenta << "L2L" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;
      std::cout << Magenta << "L2P ..." << colorReset << std::endl;
#endif
#pragma omp parallel
      for (auto& p : obj->getPoints())
#pragma omp single nowait
      {
         auto [wGPhin_wGnPhi_near, wGPhin_wGnPhi_far] = L2P(B_poles, p->X);
         p->wGPhin_wGnPhi_near = wGPhin_wGnPhi_near;
         p->wGPhin_wGnPhi_far = wGPhin_wGnPhi_far;
         p->wGPhin_wGnPhi_FMM = wGPhin_wGnPhi_near + wGPhin_wGnPhi_far;
      }
#if defined(_DEBUG_FMM_)
      std::cout << Magenta << "L2P" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;
#endif
   }

   void directIntegration() {
      // ... (省略)
   }

   void outputResults() {
      // ... (省略)
   }

   // 行列ベクトル積を計算するメソッド
   V_d apply(const V_d& x) {
      // xをFMMの形式に変換
      // ここでは簡略化のために直接使用
      for (auto& p : obj->getPoints()) {
         p->phiphin = {x[p->index], x[p->index]};
      }

      // FMMを実行
      this->multipoleExpansion();
      this->M2M_M2L_L2L_L2P();

      // 結果を取得
      V_d result(x.size());
      for (const auto& p : obj->getPoints()) {
         result[p->index] = p->wGPhin_wGnPhi_FMM[0];  // 適切な結果を設定
      }

      return result;
   }

  private:
   std::unique_ptr<Network> obj;
   int max_level;
   Buckets<sp_pole4FMM> B_poles;
   bool initialized;
};

#endif