# BEM Project

## GUI
設定ファイル編集用GUIの場所：
```
/Users/tomoaki/Library/CloudStorage/Dropbox/code/python/GUI_for_simulation_inputfiles
```

起動コマンド：
```bash
.venv/bin/python qt_pyside_gui/main.py
```

---

## Metal M2L (GPU加速M2L変換)

### 概要

FMMのM2L (Multipole-to-Local) 変換をMetal GPUで高速化。
約4億8600万項の行列ベクトル積をGPUで並列計算する。

### ファイル構成

```
metal_m2l/
├── CMakeLists.txt           # Metal ライブラリのビルド設定
├── metal_m2l.metal          # GPU シェーダー (MSL)
├── metal_m2l.mm             # Objective-C++ Metal API 実装
├── metal_m2l.h              # C API ヘッダー
├── metal_m2l_wrapper.hpp    # C++ ラッパークラス (ヘッダーオンリー)
└── build/
    └── libmetal_m2l.a       # ビルド済み静的ライブラリ
```

### ファイル間の依存関係

```
moments.hpp
    │
    ├── M2LTerm 構造体定義
    │   - AAAY (係数)
    │   - src_MM (ソースMMポインタ)
    │   - src_bucket (バケットポインタ) ← O(1)ルックアップ用
    │   - src_coeff_idx (係数インデックス) ← O(1)ルックアップ用
    │
    └── #include "metal_m2l_wrapper.hpp" (USE_METAL_M2L 有効時)
            │
            ├── setup(): M2LTerm から GPU バッファを構築
            │   - bucket_ptr_to_idx マップで O(1) ルックアップ
            │
            └── #include "metal_m2l.h" (C API)
                    │
                    └── libmetal_m2l.a (静的ライブラリ)
                            │
                            ├── metal_m2l.mm (Objective-C++ 実装)
                            │
                            └── metal_m2l.metal (GPU カーネル)
                                - m2l_float_kahan: Float + Kahan summation
                                - m2l_float_kahan_tg: Threadgroup 並列
```

### ビルド手順

#### 1. Metal ライブラリのビルド（変更時のみ）

`metal_m2l.metal` または `metal_m2l.mm` を変更した場合:

```bash
cd /Users/tomoaki/Library/CloudStorage/Dropbox/code/cpp/builds/build_bem/metal_m2l
mkdir -p build && cd build
cmake ..
make -j8
```

生成物:
- `libmetal_m2l.a` - 静的ライブラリ
- `default.metallib` - コンパイル済みGPUシェーダー

#### 2. メインBEMプログラムのビルド

ヘッダファイル (`metal_m2l_wrapper.hpp`, `moments.hpp`) を変更した場合、
またはMetal M2Lを有効にする場合:

```bash
cd /Users/tomoaki/Library/CloudStorage/Dropbox/code/cpp/builds/build_bem
make -j8
```

CMake使用時:
```bash
cmake --build . -j8
```

### コンパイルフラグ

| フラグ | 説明 |
|--------|------|
| `USE_METAL_M2L` | Metal M2L を有効化 |

### 設定（settings.json）

Metal M2LはGUIまたはsettings.jsonから設定可能。GMRES solver使用時のみ有効。

| キー | 型 | デフォルト | 説明 |
|------|-----|-----------|------|
| `use_metal_m2l` | bool | false | Metal M2L を有効化 (Float+Kahan精度) |
| `metal_m2l_threadgroup` | bool | false | Threadgroup並列カーネルを使用 |
| `metal_m2l_sort_terms` | bool | false | 項をソートしてメモリ局所性向上 |

設定例 (settings.json):
```json
{
  "solver": "GMRES",
  "use_metal_m2l": true,
  "metal_m2l_threadgroup": false,
  "metal_m2l_sort_terms": false
}
```

### カーネル選択ガイド

| カーネル | threadgroup | 特徴 |
|---------|-------------|------|
| `m2l_float_kahan` | false (デフォルト) | 1スレッド/行、シンプル |
| `m2l_float_kahan_tg` | true | Threadgroup並列、長い行向け |

M3/M4 GPUでは `threadgroup=false` を推奨（flexible on-chip memoryにより直接バッファアクセスが効率的）。

### 最適化履歴

#### O(n×m) → O(1) ルックアップ最適化 (2025-01-28)

**問題**: setup() で各M2L項のソースバケットを特定するのに O(項数 × バケット数) の計算が必要だった。
- 約4.86億項 × 639バケット ≒ 3100億回の比較
- setup時間: 約170秒

**解決**: M2LTerm構造体に `src_bucket` と `src_coeff_idx` を追加し、
set_m2l() 時点で情報を保存。setup() では O(1) のハッシュマップ検索のみ。

```cpp
// moments.hpp - M2LTerm 拡張
struct M2LTerm {
  cd AAAY;
  const std::array<cd, 2> *src_MM;
  const void* src_bucket;   // 追加: ソースバケットポインタ
  int src_coeff_idx;        // 追加: MM_配列内インデックス
};

// metal_m2l_wrapper.hpp - O(1) ルックアップ
std::unordered_map<const void*, int32_t> bucket_ptr_to_idx;
for (int32_t i = 0; i < num_buckets_; ++i) {
    bucket_ptr_to_idx[all_buckets[i]] = i;
}
// term.src_bucket から直接インデックスを取得
auto it = bucket_ptr_to_idx.find(term.src_bucket);
int32_t src_bucket_idx = it->second;
```

### 変更時の再ビルド判定表

| 変更ファイル | Metal lib | BEM main |
|-------------|:---------:|:--------:|
| `metal_m2l.metal` | ○ | ○ |
| `metal_m2l.mm` | ○ | ○ |
| `metal_m2l.h` | △* | ○ |
| `metal_m2l_wrapper.hpp` | - | ○ |
| `moments.hpp` | - | ○ |

*: API変更時のみ

---

## M2L 精度モード

| モード | 実行場所 | 精度 | 有効化方法 |
|--------|----------|------|-----------|
| Float + Kahan | GPU (Metal) | ~7桁 + 誤差補償 | `BEM_METAL_M2L=1` |
| Double | CPU | ~15桁 | デフォルト (Metal M2L無効時) |

GPU M2Lは Float + Kahan summation で計算される。
高精度が必要な場合は `use_metal_m2l=false` でCPU double計算に切り替え可能。

---

## FMM タイミング統計の注意点

### 現状の問題

報告される near/far 時間が実際のウォールクロック時間より大きい。

**原因**: OpenMP 並列領域内で各スレッドの時間を合計しているため、
CPU時間（全スレッドの計算時間の総和）が報告される。

例: 8スレッドが並列で0.1秒ずつ計算
- ウォールクロック: 0.1秒
- 報告時間: 0.8秒 (8 × 0.1)

### 改善案（検討中）

- 並列領域全体のウォールクロックを計測し、near+far 合計として報告
- スレッド毎の計測オーバーヘッドも削減

---

## 変更履歴

### 2025-01-31: Metal M2L設定をJSON/GUIベースに移行

環境変数依存を削除し、settings.json と GUI から Metal M2L を設定可能に。

**変更内容**:
- 環境変数 `BEM_METAL_M2L`, `BEM_METAL_M2L_TG`, `BEM_METAL_M2L_SORT_TERMS` を削除
- settings.json キー追加: `use_metal_m2l`, `metal_m2l_threadgroup`, `metal_m2l_sort_terms`
- GUI に Metal M2L 設定セクション追加 (GMRES solver 選択時のみ表示)
- `MTLMathModeFast` → `MTLMathModeSafe` に変更 (Kahan summation の精度保証)
- Metal シェーダーコンパイルに `-fno-fast-math` 追加

**影響ファイル**:
- `BEM_inputfile_reader.hpp`: JSON パース追加
- `BEM_solveBVP.hpp`, `BEM_solveBVP_GMRES_FMM.hpp`: グローバル変数追加
- `main_time_domain.cpp`: 設定読み込み
- `moments.hpp`: `MetalM2LSettings` 構造体、`initMetalM2L()` 更新
- `metal_m2l.mm`: `metal_m2l_init(int threadgroup_mode)` シグネチャ変更
- `metal_m2l_wrapper.hpp`: コンストラクタ引数追加
- GUI: `settings_editor.py`, `utils.py`, `tooltips.py`

### 2025-01-31: TG=2 (SIMD shuffle カーネル) 削除

Metal M2LのTG=2モード（SIMD shuffle reduction）を削除。

**削除理由**:
- 収束しない問題が発生
- SIMD shuffle操作の条件分岐での未定義動作が原因と思われる
- TG=0/TG=1で十分な性能が得られる

**削除対象**:
- `m2l_float_kahan_tg_simd` カーネル (`metal_m2l.metal`)
- `simd_reduce_add` ヘルパー関数 (`metal_m2l.metal`)
- `pipelineFloatKahanTGSimd` パイプライン (`metal_m2l.mm`)
- TG=2モード処理 (`metal_m2l.mm`)
- `BEM_METAL_M2L_TG=2` 環境変数サポート

**残存カーネル**:
- `m2l_float_kahan` (TG=0, デフォルト): 1スレッド/行
- `m2l_float_kahan_tg` (TG=1): Threadgroup並列

### 2025-01-28: M2L Double-Float精度オプション削除

M2LのDouble-Float精度オプションを削除。
理由:
- 収束が良くない
- 精度と速度のトレードオフが中途半端

削除対象:
- `BEM_METAL_M2L_DOUBLE_FLOAT` 環境変数
- `m2l_double_float`, `m2l_double_float_full_output` カーネル
- 関連するC API, C++ラッパーコード

現在のM2L精度オプション:
- **GPU (Metal)**: Float + Kahan summation (~7桁 + 誤差補償)
- **CPU**: Full double (~15桁)

### 2025-01-28: Near-field混合精度コードの完全削除

Near-field直接積分の混合精度（float）オプションを完全に削除。

**削除理由**:
- 効果が薄い
- データ保存に時間がかかる
- トレードオフの調整が難しく、将来性が感じられない

**削除対象** (`lib_multipole_expansion.hpp`):
- `use_mixed_precision` パラメータ（`setDirectIntegration`）
- Float用メンバ変数群:
  - `near_indices_double`, `near_weights_phi_double`, `near_weights_phin_double`
  - `near_indices_float`, `near_weights_phi_float`, `near_weights_phin_float`
  - `near_run_base_idx_float`, `near_run_pos_float`, `near_run_len_float`
  - `near_weights_phi_float_high/low`, `near_weights_phin_float_high/low`
- Float用メソッド:
  - `integrateNearFieldFloat`
  - `integrateNearFieldDense`
  - `integrateNearFieldDoublePart`
  - `integrateNearFieldFloatPart`

**削除対象** (`moments.hpp`):
- `initializeFMM`の`use_mixed_precision`パラメータ

**リネーム**:
- `near_weights_phi_D` → `near_weights_phi`
- `near_weights_phin_D` → `near_weights_phin`
- `integrateNearFieldDenseDouble` → `integrateNearField`

**その他の削除対象**:
- `USE_METAL_NEARFIELD` コンパイルオプション
- 関連するMetal GPU near-fieldコード

**注意**: M2L関連のMetal GPU (`USE_METAL_M2L`) は維持。

### 2025-01-29: 座標スケーリングの改善

BIE求解時の数値安定性向上のため、座標スケーリングを改善。

**変更内容**:

1. **スケール係数の変更** (`Network.hpp`):
   - 変更前: 最大辺長 `max(Lx, Ly, Lz)`
   - 変更後: 対角距離の3/4 `sqrt(Lx² + Ly² + Lz²) * 0.75`

   ```cpp
   double computeCharacteristicLength() const {
     double diagonal = std::sqrt(Lx * Lx + Ly * Ly + Lz * Lz);
     return diagonal * 0.75;
   }
   ```

2. **デフォルト有効化** (`BEM_solveBVP.hpp`, `BEM_solveBVP_GMRES_FMM.hpp`):
   - スケーリングは常にON（環境変数チェックを削除）
   - `coordinate_scale_factor_ > 1e-10` の場合のみ適用

**スケーリングの数学**:

座標変換 `x' = x/L` の場合：
- `φ' = φ` (ポテンシャルは不変)
- `∂φ/∂n' = L × ∂φ/∂n` (法線微分はスケール)

| 段階 | Neumann BC (phin既知) | Dirichlet BC (phi既知) |
|------|----------------------|----------------------|
| 入力時 | `phin' = phin × L` | `phi' = phi` |
| 出力時 | `phi'` そのまま | `phin = phin' / L` |

---

## FMM 現状まとめ

### ツリー再利用 (reuse_static_tree_cache)

`reuse_static_tree_cache=true` の場合でも、以下は毎回実行される：

| ステップ | 処理 | 再利用時 |
|---------|------|---------|
| Step 1-4 | ツリー構築、M2M/M2L/L2L設定 | スキップ |
| Step 5 | setL2P | **毎回実行** |
| Step 6 | setDirectIntegration | **毎回実行** (~0.37s) |

**理由**: ALE等で節点位置が動くと、G/Gnカーネル（距離依存）の再計算が必要。
インデックス構造 (`near_indices`, `near_run_*`) は不変だが、重み (`near_weights_*`) は変わる。

### setDirectIntegration の構造

```
setDirectIntegration (lib_multipole_expansion.hpp)
├── 各ソースに対してコールバック呼び出し
│   └── use_this_source_when_set_direct_integration(this)
│       └── 25点Gauss求積でG/Gnカーネル計算 ← 主なコスト
├── unordered_map に集約
├── ソートして連続配列に変換
└── RLE (Run-Length Encoding) 構築
```

### Metal M2L の row_info_ について

`reuse_static_fmm=true` の場合、バケットのメモリアドレスは固定されるため、
`row_info_` 内の `dst_mm` ポインタは有効なまま。再設定不要。

---

## 次のステップ

- [x] Metal M2L O(1)ルックアップ最適化
- [x] 座標スケーリングの改善（対角距離の3/4）
- [ ] FMM タイミング統計の修正（ウォールクロック計測）
- [ ] より大規模な問題での検証
