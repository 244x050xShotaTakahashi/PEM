# 安息角解析機能 - PEMシミュレーション

## 概要

PEMシミュレータに安息角解析機能を追加しました。この機能により、粒子を上部から落下させて形成される山の安息角を測定し、実験値と比較することができます。

## 主な機能

1. **ホッパーからの粒子落下**: 指定された間隔で粒子を順次落下させます
2. **ホッパー形状の定義**: 傾斜した壁を持つ漏斗形状で粒子を誘導します
3. **リアルタイム角度測定**: シミュレーション中に継続的に安息角を計算します
4. **実験値との比較**: 各種材料の文献値と自動的に比較し、妥当性を検証します
5. **結果の可視化**: アニメーションと角度履歴グラフを生成します

## 使用方法

### 1. 入力ファイルの準備

`input/input_angle_of_repose.dat` ファイルで以下のパラメータを設定します：

```
# 安息角解析モード設定
ANGLE_OF_REPOSE_MODE        1           # 安息角解析モードを有効化 (1=ON, 0=OFF)
HOPPER_OPENING_WIDTH        0.02        # ホッパー開口部の幅 [m]
HOPPER_HEIGHT               0.15        # ホッパーの高さ [m]
HOPPER_ANGLE                45.0        # ホッパーの傾斜角度 [度]
DROP_HEIGHT                 0.20        # 粒子落下開始高さ [m]
TOTAL_DROP_PARTICLES        500         # 落下させる総粒子数 [-]
DROP_INTERVAL_STEPS         200         # 粒子落下間隔 [ステップ]
TARGET_ANGLE_OF_REPOSE      32.0        # 目標安息角（実験値）[度]
BASE_PLATE_WIDTH            0.3         # 底板の幅 [m]
MEASURE_ANGLE_CONTINUOUSLY  1           # 連続的に角度を測定 (1=ON, 0=OFF)
```

### 2. シミュレーションの実行

```bash
# コンパイル
cd PEM/src
ifort -O3 pem_simulator.f90 -o pem_simulator_angle

# 実行
./pem_simulator_angle ../input/input_angle_of_repose.dat
```

または、HPCシステム（Camphor）でのジョブ投入：

```bash
pjsub job_angle_of_repose.sh
```

### 3. 結果の確認

シミュレーション完了後、以下の情報が出力されます：

- **コンソール出力**: リアルタイムの安息角測定値と最終的な妥当性検証結果
- **data/angle_history.d**: 時間毎の安息角の変化
- **data/graph11.d**: 粒子位置データ（アニメーション用）

### 4. 可視化

```bash
cd PEM/src
python3 visualize_angle_of_repose.py
```

生成される出力：
- `angle_of_repose_animation.gif`: 粒子堆積のアニメーション
- `angle_history_plot.png`: 安息角の時間変化グラフ
- `final_state_angle_of_repose.png`: 最終状態と測定角度

## パラメータ調整のガイドライン

### 安息角が低すぎる場合
- 粒子間摩擦係数（FRICTION_COEFF_PARTICLE）を増加
- 壁-粒子間摩擦係数（FRICTION_COEFF_WALL）を増加

### 安息角が高すぎる場合
- 摩擦係数を減少
- 粒子サイズ分布を調整

### 計算が不安定な場合
- 時間刻み（TIME_STEP）を減少
- ヤング率を調整して接触剛性を変更

## 実験値の参考範囲

| 材料 | 安息角範囲 [度] | 典型値 [度] |
|------|----------------|------------|
| 乾燥砂 | 30-35 | 32 |
| 湿潤砂 | 35-45 | 40 |
| 砂利 | 35-40 | 37 |
| 砕石 | 38-45 | 42 |
| ガラスビーズ | 22-28 | 25 |
| 小麦 | 25-30 | 27 |
| 米 | 28-35 | 32 |

## トラブルシューティング

### 粒子が詰まる場合
- ホッパー開口幅（HOPPER_OPENING_WIDTH）を増加
- 粒子半径を減少

### 角度測定が不安定な場合
- より多くの粒子を使用（TOTAL_DROP_PARTICLES）
- 粒子落下間隔（DROP_INTERVAL_STEPS）を調整

### メモリ不足の場合
- 最大粒子数（ni_max）の制限を確認
- 出力間隔（OUTPUT_INTERVAL_NORMAL）を増加

## 理論的背景

安息角は粒状材料の流動性を特徴づける重要なパラメータです。以下の要因が影響します：

1. **粒子形状**: 角張った粒子ほど高い安息角
2. **粒子サイズ分布**: 均一なサイズより混合サイズの方が密に充填
3. **表面粗さ**: 粗い表面ほど高い摩擦
4. **水分含有量**: 適度な水分は粒子間の付着力を増加
5. **粒子密度**: 密度差は分離を引き起こす可能性

## 参考文献

1. Mehta, A., & Barker, G. C. (1994). The dynamics of sand. Reports on Progress in Physics, 57(4), 383.
2. 粉体工学会編 (2005). 粉体シミュレーション入門. 産業図書.
3. Cundall, P. A., & Strack, O. D. (1979). A discrete numerical model for granular assemblies. Geotechnique, 29(1), 47-65.

## 開発者向け情報

### 追加されたサブルーチン

- `drop_particle_sub`: 新しい粒子を落下させる
- `calculate_angle_of_repose`: 安息角を計算
- `validate_angle_of_repose`: 実験値と比較して妥当性を検証

### 修正されたサブルーチン

- `wcont_sub`: ホッパー壁との接触を追加
- `fposit_sub`: 安息角モードの初期化を追加
- `gfout_sub`: 角度履歴の出力を追加

---

**更新日**: 2025年7月  
**バージョン**: 2.0.0 