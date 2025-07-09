# PEMプログラム粒子-壁間検証機能

このドキュメントでは、Qiitaの記事「DEM解析の妥当性を検証する方法」に基づいて実装された粒子-壁間計算妥当性検証機能について説明します。

## 概要

本検証機能は以下の理論解と数値解の比較を通じて、PEMプログラムの粒子-壁間相互作用の妥当性を検証します：

1. **自由落下からの反発** - エネルギー保存に基づく反発高さの検証
2. **摩擦がある斜面を滑る物体** - 運動方程式に基づく加速度・速度の検証

## 検証タイプ1: 自由落下からの反発

### 理論的背景
- 反発係数 e を持つ粒子が高さ h₀ から自由落下
- 理論的反発高さ: h₁ = e² × h₀

### 使用方法
```bash
./pem_simulator input_wall_validation.dat
```

### パラメータ設定
- `WALL_VALIDATION_MODE`: 1 (粒子-壁間検証を有効化)
- `WALL_VALIDATION_TYPE`: 1 (自由落下反発検証)
- `WALL_VALIDATION_DROP_HEIGHT`: 落下開始高さ [m]
- `WALL_VALIDATION_RESTITUTION_COEFF`: 反発係数
- `WALL_VALIDATION_PARTICLE_RADIUS`: 粒子半径 [m]

### 期待される結果
- 反発高さの理論値と計算値の比較
- エネルギー保存の確認

## 検証タイプ2: 摩擦がある斜面を滑る物体

### 理論的背景
- 斜面角度 θ、摩擦係数 μ の斜面を滑る粒子
- 理論加速度: a = g(sin θ - μ cos θ)
- 理論速度: v = at
- 理論位置: x = x₀ + ½at²

### 使用方法
```bash
./pem_simulator input_slope_validation.dat
```

### パラメータ設定
- `WALL_VALIDATION_MODE`: 1 (粒子-壁間検証を有効化)
- `WALL_VALIDATION_TYPE`: 2 (摩擦斜面検証)
- `WALL_VALIDATION_SLOPE_ANGLE`: 斜面角度 [ラジアン]
- `WALL_VALIDATION_SLOPE_FRICTION`: 摩擦係数
- `WALL_VALIDATION_PARTICLE_RADIUS`: 粒子半径 [m]

### 期待される結果
- 速度の理論値と計算値の比較
- 位置の理論値と計算値の比較

## 実装済み機能

### 入力ファイル対応
- `input_wall_validation.dat`: 自由落下反発検証用
- `input_slope_validation.dat`: 摩擦斜面検証用

### プログラム機能
- 検証モード自動検出
- 理論値自動計算
- 誤差率表示
- 詳細な計算過程出力

## 注意事項

### 現在の制限
1. **斜面検証**: 現在の実装では平坦な床面のみ対応。実際の斜面壁は未実装のため、摩擦斜面検証では理論値と大きな差が生じます。
2. **摩擦モデル**: 完全な3次元摩擦モデルが必要な場合、追加の実装が必要です。

### 拡張可能性
- 斜面壁の実装による真の摩擦斜面検証
- 回転エネルギーを含むエネルギー保存検証
- 複数粒子による集合的な検証

## ファイル構成

```
PEM/
├── src/
│   └── pem_simulator.f90          # メインプログラム（検証機能付き）
├── input_wall_validation.dat      # 自由落下反発検証用入力
├── input_slope_validation.dat     # 摩擦斜面検証用入力
└── VALIDATION_README.md           # このファイル
```

## 参考文献

- Qiita記事: 「DEM解析の妥当性を検証する方法」
- https://qiita.com/fujitagodai4/items/172ea4a5a056fc90057e

## コンパイルと実行

```bash
# コンパイル
ifort -o pem_simulator src/pem_simulator.f90

# 自由落下反発検証
./pem_simulator input_wall_validation.dat

# 摩擦斜面検証
./pem_simulator input_slope_validation.dat
``` 