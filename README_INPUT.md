# PEMシミュレーション inputファイルシステム

## 概要

このPEMシミュレーションは、inputファイルから設定パラメータを読み込むように改良されました。これにより、コードを再コンパイルすることなく、異なる条件でシミュレーションを実行できます。

## 使用方法

### 1. 基本的な使用方法

```bash
# デフォルトのinput.datファイルを使用
./pem_simulator

# 特定のinputファイルを指定
./pem_simulator input_validation.dat
./pem_simulator input_normal.dat
```

### 2. 提供されるinputファイル

- `input.dat`: デフォルトのinputファイル（通常モード）
- `input_normal.dat`: 通常の粒状体シミュレーション用
- `input_validation.dat`: 1次元弾性衝突の検証用

## inputファイルの書式

inputファイルは以下の形式で記述します：

```
# コメント行（#で始まる行）
キーワード    値    # 説明（任意）
```

### 利用可能なキーワード

#### シミュレーション制御パラメータ
- `TIME_STEP`: 時間刻み [s]
- `MAX_CALCULATION_STEPS`: 最大計算ステップ数

#### 材料物性値
- `YOUNG_MODULUS_PARTICLE`: 粒子のヤング率 [Pa]
- `YOUNG_MODULUS_WALL`: 壁のヤング率 [Pa]
- `POISSON_RATIO_PARTICLE`: 粒子のポアソン比 [-]
- `POISSON_RATIO_WALL`: 壁のポアソン比 [-]
- `PARTICLE_DENSITY`: 粒子の密度 [kg/m³]

#### 摩擦係数
- `FRICTION_COEFF_PARTICLE`: 粒子間摩擦係数 [-]
- `FRICTION_COEFF_WALL`: 壁-粒子間摩擦係数 [-]

#### 粒子生成パラメータ（通常モード）
- `PARTICLE_RADIUS_LARGE`: 大きな粒子の半径 [m]
- `PARTICLE_RADIUS_SMALL`: 小さな粒子の半径 [m]
- `CONTAINER_WIDTH`: 容器の幅 [m]
- `PARTICLE_GEN_LAYERS`: 初期粒子生成層数 [-]
- `RANDOM_SEED`: 乱数シード [-]

#### 検証モード設定
- `VALIDATION_MODE`: 検証モードフラグ（1=ON, 0=OFF）
- `VALIDATION_PARTICLE_RADIUS`: 検証用粒子半径 [m]
- `VALIDATION_PARTICLE1_X`: 粒子1のx座標 [m]
- `VALIDATION_PARTICLE1_Z`: 粒子1のz座標 [m]
- `VALIDATION_PARTICLE2_X`: 粒子2のx座標 [m]
- `VALIDATION_PARTICLE2_Z`: 粒子2のz座標 [m]
- `VALIDATION_PARTICLE1_VX`: 粒子1の初期x速度 [m/s]
- `VALIDATION_PARTICLE2_VX`: 粒子2の初期x速度 [m/s]

#### 出力制御パラメータ
- `OUTPUT_INTERVAL_NORMAL`: 通常モード出力間隔 [ステップ]
- `OUTPUT_INTERVAL_VALIDATION`: 検証モード出力間隔 [ステップ]

## 使用例

### 通常モードでの実行

```bash
# 通常モード用のinputファイルを使用
./pem_simulator input_normal.dat
```

このモードでは粒状体の重力落下シミュレーションが実行されます。

### 検証モードでの実行

```bash
# 検証モード用のinputファイルを使用
./pem_simulator input_validation.dat
```

このモードでは1次元弾性衝突の検証が実行され、理論値との比較結果が表示されます。

## カスタムinputファイルの作成

新しいシミュレーション条件を設定するには：

1. 既存のinputファイルをコピー
2. 必要なパラメータを修正
3. 新しいファイル名で保存
4. コマンドライン引数でファイルを指定して実行

例：
```bash
cp input_normal.dat my_simulation.dat
# my_simulation.datを編集
./pem_simulator my_simulation.dat
```

## 注意事項

- inputファイルが見つからない場合、エラーメッセージが表示され、プログラムが終了します
- 不明なキーワードがある場合は警告メッセージが表示されますが、処理は続行されます
- 指定されていないパラメータにはデフォルト値が使用されます
- 検証モードでは、摩擦係数は実行時に自動的に0にリセットされます 