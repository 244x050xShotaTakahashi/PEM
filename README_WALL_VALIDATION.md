# 壁-粒子間検証機能

粒子要素法（PEM）シミュレーターの壁-粒子間接触の精度検証機能について説明します。

## 概要

壁-粒子間検証機能は、シミュレーションにおける壁と粒子の相互作用が理論値と一致することを確認するための機能です。この機能により、接触力モデルの実装や数値計算の精度を評価できます。

## 検証タイプ

### 1. 自由落下反発検証 (Free Fall Rebound Validation)

粒子を指定した高さから自由落下させ、床との衝突後の反発高さを理論値と比較します。

**理論式:**
```
反発高さ = e² × 初期落下高さ
```
ここで、e は反発係数です。

**特徴:**
- エネルギー保存則の検証
- 反発係数の実装精度確認
- 時間刻みの適切性評価

### 2. 摩擦斜面検証 (Friction Slope Validation)

粒子を摩擦のある斜面上に配置し、滑り運動の理論値と比較します。

**理論式:**
```
加速度 = g × (sin(θ) - μ × cos(θ))
```
ここで、g は重力加速度、θ は斜面角度、μ は摩擦係数です。

**特徴:**
- 摩擦力モデルの検証
- 斜面接触の実装確認
- 動的摩擦の精度評価

### 3. パラメータ掃引検証 (Parameter Sweep Validation)

反発係数を変化させながら複数回の自由落下検証を実行し、統計的な精度評価を行います。

**特徴:**
- 複数パラメータでの統計評価
- 数値誤差の傾向分析
- システム全体の安定性確認

## 使用方法

### 1. 基本的な実行方法

```bash
# 自動実行スクリプトを使用（推奨）
./scripts/run_wall_validation.sh all

# 個別実行
./scripts/run_wall_validation.sh freefall
./scripts/run_wall_validation.sh slope
./scripts/run_wall_validation.sh sweep
```

### 2. 手動実行

```bash
# 自由落下検証
./pem_simulator input/input_wall_validation.dat

# 摩擦斜面検証
./pem_simulator input/input_slope_validation.dat

# パラメータ掃引検証
./pem_simulator input/input_wall_validation_sweep.dat
```

### 3. カスタム設定での実行

```bash
# テンプレートをコピーして編集
cp input/input_wall_validation_template.dat input/my_validation.dat
# my_validation.datを編集
./pem_simulator input/my_validation.dat
```

## 設定パラメータ

### 基本設定

| パラメータ | 説明 | デフォルト値 |
|-----------|------|-------------|
| `WALL_VALIDATION_MODE` | 壁検証モードの有効化 | 1 |
| `WALL_VALIDATION_TYPE` | 検証タイプ (1:自由落下, 2:斜面) | 1 |
| `TIME_STEP` | 時間刻み | 5.0e-7 |
| `MAX_CALCULATION_STEPS` | 最大計算ステップ数 | 500000 |

### 自由落下検証パラメータ

| パラメータ | 説明 | デフォルト値 |
|-----------|------|-------------|
| `WALL_VALIDATION_DROP_HEIGHT` | 初期落下高さ | 10.0 |
| `WALL_VALIDATION_PARTICLE_RADIUS` | 粒子半径 | 1.0 |
| `WALL_VALIDATION_RESTITUTION_COEFF` | 反発係数 | 0.8 |
| `WALL_VALIDATION_PARTICLE_X` | 粒子のx座標 | 5.0 |

### 摩擦斜面検証パラメータ

| パラメータ | 説明 | デフォルト値 |
|-----------|------|-------------|
| `WALL_VALIDATION_SLOPE_ANGLE` | 斜面角度（ラジアン） | 0.5236 (30度) |
| `WALL_VALIDATION_SLOPE_FRICTION` | 斜面摩擦係数 | 0.3 |
| `WALL_VALIDATION_FRICTION_COEFF` | 粒子摩擦係数 | 0.3 |

### パラメータ掃引設定

| パラメータ | 説明 | デフォルト値 |
|-----------|------|-------------|
| `PARAMETER_SWEEP_MODE` | パラメータ掃引モードの有効化 | 0 |
| `PARAMETER_SWEEP_COUNT` | 掃引回数 | 5 |
| `RESTITUTION_COEFF_MIN` | 反発係数の最小値 | 0.2 |
| `RESTITUTION_COEFF_MAX` | 反発係数の最大値 | 0.9 |

## 出力ファイル

### 検証レポート

- `data/wall_validation_freefall_report.txt` - 自由落下検証の詳細レポート
- `data/wall_validation_slope_report.txt` - 摩擦斜面検証の詳細レポート
- `data/wall_validation_summary.txt` - 総合分析レポート

### データファイル

- `data/parameter_sweep_results.csv` - パラメータ掃引結果
- `data/graph11.d`, `data/graph21.d` - 可視化用データ
- `data/backl.d` - バックアップデータ

### 分析結果

- `data/validation_analysis_results.json` - 分析結果（JSON形式）
- `data/plots/freefall_validation_analysis.png` - 可視化グラフ

## 結果分析

### 自動分析機能

検証結果は自動的に分析され、以下の統計値が計算されます：

- **平均誤差**: 全データ点での相対誤差の平均
- **最大誤差**: 最も大きな相対誤差
- **最小誤差**: 最も小さな相対誤差
- **標準偏差**: 誤差のばらつき

### Python分析スクリプト

```bash
python3 scripts/wall_validation_analyzer.py
```

このスクリプトは以下の機能を提供します：

- 検証結果の統計分析
- 可視化グラフの生成
- 総合評価レポートの作成
- JSON形式での結果出力

### 精度評価基準

| 平均誤差 | 評価 | 説明 |
|---------|-----|------|
| < 1% | 高精度 | 非常に良好な実装 |
| 1-5% | 良好 | 実用的な精度 |
| > 5% | 要改善 | パラメータ調整が必要 |

## トラブルシューティング

### よくある問題

1. **数値不安定性**
   - 時間刻みを小さくする (`TIME_STEP` を調整)
   - 材料パラメータの妥当性を確認

2. **過大な誤差**
   - 接触力モデルパラメータの確認
   - セル法アルゴリズムの設定確認

3. **コンパイルエラー**
   - Fortranコンパイラの確認
   - 依存関係の確認

### デバッグ方法

1. **詳細ログの有効化**
   ```bash
   # デバッグモードでの実行（詳細な出力）
   ./pem_simulator input/input_wall_validation.dat > debug.log 2>&1
   ```

2. **短時間実行での確認**
   ```
   MAX_CALCULATION_STEPS 1000  # 短時間実行
   OUTPUT_INTERVAL_VALIDATION 1  # 詳細出力
   ```

3. **可視化による確認**
   ```bash
   python3 src/visualize_results.py
   ```

## 理論的背景

### Hertz接触理論

壁-粒子間の法線接触力は、Hertz接触理論に基づいて計算されます：

```
Fn = Kn × δ^(3/2)
```

ここで、`Kn`は法線剛性、`δ`は重なり量です。

### 粘性ダンピング

反発係数を再現するため、粘性ダンピングが適用されます：

```
Fn_total = Fn_elastic + ηn × vn
```

ここで、`ηn`は粘性係数、`vn`は法線方向相対速度です。

### クーロン摩擦

せん断方向の力は、クーロン摩擦法則に従います：

```
|Fs| ≤ μ × |Fn|
```

ここで、`μ`は摩擦係数です。

## 参考文献

1. Cundall, P.A. and Strack, O.D.L. (1979). "A discrete numerical model for granular assemblies"
2. Hertz, H. (1882). "Über die Berührung fester elastischer Körper"
3. 粉体工学会編 (2014). "粉体の数値シミュレーション"

## 関連ファイル

- `src/pem_simulator.f90` - メインシミュレータ
- `input/input_wall_validation_template.dat` - 設定テンプレート
- `scripts/run_wall_validation.sh` - 実行スクリプト
- `scripts/wall_validation_analyzer.py` - 分析スクリプト
