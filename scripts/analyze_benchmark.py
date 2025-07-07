#!/usr/bin/env python3
"""
PEM Algorithm Benchmark Result Analysis
セル法アルゴリズムの効果を分析・可視化する
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import sys

def analyze_benchmark_results(csv_file="benchmark_results.csv"):
    """
    ベンチマーク結果を分析し、グラフを作成する
    """
    
    # CSVファイルの読み込み
    try:
        df = pd.read_csv(csv_file)
        print(f"ベンチマーク結果を読み込みました: {csv_file}")
        print(f"データポイント数: {len(df)}")
    except FileNotFoundError:
        print(f"エラー: ファイル '{csv_file}' が見つかりません")
        return
    except pd.errors.EmptyDataError:
        print(f"エラー: ファイル '{csv_file}' が空です")
        return
    except pd.errors.ParserError as e:
        print(f"CSVファイル解析エラー: {e}")
        return
    except Exception as e:
        print(f"ファイル読み込みエラー: {e}")
        return
    
    # データの基本検証
    required_columns = ['Particle_Layers', 'Particles', 'Algorithm', 'Time_seconds', 'Steps', 'Time_per_step']
    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        print(f"エラー: 必要な列が不足しています: {missing_columns}")
        return
    
    # データ型の検証
    try:
        df['Particle_Layers'] = pd.to_numeric(df['Particle_Layers'], errors='coerce')
        df['Particles'] = pd.to_numeric(df['Particles'], errors='coerce')
        df['Time_seconds'] = pd.to_numeric(df['Time_seconds'], errors='coerce')
        df['Steps'] = pd.to_numeric(df['Steps'], errors='coerce')
        df['Time_per_step'] = pd.to_numeric(df['Time_per_step'], errors='coerce')
    except Exception as e:
        print(f"データ型変換エラー: {e}")
        return
    
    # NaN値のチェック
    if df[required_columns].isnull().any().any():
        print("警告: データに無効な値（NaN）が含まれています")
        print("無効な行:")
        print(df[df[required_columns].isnull().any(axis=1)])
        # 無効な行を除去
        df = df.dropna(subset=required_columns)
        print(f"有効なデータ行数: {len(df)}")
    
    if len(df) == 0:
        print("エラー: 有効なデータがありません")
        return
    
    # データの表示
    print("\n=== ベンチマーク結果 ===")
    print(df.to_string(index=False))
    
    # セル法ON/OFFでデータを分離
    df_on = df[df['Algorithm'] == 'ON']
    df_off = df[df['Algorithm'] == 'OFF']
    
    if len(df_on) == 0:
        print("警告: セル法ONのデータがありません")
        return
    if len(df_off) == 0:
        print("警告: セル法OFFのデータがありません")
        print("セル法ONのデータのみで表示します")
        # セル法ONのみのデータを表示
        create_single_algorithm_plots(df_on)
        return
    
    # マージして比較用データフレームを作成
    merged = pd.merge(df_on, df_off, on='Particle_Layers', suffixes=('_ON', '_OFF'))
    
    if len(merged) == 0:
        print("警告: 対応するON/OFFデータが見つかりません")
        print("個別にデータを表示します")
        create_separate_plots(df_on, df_off)
        return
    
    # 加速比の計算
    merged['Speedup'] = merged['Time_seconds_OFF'] / merged['Time_seconds_ON']
    merged['Efficiency_Ratio'] = merged['Time_per_step_OFF'] / merged['Time_per_step_ON']
    
    print("\n=== 加速比分析 ===")
    for _, row in merged.iterrows():
        print(f"Layers: {row['Particle_Layers']}, "
              f"Particles: {row['Particles_ON']}, "
              f"Speedup: {row['Speedup']:.2f}x, "
              f"Time ON: {row['Time_seconds_ON']:.2f}s, "
              f"Time OFF: {row['Time_seconds_OFF']:.2f}s")
    
    # 可視化
    create_benchmark_plots(merged, df_on, df_off)
    
    # 統計情報の出力
    print_statistics(merged)

def create_benchmark_plots(merged, df_on, df_off):
    """
    ベンチマーク結果のグラフを作成する
    """
    
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle('PEM Algorithm Benchmark Results', fontsize=16)
    
    # 1. 実行時間の比較
    ax1 = axes[0, 0]
    ax1.plot(merged['Particles_ON'], merged['Time_seconds_ON'], 'o-', label='Cell Algorithm ON', linewidth=2)
    ax1.plot(merged['Particles_OFF'], merged['Time_seconds_OFF'], 's-', label='Cell Algorithm OFF', linewidth=2)
    ax1.set_xlabel('Number of Particles')
    ax1.set_ylabel('Execution Time (seconds)')
    ax1.set_title('Execution Time Comparison')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.set_yscale('log')
    
    # 2. 加速比
    ax2 = axes[0, 1]
    ax2.plot(merged['Particles_ON'], merged['Speedup'], 'ro-', linewidth=2, markersize=8)
    ax2.set_xlabel('Number of Particles')
    ax2.set_ylabel('Speedup Ratio')
    ax2.set_title('Speedup Ratio (OFF/ON)')
    ax2.grid(True, alpha=0.3)
    
    # 理論線 (O(N^2) vs O(N))
    particles = merged['Particles_ON']
    theoretical_speedup = (particles / particles.min()) ** 1.5  # 近似的な理論値
    ax2.plot(particles, theoretical_speedup, '--', alpha=0.7, label='Theoretical O(N^1.5)')
    ax2.legend()
    
    # 3. 1ステップあたりの時間
    ax3 = axes[1, 0]
    ax3.plot(merged['Particles_ON'], merged['Time_per_step_ON'], 'o-', label='Cell Algorithm ON', linewidth=2)
    ax3.plot(merged['Particles_OFF'], merged['Time_per_step_OFF'], 's-', label='Cell Algorithm OFF', linewidth=2)
    ax3.set_xlabel('Number of Particles')
    ax3.set_ylabel('Time per Step (seconds)')
    ax3.set_title('Time per Step Comparison')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    ax3.set_yscale('log')
    
    # 4. 計算複雑度の分析
    ax4 = axes[1, 1]
    particles = merged['Particles_ON']
    
    # 正規化された時間 (最小値で正規化)
    norm_time_on = merged['Time_per_step_ON'] / merged['Time_per_step_ON'].min()
    norm_time_off = merged['Time_per_step_OFF'] / merged['Time_per_step_OFF'].min()
    
    ax4.loglog(particles, norm_time_on, 'o-', label='Cell Algorithm ON', linewidth=2)
    ax4.loglog(particles, norm_time_off, 's-', label='Cell Algorithm OFF', linewidth=2)
    
    # 理論線
    n_min = particles.min()
    linear = particles / n_min
    quadratic = (particles / n_min) ** 2
    
    ax4.loglog(particles, linear, '--', alpha=0.7, label='O(N) Linear')
    ax4.loglog(particles, quadratic, '--', alpha=0.7, label='O(N²) Quadratic')
    
    ax4.set_xlabel('Number of Particles')
    ax4.set_ylabel('Normalized Time per Step')
    ax4.set_title('Computational Complexity Analysis')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('benchmark_results.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    print("\nグラフを保存しました: benchmark_results.png")

def print_statistics(merged):
    """
    統計情報を出力する
    """
    print("\n=== 統計情報 ===")
    print(f"平均加速比: {merged['Speedup'].mean():.2f}x")
    print(f"最大加速比: {merged['Speedup'].max():.2f}x")
    print(f"最小加速比: {merged['Speedup'].min():.2f}x")
    
    # 最大粒子数での効果
    max_particles_row = merged.loc[merged['Particles_ON'].idxmax()]
    print(f"\n最大粒子数({max_particles_row['Particles_ON']})での効果:")
    print(f"  セル法ON: {max_particles_row['Time_seconds_ON']:.2f}秒")
    print(f"  セル法OFF: {max_particles_row['Time_seconds_OFF']:.2f}秒")
    print(f"  加速比: {max_particles_row['Speedup']:.2f}x")
    
    # 計算複雑度の推定
    particles = merged['Particles_ON'].values
    time_on = merged['Time_per_step_ON'].values
    time_off = merged['Time_per_step_OFF'].values
    
    # 線形回帰で計算複雑度を推定
    log_particles = np.log(particles)
    log_time_on = np.log(time_on)
    log_time_off = np.log(time_off)
    
    slope_on = np.polyfit(log_particles, log_time_on, 1)[0]
    slope_off = np.polyfit(log_particles, log_time_off, 1)[0]
    
    print(f"\n計算複雑度の推定:")
    print(f"  セル法ON: O(N^{slope_on:.2f})")
    print(f"  セル法OFF: O(N^{slope_off:.2f})")

def create_single_algorithm_plots(df):
    """
    単一アルゴリズム（セル法ONのみ）のグラフを作成する
    """
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle('PEM Single Algorithm Results (Cell Algorithm ON)', fontsize=16)
    
    # 1. 実行時間
    ax1 = axes[0]
    ax1.plot(df['Particles'], df['Time_seconds'], 'o-', label='Cell Algorithm ON', linewidth=2, color='blue')
    ax1.set_xlabel('Number of Particles')
    ax1.set_ylabel('Execution Time (seconds)')
    ax1.set_title('Execution Time')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # 2. 1ステップあたりの時間
    ax2 = axes[1]
    ax2.plot(df['Particles'], df['Time_per_step'], 'o-', label='Cell Algorithm ON', linewidth=2, color='blue')
    ax2.set_xlabel('Number of Particles')
    ax2.set_ylabel('Time per Step (seconds)')
    ax2.set_title('Time per Step')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('single_algorithm_results.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    print("\nグラフを保存しました: single_algorithm_results.png")

def create_separate_plots(df_on, df_off):
    """
    対応しないON/OFFデータを個別に表示する
    """
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    fig.suptitle('PEM Separate Algorithm Results', fontsize=16)
    
    # セル法ON
    ax1 = axes[0, 0]
    ax1.plot(df_on['Particles'], df_on['Time_seconds'], 'o-', label='Cell Algorithm ON', linewidth=2, color='blue')
    ax1.set_xlabel('Number of Particles')
    ax1.set_ylabel('Execution Time (seconds)')
    ax1.set_title('Cell Algorithm ON - Execution Time')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    ax2 = axes[0, 1]
    ax2.plot(df_on['Particles'], df_on['Time_per_step'], 'o-', label='Cell Algorithm ON', linewidth=2, color='blue')
    ax2.set_xlabel('Number of Particles')
    ax2.set_ylabel('Time per Step (seconds)')
    ax2.set_title('Cell Algorithm ON - Time per Step')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # セル法OFF
    ax3 = axes[1, 0]
    ax3.plot(df_off['Particles'], df_off['Time_seconds'], 's-', label='Cell Algorithm OFF', linewidth=2, color='red')
    ax3.set_xlabel('Number of Particles')
    ax3.set_ylabel('Execution Time (seconds)')
    ax3.set_title('Cell Algorithm OFF - Execution Time')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    ax4 = axes[1, 1]
    ax4.plot(df_off['Particles'], df_off['Time_per_step'], 's-', label='Cell Algorithm OFF', linewidth=2, color='red')
    ax4.set_xlabel('Number of Particles')
    ax4.set_ylabel('Time per Step (seconds)')
    ax4.set_title('Cell Algorithm OFF - Time per Step')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('separate_algorithm_results.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    print("\nグラフを保存しました: separate_algorithm_results.png")

def main():
    """
    メイン関数
    """
    # コマンドライン引数でCSVファイルを指定可能
    csv_file = sys.argv[1] if len(sys.argv) > 1 else "benchmark_results.csv"
    
    print("=== PEM Algorithm Benchmark Analysis ===")
    print(f"使用するCSVファイル: {csv_file}")
    
    analyze_benchmark_results(csv_file)

if __name__ == "__main__":
    main() 