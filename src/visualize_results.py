#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
改良版自由落下シミュレーション結果可視化スクリプト
dem_calc.py準拠の詳細解析機能

このスクリプトは以下の機能を提供します：
1. 粒子の軌跡アニメーション
2. 速度・位置の時系列プロット
3. エネルギー保存則の検証プロット
4. 理論解との比較プロット
5. 接触解析結果の可視化
"""

import os
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import pandas as pd
from pathlib import Path
from typing import Optional, Tuple, List

# 日本語フォント設定
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 10

class FreeFallAnalyzer:
    """自由落下シミュレーション結果解析クラス"""
    
    def __init__(self, data_dir='data'):
        """
        初期化
        
        Args:
            data_dir (str): データディレクトリのパス
        """
        self.data_dir = Path(data_dir)
        self.state_files: List[str] = []
        self.detailed_data: Optional[pd.DataFrame] = None
        self.contact_data: Optional[pd.DataFrame] = None
        self.theory_data: Optional[pd.DataFrame] = None
        self.summary_data: Optional[pd.DataFrame] = None
        
        print("========================================")
        print("  自由落下シミュレーション結果解析")
        print("  (dem_calc.py準拠の詳細解析機能)")
        print("========================================")
        
        self.load_data()
    
    def load_data(self):
        """データファイルの読み込み"""
        try:
            # 状態ファイルの読み込み
            state_pattern = str(self.data_dir / 'state_step_*.csv')
            self.state_files = sorted(glob.glob(state_pattern))
            print(f"状態ファイル数: {len(self.state_files)}")
            
            # 詳細解析データの読み込み
            energy_file = self.data_dir / 'detailed' / 'energy_analysis.csv'
            if energy_file.exists():
                data = pd.read_csv(energy_file)
                assert isinstance(data, pd.DataFrame), "エネルギー解析データの読み込みに失敗"
                self.detailed_data = data
                print(f"エネルギー解析データ: {len(self.detailed_data)} 点")
            
            # 接触解析データの読み込み
            contact_file = self.data_dir / 'analysis' / 'contact_analysis.csv'
            if contact_file.exists():
                data = pd.read_csv(contact_file)
                assert isinstance(data, pd.DataFrame), "接触解析データの読み込みに失敗"
                self.contact_data = data
                print(f"接触解析データ: {len(self.contact_data)} 回")
            
            # 理論解比較データの読み込み
            theory_file = self.data_dir / 'theory_comparison.csv'
            if theory_file.exists():
                data = pd.read_csv(theory_file)
                assert isinstance(data, pd.DataFrame), "理論解比較データの読み込みに失敗"
                self.theory_data = data
                print("理論解比較データを読み込みました")
            
            # シミュレーション概要データの読み込み
            summary_file = self.data_dir / 'simulation_summary.csv'
            if summary_file.exists():
                data = pd.read_csv(summary_file)
                assert isinstance(data, pd.DataFrame), "シミュレーション概要データの読み込みに失敗"
                self.summary_data = data
                print("シミュレーション概要データを読み込みました")
            
        except Exception as e:
            print(f"データ読み込みエラー: {e}")
            sys.exit(1)
    
    def extract_trajectory_data(self) -> Tuple[Optional[np.ndarray], Optional[np.ndarray], Optional[np.ndarray], Optional[np.ndarray]]:
        """状態ファイルから軌跡データを抽出"""
        if not self.state_files:
            print("警告: 状態ファイルが見つかりません")
            return None, None, None, None
        
        times = []
        x_positions = []
        z_positions = []
        z_velocities = []
        
        for file_path in self.state_files:
            try:
                df = pd.read_csv(file_path)
                assert isinstance(df, pd.DataFrame), f"ファイル {file_path} の読み込みに失敗"
                if len(df) > 0:
                    times.append(float(df['Time'].iloc[0]))  # type: ignore
                    x_positions.append(float(df['X'].iloc[0]))  # type: ignore
                    z_positions.append(float(df['Z'].iloc[0]))  # type: ignore
                    z_velocities.append(float(df['Vz'].iloc[0]))  # type: ignore
            except Exception as e:
                print(f"ファイル読み込みエラー {file_path}: {e}")
                continue
        
        if not times:
            return None, None, None, None
        
        return np.array(times), np.array(x_positions), np.array(z_positions), np.array(z_velocities)
    
    def plot_trajectory(self):
        """軌跡プロット"""
        times, x_pos, z_pos, z_vel = self.extract_trajectory_data()
        
        if times is None or x_pos is None or z_pos is None or z_vel is None:
            print("軌跡データが不足しています")
            return
        
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 10))
        
        # 高さの時系列
        ax1.plot(times, z_pos, 'b-', linewidth=2, label='計算値')
        ax1.set_xlabel('時間 [s]')
        ax1.set_ylabel('高さ [m]')
        ax1.set_title('粒子の高さ変化')
        ax1.grid(True, alpha=0.3)
        ax1.legend()
        
        # 速度の時系列
        ax2.plot(times, z_vel, 'r-', linewidth=2, label='z方向速度')
        ax2.set_xlabel('時間 [s]')
        ax2.set_ylabel('速度 [m/s]')
        ax2.set_title('粒子の速度変化')
        ax2.grid(True, alpha=0.3)
        ax2.legend()
        
        # 位相平面プロット
        ax3.plot(z_pos, z_vel, 'g-', linewidth=2)
        ax3.set_xlabel('高さ [m]')
        ax3.set_ylabel('速度 [m/s]')
        ax3.set_title('位相平面プロット')
        ax3.grid(True, alpha=0.3)
        
        # 2D軌跡プロット
        ax4.plot(x_pos, z_pos, 'purple', linewidth=2, marker='o', markersize=2)
        ax4.set_xlabel('x座標 [m]')
        ax4.set_ylabel('z座標 [m]')
        ax4.set_title('2D軌跡')
        ax4.grid(True, alpha=0.3)
        ax4.set_aspect('equal')
        
        plt.tight_layout()
        plt.savefig(self.data_dir / 'trajectory_analysis.png', dpi=300, bbox_inches='tight')
        plt.show()
        print("軌跡解析図を保存しました: trajectory_analysis.png")
    
    def plot_energy_analysis(self):
        """エネルギー解析プロット"""
        if self.detailed_data is None:
            print("エネルギー解析データがありません")
            return
        
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 10))
        
        # エネルギー時系列
        ax1.plot(self.detailed_data['Time_s'], self.detailed_data['Potential_Energy_J'], 
                'b-', linewidth=2, label='位置エネルギー')
        ax1.plot(self.detailed_data['Time_s'], self.detailed_data['Kinetic_Energy_J'], 
                'r-', linewidth=2, label='運動エネルギー')
        ax1.plot(self.detailed_data['Time_s'], self.detailed_data['Total_Energy_J'], 
                'g-', linewidth=2, label='全エネルギー')
        ax1.set_xlabel('時間 [s]')
        ax1.set_ylabel('エネルギー [J]')
        ax1.set_title('エネルギー保存則の検証')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # エネルギー誤差
        ax2.plot(self.detailed_data['Time_s'], self.detailed_data['Energy_Error_J'], 
                'k-', linewidth=2)
        ax2.set_xlabel('時間 [s]')
        ax2.set_ylabel('エネルギー誤差 [J]')
        ax2.set_title('エネルギー保存誤差')
        ax2.grid(True, alpha=0.3)
        
        # エネルギー比率
        total_energy_initial = self.detailed_data['Total_Energy_J'].iloc[0]  # type: ignore
        potential_ratio = self.detailed_data['Potential_Energy_J'] / total_energy_initial
        kinetic_ratio = self.detailed_data['Kinetic_Energy_J'] / total_energy_initial
        
        ax3.plot(self.detailed_data['Time_s'], potential_ratio, 'b-', linewidth=2, label='位置エネルギー比')
        ax3.plot(self.detailed_data['Time_s'], kinetic_ratio, 'r-', linewidth=2, label='運動エネルギー比')
        ax3.set_xlabel('時間 [s]')
        ax3.set_ylabel('エネルギー比 [-]')
        ax3.set_title('エネルギー比率変化')
        ax3.legend()
        ax3.grid(True, alpha=0.3)
        
        # エネルギー散布図
        ax4.scatter(self.detailed_data['Potential_Energy_J'], self.detailed_data['Kinetic_Energy_J'], 
                   c=self.detailed_data['Time_s'], cmap='viridis', alpha=0.6)
        ax4.set_xlabel('位置エネルギー [J]')
        ax4.set_ylabel('運動エネルギー [J]')
        ax4.set_title('エネルギー相関')
        cb = plt.colorbar(ax4.collections[0], ax=ax4)
        cb.set_label('時間 [s]')
        ax4.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(self.data_dir / 'energy_analysis.png', dpi=300, bbox_inches='tight')
        plt.show()
        print("エネルギー解析図を保存しました: energy_analysis.png")
    
    def plot_theory_comparison(self):
        """理論解との比較プロット"""
        if self.theory_data is None:
            print("理論解比較データがありません")
            return
        
        fig, ax = plt.subplots(2, 2, figsize=(15, 10))
        
        # パラメータ名と値を抽出
        parameters = self.theory_data['Parameter'].tolist()  # type: ignore
        theoretical = self.theory_data['Theoretical'].tolist()  # type: ignore
        actual = self.theory_data['Actual'].tolist()  # type: ignore
        errors = self.theory_data['Error_%'].tolist()  # type: ignore
        
        # 比較棒グラフ
        x = np.arange(len(parameters))
        width = 0.35
        
        ax[0,0].bar(x - width/2, theoretical, width, label='理論値', alpha=0.8)
        ax[0,0].bar(x + width/2, actual, width, label='計算値', alpha=0.8)
        ax[0,0].set_xlabel('パラメータ')
        ax[0,0].set_ylabel('値')
        ax[0,0].set_title('理論値vs計算値')
        ax[0,0].set_xticks(x)
        ax[0,0].set_xticklabels([p.replace('_', '\n') for p in parameters], rotation=45, ha='right')
        ax[0,0].legend()
        ax[0,0].grid(True, alpha=0.3)
        
        # 相対誤差
        ax[0,1].bar(range(len(parameters)), errors, alpha=0.8, color='red')
        ax[0,1].set_xlabel('パラメータ')
        ax[0,1].set_ylabel('相対誤差 [%]')
        ax[0,1].set_title('理論解との相対誤差')
        ax[0,1].set_xticks(range(len(parameters)))
        ax[0,1].set_xticklabels([p.replace('_', '\n') for p in parameters], rotation=45, ha='right')
        ax[0,1].grid(True, alpha=0.3)
        
        # 理論値と計算値の散布図
        ax[1,0].scatter(theoretical, actual, alpha=0.8, s=100)
        min_val = min(min(theoretical), min(actual))
        max_val = max(max(theoretical), max(actual))
        ax[1,0].plot([min_val, max_val], [min_val, max_val], 'r--', linewidth=2, label='理想線')
        ax[1,0].set_xlabel('理論値')
        ax[1,0].set_ylabel('計算値')
        ax[1,0].set_title('理論値vs計算値 散布図')
        ax[1,0].legend()
        ax[1,0].grid(True, alpha=0.3)
        
        # 誤差の円グラフ
        if len(errors) > 1:
            ax[1,1].pie(errors, labels=[p.replace('_', '\n') for p in parameters], autopct='%1.1f%%')
            ax[1,1].set_title('相対誤差の分布')
        else:
            ax[1,1].text(0.5, 0.5, '誤差データが\n不足しています', ha='center', va='center', transform=ax[1,1].transAxes)
        
        plt.tight_layout()
        plt.savefig(self.data_dir / 'theory_comparison.png', dpi=300, bbox_inches='tight')
        plt.show()
        print("理論解比較図を保存しました: theory_comparison.png")
    
    def plot_contact_analysis(self):
        """接触解析プロット"""
        if self.contact_data is None:
            print("接触解析データがありません")
            return
        
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 10))
        
        # 接触継続時間
        durations = self.contact_data['Duration_s'].tolist()  # type: ignore
        ax1.hist(durations, bins=20, alpha=0.7, edgecolor='black')
        ax1.set_xlabel('接触継続時間 [s]')
        ax1.set_ylabel('頻度')
        ax1.set_title('接触継続時間の分布')
        ax1.grid(True, alpha=0.3)
        
        # 衝突前後速度
        pre_vel = self.contact_data['Pre_Velocity_m/s'].tolist()  # type: ignore
        post_vel = self.contact_data['Post_Velocity_m/s'].tolist()  # type: ignore
        ax2.scatter(pre_vel, post_vel, alpha=0.8, s=100)
        max_vel = max(max(pre_vel), max(post_vel))
        ax2.plot([0, max_vel], [0, max_vel], 'r--', linewidth=2, label='v_post = v_pre')
        ax2.set_xlabel('衝突前速度 [m/s]')
        ax2.set_ylabel('衝突後速度 [m/s]')
        ax2.set_title('衝突前後速度の関係')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        # オーバーラップ解析
        overlap_ratios = self.contact_data['Overlap_Ratio_%'].tolist()  # type: ignore
        ax3.hist(overlap_ratios, bins=20, alpha=0.7, edgecolor='black')
        ax3.set_xlabel('オーバーラップ比 [%]')
        ax3.set_ylabel('頻度')
        ax3.set_title('最大オーバーラップ比の分布')
        ax3.grid(True, alpha=0.3)
        
        # 反発係数の検証
        calc_rest = self.contact_data['Calculated_Restitution'].tolist()  # type: ignore
        set_rest = self.contact_data['Set_Restitution'].tolist()  # type: ignore
        ax4.scatter(set_rest, calc_rest, alpha=0.8, s=100)
        min_rest = min(min(calc_rest), min(set_rest))
        max_rest = max(max(calc_rest), max(set_rest))
        ax4.plot([min_rest, max_rest], [min_rest, max_rest], 'r--', linewidth=2, label='理想線')
        ax4.set_xlabel('設定反発係数')
        ax4.set_ylabel('計算反発係数')
        ax4.set_title('反発係数の検証')
        ax4.legend()
        ax4.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(self.data_dir / 'contact_analysis.png', dpi=300, bbox_inches='tight')
        plt.show()
        print("接触解析図を保存しました: contact_analysis.png")
    
    def create_animation(self, interval=50, save_gif=True):
        """粒子の軌跡アニメーション作成"""
        times, x_pos, z_pos, z_vel = self.extract_trajectory_data()
        
        if times is None or x_pos is None or z_pos is None or z_vel is None:
            print("アニメーション作成に必要なデータが不足しています")
            return
        
        if len(times) < 2:
            print("アニメーション作成に必要なデータ点が不足しています")
            return
        
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # 軌跡の設定
        trail_length = min(50, len(z_pos) // 10)  # 軌跡の長さ
        particle_radius = 5.0  # 粒子半径（仮定）
        
        # プロット範囲の設定
        x_margin = particle_radius * 2
        z_margin = particle_radius * 2
        ax.set_xlim(min(x_pos) - x_margin, max(x_pos) + x_margin)
        ax.set_ylim(min(z_pos) - z_margin, max(z_pos) + z_margin)
        
        # 床の描画
        floor_y = 0
        ax.axhline(y=floor_y, color='brown', linewidth=3, label='床')
        
        # 初期化
        particle_circle = plt.Circle((x_pos[0], z_pos[0]), int(particle_radius), 
                                   color='blue', alpha=0.8)
        ax.add_patch(particle_circle)
        
        trail_line, = ax.plot([], [], 'r-', alpha=0.5, linewidth=2, label='軌跡')
        velocity_text = ax.text(0.02, 0.98, '', transform=ax.transAxes, 
                               verticalalignment='top', fontsize=12,
                               bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        ax.set_xlabel('x座標 [m]')
        ax.set_ylabel('z座標 [m]')
        ax.set_title('自由落下シミュレーション アニメーション')
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_aspect('equal')
        
        def animate(frame):
            if frame >= len(times):
                return particle_circle, trail_line, velocity_text
            
            # 粒子位置の更新
            particle_circle.center = (x_pos[frame], z_pos[frame])
            
            # 軌跡の更新
            start_idx = max(0, frame - trail_length)
            trail_line.set_data(x_pos[start_idx:frame+1], z_pos[start_idx:frame+1])
            
            # 速度情報の更新
            velocity_text.set_text(f'時刻: {times[frame]:.3f} s\n'
                                 f'高さ: {z_pos[frame]:.2f} m\n'
                                 f'速度: {z_vel[frame]:.2f} m/s')
            
            return particle_circle, trail_line, velocity_text
        
        # アニメーション作成
        anim = animation.FuncAnimation(fig, animate, frames=len(times), 
                                     interval=interval, blit=True, repeat=True)
        
        if save_gif:
            gif_path = self.data_dir / 'free_fall_animation.gif'
            print(f"アニメーションを保存中... ({gif_path})")
            anim.save(gif_path, writer='pillow', fps=20)
            print("アニメーション保存完了")
        
        plt.show()
        
    def generate_report(self):
        """総合レポートの生成"""
        print("\n========================================")
        print("         総合解析レポート")
        print("========================================")
        
        if self.summary_data is not None:
            print("シミュレーション概要:")
            for _, row in self.summary_data.iterrows():
                param = row['Parameter']
                value = row['Value']
                print(f"  {param}: {value}")
        
        if self.theory_data is not None:
            print("\n理論解との比較:")
            for _, row in self.theory_data.iterrows():
                param = row['Parameter']
                theoretical = row['Theoretical']
                actual = row['Actual']
                error = row['Error_%']
                print(f"  {param}:")
                print(f"    理論値: {theoretical:.6f}")
                print(f"    計算値: {actual:.6f}")
                print(f"    誤差: {error:.2f}%")
        
        if self.contact_data is not None:
            print("\n接触解析結果:")
            print(f"  接触回数: {len(self.contact_data)}")
            if len(self.contact_data) > 0:
                avg_duration = self.contact_data['Duration_s'].mean()  # type: ignore
                avg_overlap = self.contact_data['Overlap_Ratio_%'].mean()  # type: ignore
                print(f"  平均接触時間: {avg_duration:.6f} s")
                print(f"  平均オーバーラップ比: {avg_overlap:.2f}%")
        
        print("========================================")

def main():
    """メイン実行関数"""
    analyzer = FreeFallAnalyzer()
    
    print("\n解析メニュー:")
    print("1. 軌跡解析")
    print("2. エネルギー解析")
    print("3. 理論解比較")
    print("4. 接触解析")
    print("5. アニメーション作成")
    print("6. 全解析実行")
    print("7. レポート生成")
    
    choice = input("選択してください (1-7): ").strip()
    
    if choice == '1':
        analyzer.plot_trajectory()
    elif choice == '2':
        analyzer.plot_energy_analysis()
    elif choice == '3':
        analyzer.plot_theory_comparison()
    elif choice == '4':
        analyzer.plot_contact_analysis()
    elif choice == '5':
        analyzer.create_animation()
    elif choice == '6':
        print("全解析を実行します...")
        analyzer.plot_trajectory()
        analyzer.plot_energy_analysis()
        analyzer.plot_theory_comparison()
        analyzer.plot_contact_analysis()
        analyzer.create_animation()
        analyzer.generate_report()
    elif choice == '7':
        analyzer.generate_report()
    else:
        print("無効な選択です")

if __name__ == "__main__":
    main() 