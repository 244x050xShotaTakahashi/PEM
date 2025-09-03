#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
壁-粒子間検証結果の自動分析スクリプト

このスクリプトは、PEMシミュレーターによる壁-粒子間検証の結果を分析し、
詳細なレポートと可視化グラフを生成します。

使用方法:
    python scripts/wall_validation_analyzer.py [data_directory]

機能:
    - 自由落下反発検証の結果分析
    - 摩擦斜面検証の結果分析
    - パラメータ掃引結果の分析
    - エラー統計の計算
    - 結果の可視化（グラフ生成）
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
import argparse
import json
from datetime import datetime


class WallValidationAnalyzer:
    """壁-粒子間検証結果の分析クラス"""
    
    def __init__(self, data_dir="data"):
        self.data_dir = Path(data_dir)
        self.results = {}
        self.plots_dir = self.data_dir / "plots"
        self.plots_dir.mkdir(exist_ok=True)
        
        # 日本語フォント設定
        plt.rcParams['font.family'] = 'DejaVu Sans'
        plt.rcParams['axes.unicode_minus'] = False
        
    def analyze_freefall_validation(self):
        """自由落下反発検証の結果を分析"""
        print("自由落下反発検証の分析を開始...")
        
        # パラメータ掃引結果の読み込み
        sweep_file = self.data_dir / "parameter_sweep_results.csv"
        if sweep_file.exists():
            df = pd.read_csv(sweep_file)
            
            # エラー統計の計算
            mean_error = df['Error_Percent'].mean()
            max_error = df['Error_Percent'].max()
            min_error = df['Error_Percent'].min()
            std_error = df['Error_Percent'].std()
            
            self.results['freefall'] = {
                'mean_error': mean_error,
                'max_error': max_error,
                'min_error': min_error,
                'std_error': std_error,
                'data_points': len(df)
            }
            
            # 可視化
            self._plot_parameter_sweep_results(df)
            
            print(f"  平均誤差: {mean_error:.4f}%")
            print(f"  最大誤差: {max_error:.4f}%")
            print(f"  最小誤差: {min_error:.4f}%")
            print(f"  標準偏差: {std_error:.4f}%")
            
        else:
            print("  パラメータ掃引結果ファイルが見つかりません")
    
    def analyze_slope_validation(self):
        """摩擦斜面検証の結果を分析"""
        print("摩擦斜面検証の分析を開始...")
        
        # レポートファイルから情報を抽出
        report_file = self.data_dir / "wall_validation_slope_report.txt"
        if report_file.exists():
            with open(report_file, 'r', encoding='utf-8') as f:
                content = f.read()
                print("  斜面検証レポートを発見しました")
                self.results['slope'] = {'report_available': True}
        else:
            print("  斜面検証レポートが見つかりません")
            self.results['slope'] = {'report_available': False}
    
    def _plot_parameter_sweep_results(self, df):
        """パラメータ掃引結果の可視化"""
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
        
        # 反発係数vs誤差のプロット
        ax1.plot(df['Restitution_Coefficient'], df['Error_Percent'], 'bo-', linewidth=2, markersize=8)
        ax1.set_xlabel('反発係数', fontsize=12)
        ax1.set_ylabel('相対誤差 (%)', fontsize=12)
        ax1.set_title('反発係数vs相対誤差', fontsize=14)
        ax1.grid(True, alpha=0.3)
        ax1.set_ylim(bottom=0)
        
        # 理論値vs計算値の比較
        ax2.plot(df['Theoretical_Height'], df['Actual_Height'], 'ro', markersize=8, label='計算値')
        min_val = min(df['Theoretical_Height'].min(), df['Actual_Height'].min())
        max_val = max(df['Theoretical_Height'].max(), df['Actual_Height'].max())
        ax2.plot([min_val, max_val], [min_val, max_val], 'k--', linewidth=2, label='理論値=計算値')
        ax2.set_xlabel('理論反発高さ', fontsize=12)
        ax2.set_ylabel('計算反発高さ', fontsize=12)
        ax2.set_title('理論値vs計算値', fontsize=14)
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        ax2.axis('equal')
        
        plt.tight_layout()
        plt.savefig(self.plots_dir / "freefall_validation_analysis.png", dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"  分析結果を保存: {self.plots_dir / 'freefall_validation_analysis.png'}")
    
    def generate_summary_report(self):
        """総合的な分析レポートを生成"""
        report_file = self.data_dir / "wall_validation_summary.txt"
        
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write("=" * 60 + "\n")
            f.write("PEM壁-粒子間検証 総合分析レポート\n")
            f.write("=" * 60 + "\n")
            f.write(f"生成日時: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            
            # 自由落下反発検証の結果
            if 'freefall' in self.results:
                f.write("--- 自由落下反発検証 ---\n")
                freefall = self.results['freefall']
                f.write(f"データ点数: {freefall['data_points']}\n")
                f.write(f"平均誤差: {freefall['mean_error']:.4f}%\n")
                f.write(f"最大誤差: {freefall['max_error']:.4f}%\n")
                f.write(f"最小誤差: {freefall['min_error']:.4f}%\n")
                f.write(f"標準偏差: {freefall['std_error']:.4f}%\n")
                
                # 精度評価
                if freefall['mean_error'] < 1.0:
                    f.write("評価: 高精度 (平均誤差 < 1%)\n")
                elif freefall['mean_error'] < 5.0:
                    f.write("評価: 良好 (平均誤差 < 5%)\n")
                else:
                    f.write("評価: 要改善 (平均誤差 >= 5%)\n")
                f.write("\n")
            
            # 斜面検証の結果
            if 'slope' in self.results:
                f.write("--- 摩擦斜面検証 ---\n")
                if self.results['slope']['report_available']:
                    f.write("検証レポート: 利用可能\n")
                else:
                    f.write("検証レポート: 未実行\n")
                f.write("\n")
            
            f.write("--- 推奨事項 ---\n")
            f.write("1. 時間刻みを適切に設定してください（安定性条件を満たす）\n")
            f.write("2. 材料物性値が実験値と一致することを確認してください\n")
            f.write("3. 数値誤差を最小化するため、十分な数値精度を使用してください\n")
            f.write("4. 複数の検証ケースで結果の一貫性を確認してください\n")
            f.write("\n")
            f.write("=" * 60 + "\n")
        
        print(f"総合分析レポートを生成: {report_file}")
    
    def run_analysis(self):
        """全ての分析を実行"""
        print("壁-粒子間検証の分析を開始します...")
        print(f"データディレクトリ: {self.data_dir}")
        
        self.analyze_freefall_validation()
        self.analyze_slope_validation()
        self.generate_summary_report()
        
        # 結果をJSONでも保存
        results_file = self.data_dir / "validation_analysis_results.json"
        with open(results_file, 'w', encoding='utf-8') as f:
            json.dump(self.results, f, indent=2, ensure_ascii=False)
        
        print(f"分析結果を保存: {results_file}")
        print("分析完了！")


def main():
    """メイン関数"""
    parser = argparse.ArgumentParser(description="壁-粒子間検証結果の分析")
    parser.add_argument('data_dir', nargs='?', default='data',
                       help='データディレクトリのパス (デフォルト: data)')
    
    args = parser.parse_args()
    
    # 分析器を作成して実行
    analyzer = WallValidationAnalyzer(args.data_dir)
    analyzer.run_analysis()


if __name__ == "__main__":
    main()
