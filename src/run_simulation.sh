#!/bin/bash

# 改良版自由落下シミュレーション実行スクリプト
# dem_calc.py準拠の詳細解析機能付き

echo "========================================"
echo "  改良版自由落下シミュレーション"
echo "  (dem_calc.py準拠の詳細解析機能付き)"
echo "========================================"

# スクリプトのディレクトリに移動
cd "$(dirname "$0")"

# 必要なディレクトリの作成
echo "出力ディレクトリを準備中..."
mkdir -p data
mkdir -p data/detailed
mkdir -p data/analysis

# 既存の出力ファイルをクリーンアップ
echo "既存の出力ファイルをクリーンアップ中..."
rm -f data/*.csv
rm -f data/detailed/*.csv
rm -f data/analysis/*.csv
rm -f data/*.png
rm -f data/*.gif

# Fortranプログラムのコンパイル
echo "Fortranプログラムをコンパイル中..."
if [ -f "Makefile_improved" ]; then
    make -f Makefile_improved clean
    make -f Makefile_improved all
else
    # Makefileがない場合は直接コンパイル
    gfortran -O3 -fdefault-real-8 -fbacktrace -g -Wall -o pem_simulator_improved pem_simulator_improved.f90
fi

# コンパイルが成功したかチェック
if [ ! -f "pem_simulator_improved" ]; then
    echo "エラー: コンパイルに失敗しました"
    exit 1
fi

echo "コンパイル完了"

# 入力ファイルの確認
if [ ! -f "../input/input_free_fall.dat" ]; then
    echo "警告: ../input/input_free_fall.dat が見つかりません"
    echo "デフォルト設定で実行します"
fi

# シミュレーション実行
echo ""
echo "シミュレーションを開始します..."
echo "========================================"

./pem_simulator_improved

# 実行結果の確認
if [ $? -eq 0 ]; then
    echo ""
    echo "========================================"
    echo "シミュレーション完了"
    echo "========================================"
    
    # 生成されたファイルの確認
    echo "生成されたファイル:"
    find data -name "*.csv" -type f | head -10
    
    # 可視化の実行（Pythonがインストールされている場合）
    if command -v python3 &> /dev/null; then
        echo ""
        echo "結果の可視化を実行しますか？ (y/n)"
        read -r response
        if [[ "$response" =~ ^[Yy]$ ]]; then
            echo "可視化スクリプトを実行中..."
            python3 visualize_results.py
        fi
    else
        echo ""
        echo "Python3が見つかりません。手動で可視化スクリプトを実行してください:"
        echo "python3 visualize_results.py"
    fi
    
    echo ""
    echo "全ての処理が完了しました"
    
else
    echo ""
    echo "エラー: シミュレーションの実行に失敗しました"
    exit 1
fi 