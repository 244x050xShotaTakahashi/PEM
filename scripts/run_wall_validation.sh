#!/bin/bash

# 壁-粒子間検証実行スクリプト
# 使用方法: ./scripts/run_wall_validation.sh [validation_type] [input_file]
#
# validation_type:
#   freefall - 自由落下反発検証
#   slope    - 摩擦斜面検証
#   sweep    - パラメータ掃引検証
#   all      - 全ての検証を実行
#
# 例:
#   ./scripts/run_wall_validation.sh freefall
#   ./scripts/run_wall_validation.sh slope input/input_slope_validation.dat
#   ./scripts/run_wall_validation.sh all

set -e  # エラー時に停止

# スクリプトのディレクトリを取得
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

# プロジェクトディレクトリに移動
cd "$PROJECT_DIR"

# 関数定義
print_header() {
    echo "=================================="
    echo "$1"
    echo "=================================="
}

print_info() {
    echo "[INFO] $1"
}

print_error() {
    echo "[ERROR] $1" >&2
}

check_executable() {
    if [[ ! -x "./pem_simulator" ]]; then
        print_error "実行ファイル pem_simulator が見つかりません"
        print_info "make コマンドでコンパイルしてください"
        exit 1
    fi
}

create_data_dir() {
    if [[ ! -d "data" ]]; then
        mkdir -p data
        print_info "dataディレクトリを作成しました"
    fi
}

run_freefall_validation() {
    print_header "自由落下反発検証を実行中..."
    
    local input_file="${1:-input/input_wall_validation.dat}"
    
    if [[ ! -f "$input_file" ]]; then
        print_error "入力ファイルが見つかりません: $input_file"
        return 1
    fi
    
    print_info "入力ファイル: $input_file"
    
    # シミュレーション実行
    ./pem_simulator "$input_file"
    
    print_info "自由落下反発検証が完了しました"
    
    # 結果ファイルの確認
    if [[ -f "data/wall_validation_freefall_report.txt" ]]; then
        print_info "検証レポートが生成されました: data/wall_validation_freefall_report.txt"
    fi
}

run_slope_validation() {
    print_header "摩擦斜面検証を実行中..."
    
    local input_file="${1:-input/input_slope_validation.dat}"
    
    if [[ ! -f "$input_file" ]]; then
        print_error "入力ファイルが見つかりません: $input_file"
        return 1
    fi
    
    print_info "入力ファイル: $input_file"
    
    # シミュレーション実行
    ./pem_simulator "$input_file"
    
    print_info "摩擦斜面検証が完了しました"
    
    # 結果ファイルの確認
    if [[ -f "data/wall_validation_slope_report.txt" ]]; then
        print_info "検証レポートが生成されました: data/wall_validation_slope_report.txt"
    fi
}

run_parameter_sweep() {
    print_header "パラメータ掃引検証を実行中..."
    
    local input_file="${1:-input/input_wall_validation_sweep.dat}"
    
    if [[ ! -f "$input_file" ]]; then
        print_error "入力ファイルが見つかりません: $input_file"
        return 1
    fi
    
    print_info "入力ファイル: $input_file"
    
    # シミュレーション実行
    ./pem_simulator "$input_file"
    
    print_info "パラメータ掃引検証が完了しました"
    
    # 結果ファイルの確認
    if [[ -f "data/parameter_sweep_results.csv" ]]; then
        print_info "掃引結果が生成されました: data/parameter_sweep_results.csv"
    fi
}

run_analysis() {
    print_header "結果分析を実行中..."
    
    # Python解析スクリプトの実行
    if command -v python3 &> /dev/null; then
        if [[ -f "scripts/wall_validation_analyzer.py" ]]; then
            python3 scripts/wall_validation_analyzer.py
        else
            print_error "分析スクリプトが見つかりません: scripts/wall_validation_analyzer.py"
        fi
    else
        print_info "Python3が見つかりません。分析スクリプトをスキップします"
    fi
}

# メイン処理
main() {
    local validation_type="${1:-all}"
    local input_file="$2"
    
    print_header "PEM 壁-粒子間検証実行スクリプト"
    print_info "検証タイプ: $validation_type"
    
    # 前提条件チェック
    check_executable
    create_data_dir
    
    case "$validation_type" in
        "freefall")
            run_freefall_validation "$input_file"
            ;;
        "slope")
            run_slope_validation "$input_file"
            ;;
        "sweep")
            run_parameter_sweep "$input_file"
            ;;
        "all")
            run_freefall_validation "$input_file"
            echo ""
            run_slope_validation "$input_file"
            echo ""
            run_parameter_sweep "$input_file"
            ;;
        *)
            print_error "不明な検証タイプ: $validation_type"
            echo "使用可能なタイプ: freefall, slope, sweep, all"
            exit 1
            ;;
    esac
    
    echo ""
    run_analysis
    
    print_header "全ての検証が完了しました"
    print_info "結果は data/ ディレクトリに保存されています"
}

# スクリプト実行
main "$@"
