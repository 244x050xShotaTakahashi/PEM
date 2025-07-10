# 改良版自由落下シミュレーション用Makefile

# コンパイラとフラグの設定
FC = gfortran
FFLAGS = -O3 -fdefault-real-8 -fbacktrace -g -Wall
TARGET = pem_simulator_improved
SOURCE = pem_simulator_improved.f90

# デフォルトターゲット
all: $(TARGET)

# メインターゲット
$(TARGET): $(SOURCE)
	$(FC) $(FFLAGS) -o $(TARGET) $(SOURCE)

# 実行ターゲット
run: $(TARGET)
	./$(TARGET)

# クリーンアップ
clean:
	rm -f $(TARGET) *.mod *.o
	rm -rf data/

# データディレクトリの作成
prepare:
	mkdir -p data data/detailed data/analysis

# 完全なビルドと実行
build-run: clean prepare $(TARGET) run

# ヘルプ
help:
	@echo "利用可能なターゲット:"
	@echo "  all       - プログラムをコンパイル"
	@echo "  run       - プログラムを実行"
	@echo "  clean     - 生成ファイルを削除"
	@echo "  prepare   - データディレクトリを作成"
	@echo "  build-run - クリーン、コンパイル、実行を一括実行"
	@echo "  help      - このヘルプを表示"

.PHONY: all run clean prepare build-run help 