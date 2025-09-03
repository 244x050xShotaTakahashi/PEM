# 改良版自由落下シミュレーション用Makefile

# コンパイラとフラグの設定
FC = ifort
FFLAGS = -O3 -r8 -traceback -g -warn all
TARGET = pem_simulator
SOURCE = src/pem_simulator.f90

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
	rm -f $(TARGET) src/*.mod src/*.o *.mod *.o
	rm -f data/state_step_*.csv
	
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