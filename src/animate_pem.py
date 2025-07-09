import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
from matplotlib.patches import Circle
from pathlib import Path
from typing import Union

DATA_FILE_DEFAULT = Path(__file__).resolve().parent.parent / "data" / "graph11.d"

def read_simulation_data(filename: Union[str, Path] = DATA_FILE_DEFAULT):
    """graph11.d を読み込み、各タイムステップのデータを辞書のリストとして返す。

    今回はヘッダーが 2 行に分かれているケース（例: 3 つの値 + 次行に rmax のみ）にも対応する。
    読み込みはトークンストリームとして処理し、必要な個数だけ値を順次取り出す方式に変更した。
    """

    def _read_floats(token_buffer, fh, n):
        """token_buffer (list[str]) から n 個の float を取り出して返す。
        不足している場合はファイルから行を読み、トークンを補充する。"""
        vals = []
        while len(vals) < n:
            if not token_buffer:
                line = fh.readline()
                if not line:
                    raise EOFError("予期せぬ EOF")
                token_buffer.extend(line.split())
            # token_buffer が空でなければ pop
            vals.append(float(token_buffer.pop(0)))
        return vals

    frames_data = []
    try:
        with open(filename, "r") as fh:
            tokens: list[str] = []  # 行を跨いで残ったトークンを一時保持
            while True:
                try:
                    # ヘッダー: num_particles, time, container_width
                    num_particles_f, time_val, container_width = _read_floats(tokens, fh, 3)
                    num_particles = int(num_particles_f)
                    # rmax は次の 1 値
                    rmax_val, = _read_floats(tokens, fh, 1) # 使わないが読み飛ばす
                except EOFError:
                    break  # 正常終了

                if num_particles < 0:
                    print("警告: num_particles が負の値です。スキップします。")
                    continue

                particles_data_to_read = num_particles * 3
                
                # 粒子データ (x,z,r) × num_particles
                particle_values = []
                if num_particles > 0:
                    particle_values = _read_floats(tokens, fh, particles_data_to_read)

                # 速度データも同数だけ存在する。読み飛ばす。
                velocity_values = []
                if num_particles > 0:
                    velocity_values = _read_floats(tokens, fh, particles_data_to_read)

                # 回転角度データを読み込む
                rotation_angles = []
                if num_particles > 0:
                    rotation_angles = _read_floats(tokens, fh, num_particles)

                particles = [
                    {
                        "x": particle_values[i * 3 + 0],
                        "z": particle_values[i * 3 + 1],
                        "r": particle_values[i * 3 + 2],
                        "rotation_angle": rotation_angles[i] if i < len(rotation_angles) else 0.0,
                    }
                    for i in range(num_particles)
                ]

                frames_data.append(
                    {
                        "time": time_val,
                        "num_particles": num_particles,
                        "container_width": container_width,
                        "particles": particles,
                    }
                )
    except FileNotFoundError:
        print(f"エラー: ファイル '{filename}' が見つかりません。")
        return None
    except Exception as e:
        print(f"ファイル読み込み中にエラーが発生しました: {e}")
        return None
        
    return frames_data

def animate(frames_data, output_filename="pem_animation.gif"):
    if not frames_data:
        print("アニメーションするデータがありません。")
        return

    # フォント設定を無効化して英語のみ使用
    import matplotlib
    matplotlib.rcParams['font.family'] = 'DejaVu Sans'
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # 描画範囲とアスペクト比は最初に一度だけ設定
    container_width = frames_data[0]['container_width']
    
    # 検証モード用の表示設定
    is_validation_mode = False
    is_slope_validation = False
    if frames_data[0]['num_particles'] == 1 or frames_data[0]['num_particles'] == 2:
        is_validation_mode = True
        
        # 斜面検証かどうかを判定（粒子のx方向移動量で判定）
        if len(frames_data) > 10 and frames_data[0]['num_particles'] == 1:
            initial_x = frames_data[0]['particles'][0]['x']
            later_x = frames_data[min(10, len(frames_data)-1)]['particles'][0]['x']
            x_movement = abs(later_x - initial_x)
            if x_movement > 0.5:  # x方向に大きく移動している場合は斜面検証
                is_slope_validation = True
    
    max_z_overall = 0.0
    if frames_data[0]['num_particles'] > 0 and frames_data[0]['particles']:
         max_z_overall = max(p['z'] + p['r'] for p in frames_data[0]['particles']) 
    else: 
         max_z_overall = container_width * 0.5 if container_width > 0 else 1.0
    
    for frame_data in frames_data:
        if frame_data['num_particles'] > 0 and frame_data['particles']:
            current_max_z = max(p['z'] + p['r'] for p in frame_data['particles'])
            if current_max_z > max_z_overall:
                max_z_overall = current_max_z
                
    if max_z_overall <= 0: # フォールバック
        max_z_overall = container_width * 0.5 if container_width > 0 else 1.0
    
    # 検証モード用の表示範囲調整
    if is_validation_mode:
        max_z_overall = max(max_z_overall, 12.0)  # 自由落下検証用に高さを確保

    # update_frameの外で固定のテキストオブジェクトを一度だけ作成
    time_text_obj = ax.text(0.05, 0.95, '', transform=ax.transAxes, ha="left", va="top", fontsize=10)


    def update_frame(frame_idx):
        ax.cla() # 現在のアックスの内容をすべてクリア

        # selected_framesから実際のフレームデータを取得
        data = selected_frames[frame_idx]
        
        # クリア後、軸の範囲やラベル、タイトルを再設定
        current_container_width_for_xlim = data['container_width']
        ax.set_xlim(-0.1 * current_container_width_for_xlim, 1.1 * current_container_width_for_xlim)
        ax.set_ylim(-0.1 * max_z_overall , 1.2 * max_z_overall)
        ax.set_aspect('equal', adjustable='box') 
        ax.set_xlabel("X coordinate")
        ax.set_ylabel("Z coordinate")
        fig.suptitle(f"PEM Validation (Frame {frame_idx+1}/{len(selected_frames)})", fontsize=12)
        current_container_width = data['container_width']

        # 壁の描画
        ax.plot([0, 0], [0, max_z_overall * 1.15], 'k-', lw=2) # 左壁
        ax.plot([0, current_container_width], [0, 0], 'k-', lw=2) # 床
        ax.plot([current_container_width, current_container_width], [0, max_z_overall * 1.15], 'k-', lw=2) # 右壁
        
        # 斜面壁の描画 (摩擦斜面検証用のみ)
        if is_validation_mode and is_slope_validation and data['num_particles'] == 1:
            # 30度の斜面を描画
            slope_angle = 30.0 * np.pi / 180.0  # 30度をラジアンに変換
            slope_x_end = min(current_container_width, max_z_overall * 1.15 / np.tan(slope_angle))
            slope_z_end = slope_x_end * np.tan(slope_angle)
            ax.plot([0, slope_x_end], [0, slope_z_end], 'r-', lw=2, alpha=0.7, label='Slope Wall')
            if frame_idx == 0:  # 最初のフレームでのみ凡例を追加
                ax.legend(loc='upper right')

        for p_data in data['particles']:
            circle = Circle((p_data['x'], p_data['z']), p_data['r'], facecolor='white', edgecolor='black', linewidth=1.5, alpha=0.9)
            ax.add_patch(circle)
            
            # 回転を示す線を描画
            x_center = p_data['x']
            z_center = p_data['z']
            radius = p_data['r']
            rotation_angle = p_data['rotation_angle']
            
            # 回転角度の方向に線を描画（0度は右方向、反時計回りが正）
            x_end = x_center + radius * np.cos(rotation_angle)
            z_end = z_center + radius * np.sin(rotation_angle)
            
            ax.plot([x_center, x_end], [z_center, z_end], 'k-', linewidth=1, alpha=0.8)
        
        time_text_obj.set_text(f"Time: {data['time']:.6f} s")
        ax.add_artist(time_text_obj) 
        
        return [] 

    # フレーム数を制限してパフォーマンスを向上
    max_frames = min(len(frames_data), 200)  # 最大200フレームに制限
    frame_step = max(1, len(frames_data) // max_frames)
    
    selected_frames = frames_data[::frame_step]
    
    ani = animation.FuncAnimation(fig, update_frame, frames=len(selected_frames), blit=False, interval=150)

    try:
        print(f"Animation saving to '{output_filename}' (frames: {len(selected_frames)})...")
        writer = animation.PillowWriter(fps=10)  # fpsを下げてファイルサイズを削減
        ani.save(output_filename, writer=writer, dpi=80)  # dpiを下げて軽量化
        
        print("Animation saved successfully.")
    except KeyboardInterrupt:
        print("Animation creation was interrupted.")
        return
    except Exception as e:
        print(f"Error during animation creation: {e}")
        print("Please check if Pillow is installed correctly.")
        print("Example: pip install Pillow")


if __name__ == "__main__":
    # PEM/src から 1 つ上の PEM へ移動し data/graph11.d を指す
    data_file = DATA_FILE_DEFAULT  # Path オブジェクト
    
    all_frames_data = read_simulation_data(data_file)
    
    if all_frames_data:
        animate(all_frames_data, "pem_animation.gif")
    else:
        print(f"{data_file} からデータを読み込めませんでした。")
