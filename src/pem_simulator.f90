! メインプログラムおよびモジュール: 2次元粒子要素法 (球体モデル)
! 提供されたPDFのi-pem.fに基づく。

! モジュール: シミュレーション定数 (配列サイズ、数学定数)
module simulation_constants_mod
    implicit none
    integer, parameter :: ni_max = 1000  ! ni: 最大粒子数
    integer, parameter :: nj_max = 13    ! nj: 粒子ごとの最大接触点数 (粒子間10 + 壁3)
    integer, parameter :: nc_max = 20000 ! nc: グリッド内の最大セル数
    real(8), parameter :: PI_VAL = 3.141592653589793d0 ! pi: 円周率
    real(8), parameter :: GRAVITY_ACCEL = 9.80665d0    ! g: 重力加速度
end module simulation_constants_mod

! モジュール: シミュレーション制御パラメータと材料物性値
module simulation_parameters_mod
    use simulation_constants_mod, only: ni_max
    implicit none

    ! シミュレーション制御パラメータ
    real(8) :: time_step                      ! dt: 時間刻み
    real(8) :: friction_coeff_particle        ! fri: 粒子間摩擦係数
    real(8) :: friction_coeff_wall            ! frw: 壁-粒子間摩擦係数
    real(8) :: young_modulus_particle         ! e: 粒子のヤング率
    real(8) :: young_modulus_wall             ! ew: 壁のヤング率
    real(8) :: poisson_ratio_particle         ! po: 粒子のポアソン比
    real(8) :: poisson_ratio_wall             ! pow: 壁のポアソン比
    real(8) :: shear_to_normal_stiffness_ratio ! so: せん断弾性係数と法線方向弾性係数の比
    real(8) :: particle_density               ! de: 粒子の密度
    
    ! 粒子生成パラメータ
    real(8) :: particle_radius_large          ! r1: 大きな粒子の半径
    real(8) :: particle_radius_small          ! r2: 小さな粒子の半径
    real(8) :: container_width                ! w: 容器の幅
    integer :: particle_gen_layers            ! ipz: 初期粒子生成層数
    integer :: random_seed                    ! 乱数シード
    
    ! セル法アルゴリズム制御パラメータ
    logical :: disable_cell_algorithm         ! セル法アルゴリズムを無効化するフラグ
    real(8) :: cell_size_override            ! セルサイズの手動設定値 (0.0=自動計算)
    
    ! 検証モード設定
    logical :: validation_mode                ! 検証モードフラグ
    real(8) :: validation_particle1_x        ! 検証モード: 粒子1のx座標
    real(8) :: validation_particle1_z        ! 検証モード: 粒子1のz座標
    real(8) :: validation_particle2_x        ! 検証モード: 粒子2のx座標
    real(8) :: validation_particle2_z        ! 検証モード: 粒子2のz座標
    real(8) :: validation_particle1_vx       ! 検証モード: 粒子1の初期x速度
    real(8) :: validation_particle2_vx       ! 検証モード: 粒子2の初期x速度
    real(8) :: validation_particle_radius    ! 検証モード: 粒子半径
    
    ! 出力制御パラメータ
    integer :: output_interval_normal         ! 通常モード出力間隔
    integer :: output_interval_validation     ! 検証モード出力間隔
    integer :: max_calculation_steps          ! 最大計算ステップ数

    ! 粒子-壁間検証モード設定
    logical :: wall_validation_mode           ! 粒子-壁間検証モードフラグ
    integer :: wall_validation_type           ! 検証タイプ (1: 自由落下反発, 2: 摩擦斜面)
    real(8) :: wall_validation_drop_height    ! 落下開始高さ
    real(8) :: wall_validation_particle_radius ! 検証用粒子半径
    real(8) :: wall_validation_particle_x     ! 検証用粒子のx座標
    real(8) :: wall_validation_restitution_coeff ! 理論計算用反発係数
    real(8) :: wall_validation_friction_coeff    ! 理論計算用摩擦係数
    real(8) :: wall_validation_slope_angle      ! 斜面角度（ラジアン）
    real(8) :: wall_validation_slope_friction   ! 斜面摩擦係数
    
    ! パラメータ掃引設定
    logical :: parameter_sweep_mode           ! パラメータ掃引モードフラグ
    integer :: parameter_sweep_count          ! 掃引回数
    real(8) :: restitution_coeff_min          ! 反発係数の最小値
    real(8) :: restitution_coeff_max          ! 反発係数の最大値

    save
end module simulation_parameters_mod

! モジュール: 粒子固有データ (物理特性、運動学、力)
module particle_data_mod
    use simulation_constants_mod, only: ni_max, nj_max
    implicit none

    ! 物理特性
    real(8), dimension(ni_max) :: radius         ! rr(ni): 粒子半径
    real(8), dimension(ni_max) :: mass           ! wei(ni): 粒子質量
    real(8), dimension(ni_max) :: moment_inertia ! pmi(ni): 粒子の慣性モーメント

    ! 位置と向き
    real(8), dimension(ni_max) :: x_coord        ! x0(ni): 粒子中心のx座標
    real(8), dimension(ni_max) :: z_coord        ! z0(ni): 粒子中心のz座標 (原文ではy、コードではz)
    real(8), dimension(ni_max) :: rotation_angle ! qq(ni): 粒子の回転変位 (角度)

    ! 速度 (並進および回転)
    real(8), dimension(ni_max) :: x_vel          ! u0(ni): 粒子のx方向速度
    real(8), dimension(ni_max) :: z_vel          ! v0(ni): 粒子のz方向速度
    real(8), dimension(ni_max) :: rotation_vel   ! f0(ni): 粒子の回転速度

    ! 合力とモーメント
    real(8), dimension(ni_max) :: x_force_sum    ! xf(ni): 粒子に働くx方向の合力
    real(8), dimension(ni_max) :: z_force_sum    ! zf(ni): 粒子に働くz方向の合力
    real(8), dimension(ni_max) :: moment_sum     ! mf(ni): 粒子に働くモーメント

    ! 接触力の成分と接触相手のインデックス
    real(8), dimension(ni_max, nj_max) :: normal_force_contact  ! en(ni,nj): 法線方向接触力
    real(8), dimension(ni_max, nj_max) :: shear_force_contact   ! es(ni,nj): せん断方向接触力
    integer, dimension(ni_max, nj_max) :: contact_partner_idx ! je(ni,nj): 接触点番号配列 (接触している粒子/壁のインデックスを格納)

    ! 現時間ステップにおける増分変位 (common/dpm/ より)
    real(8), dimension(ni_max) :: x_disp_incr    ! u(ni): x方向変位増分
    real(8), dimension(ni_max) :: z_disp_incr    ! v(ni): z方向変位増分
    real(8), dimension(ni_max) :: rot_disp_incr  ! f(ni): 回転変位増分
    
    save
end module particle_data_mod

! モジュール: セル格子システムデータ
module cell_system_mod
    use simulation_constants_mod, only: ni_max, nc_max
    implicit none

    integer :: num_particles          ! n: 粒子数
    integer :: cells_x_dir            ! idx: x方向のセル数
    integer :: cells_z_dir            ! idz: z方向のセル数 (使用状況から推測)

    real(8) :: cell_size              ! c: セルの幅/サイズ

    ! 各セルに属する粒子を連結リストで保持する（セル先頭 index → next → ...）
    integer, dimension(nc_max) :: cell_head        ! そのセルで最初に登録された粒子インデックス (空=0)
    integer, dimension(ni_max) :: particle_cell_next ! 次の粒子インデックス (0 ならリスト終端)

    ! 後方互換用に「最後に登録された粒子」を保持（デバッグ用途）
    integer, dimension(nc_max) :: cell_particle_map

    integer, dimension(ni_max) :: particle_cell_idx ! 粒子iが格納されているセル番号

    save
end module cell_system_mod

! メインプログラム
program two_dimensional_pem
    use simulation_constants_mod
    use simulation_parameters_mod
    use particle_data_mod
    use cell_system_mod
    implicit none

    integer :: it_step, static_judge_flag          ! static_judge_flag: 静止判定フラグ
    integer :: i ! ループカウンタ用にiを宣言
    real(8) :: current_time                        ! 現在時刻
    real(8) :: rmax_particle_radius                ! fpositから返される最大粒子半径
    
    ! 計算時間計測用変数
    integer :: start_time, end_time, clock_rate
    real(8) :: elapsed_time
    
    ! --- ▼ 1D 衝突検証用追加変数 ▼ ---
    logical :: collision_started = .false.  ! 衝突が開始したか
    logical :: collision_finished = .false. ! 衝突が終了したか
    real(8) :: initial_v1, initial_v2        ! 衝突前速度
    real(8) :: final_v1, final_v2            ! 衝突後速度
    real(8) :: dist, sumr                    ! 接触判定用
    real(8) :: m1, m2, e_coeff, v1_theo, v2_theo, err1, err2 ! 評価用
    real(8) :: ke_initial, ke_final, ke_theo ! エネルギー保存チェック用
    
    ! --- ▼ 粒子-壁間検証用追加変数 ▼ ---
    logical :: wall_collision_started = .false.  ! 壁衝突が開始したか
    logical :: wall_collision_finished = .false. ! 壁衝突が終了したか
    logical :: particle_falling = .true.         ! 粒子が落下中か
    logical :: particle_rising = .false.         ! 粒子が上昇中か
    real(8) :: initial_drop_height              ! 初期落下高さ
    real(8) :: pre_collision_velocity           ! 衝突直前の速度
    real(8) :: post_collision_velocity          ! 衝突直後の速度
    real(8) :: max_rebound_height               ! 最大反発高さ
    real(8) :: theoretical_rebound_height       ! 理論的反発高さ
    real(8) :: height_error                     ! 高さの誤差
    
    ! 摩擦斜面検証用変数
    real(8) :: theoretical_velocity, theoretical_position
    real(8) :: actual_velocity, actual_position
    real(8) :: velocity_error, position_error

    ! 計算時間計測開始
    call system_clock(start_time, clock_rate)
    
    ! inputファイルからパラメータを読み込み
    call read_input_file
    
    ! パラメータ掃引モードかどうかの判定と実行
    if (parameter_sweep_mode .and. wall_validation_mode .and. wall_validation_type == 1) then
        call parameter_sweep_validation
        stop
    end if
    
    ! 初期位置と初期条件の設定
    call fposit_sub(rmax_particle_radius)
    call inmat_sub
    call init_sub

    ! -----------------------------------------------
    ! 検証モードの場合：初期水平速度と摩擦係数を設定
    ! -----------------------------------------------
    if (validation_mode) then
        ! 粒子1に指定された速度、粒子2は静止
        x_vel(1) = validation_particle1_vx
        x_vel(2) = validation_particle2_vx
        z_vel(1) = 0.0d0
        z_vel(2) = 0.0d0

        ! 理論値計算用の初期速度保存
        initial_v1 = x_vel(1)
        initial_v2 = x_vel(2)

        ! 摩擦を無効化
        friction_coeff_particle = 0.0d0
        friction_coeff_wall     = 0.0d0
    end if

    ! -----------------------------------------------
    ! 粒子-壁間検証モードの場合：初期速度と摩擦係数を設定
    ! -----------------------------------------------
    if (wall_validation_mode) then
        if (wall_validation_type == 1) then
            ! 自由落下反発検証: 初速ゼロ、重力あり
            x_vel(1) = 0.0d0
            z_vel(1) = 0.0d0
            rotation_vel(1) = 0.0d0
            ! 摩擦を無効化してエネルギー保存に集中
            friction_coeff_particle = 0.0d0
            friction_coeff_wall = 0.0d0
        else if (wall_validation_type == 2) then
            ! 摩擦斜面検証: 初速ゼロ、重力あり、摩擦あり
            x_vel(1) = 0.0d0
            z_vel(1) = 0.0d0
            rotation_vel(1) = 0.0d0
            ! 指定された摩擦係数を使用
            friction_coeff_particle = wall_validation_friction_coeff
            friction_coeff_wall = wall_validation_slope_friction
        end if
    end if

    ! -----------------------------------------------
    ! 粒子-壁間検証モードの初期設定
    ! -----------------------------------------------
    if (wall_validation_mode .and. wall_validation_type == 1) then
        initial_drop_height = wall_validation_drop_height
        theoretical_rebound_height = wall_validation_restitution_coeff**2 * initial_drop_height
        write(*,*) '自由落下反発検証パラメータ:'
        write(*,*) '  初期落下高さ: ', initial_drop_height
        write(*,*) '  反発係数: ', wall_validation_restitution_coeff
        write(*,*) '  理論反発高さ: ', theoretical_rebound_height
    end if

    current_time = 0.0d0

    ! 各ステップの繰り返し計算
    do it_step = 1, max_calculation_steps
        current_time = current_time + time_step

        ! 粒子をセルに格納
        call ncel_sub

        ! 全粒子の合力をクリア
        do i = 1, num_particles
            x_force_sum(i) = 0.0d0
            z_force_sum(i) = 0.0d0
            moment_sum(i) = 0.0d0
        end do
        
        do i = 1, num_particles
            ! 粒子と壁との接触力計算
            call wcont_sub(i)
            ! 粒子間の接触力計算
            call pcont_sub(i, rmax_particle_radius)
        end do

        ! 増分変位の重ね合わせ (運動方程式の積分)
        call nposit_sub(static_judge_flag)

        ! === ▼ 1D 衝突検証ロジック ▼ ===
        if (validation_mode .and. .not. collision_finished) then
            dist = dabs(x_coord(1) - x_coord(2))
            sumr = radius(1) + radius(2)

            ! 衝突開始判定
            if (.not. collision_started .and. dist < sumr) then
                collision_started = .true.
                write(*,*) '衝突開始: 時刻 = ', current_time
                write(*,*) '衝突前速度: v1 = ', x_vel(1), ' v2 = ', x_vel(2)
            end if

            ! 衝突終了判定: 一度接触後、再び離れた
            if (collision_started .and. dist > sumr + 1.0d-6) then
                collision_finished = .true.

                final_v1 = x_vel(1)
                final_v2 = x_vel(2)

                ! 正確な質量を使った理論計算 (完全弾性衝突)
                m1 = mass(1)
                m2 = mass(2)
                e_coeff = 1.0d0  ! 完全弾性衝突

                ! 運動量保存とエネルギー保存から導出される公式
                v1_theo = ((m1 - e_coeff*m2)*initial_v1 + (1.0d0+e_coeff)*m2*initial_v2)/(m1 + m2)
                v2_theo = ((1.0d0+e_coeff)*m1*initial_v1 + (m2 - e_coeff*m1)*initial_v2)/(m1 + m2)

                ! 同質量の場合の簡単チェック
                if (dabs(m1 - m2) < 1.0d-12) then
                    ! 同質量・完全弾性衝突では速度が交換される
                    v1_theo = initial_v2
                    v2_theo = initial_v1
                end if

                ! 相対誤差計算
                if (dabs(v1_theo) > 1.0d-12) then
                    err1 = dabs((final_v1 - v1_theo)/v1_theo) * 100.0d0
                else
                    err1 = dabs(final_v1 - v1_theo) * 100.0d0
                end if
                if (dabs(v2_theo) > 1.0d-12) then
                    err2 = dabs((final_v2 - v2_theo)/v2_theo) * 100.0d0
                else
                    err2 = dabs(final_v2 - v2_theo) * 100.0d0
                end if

                write(*,*) '================================='
                write(*,*) '一次元衝突検証結果'
                write(*,*) '================================='
                write(*,'(A,ES12.5,A,ES12.5)') '質量: m1 = ', m1, ' m2 = ', m2
                write(*,'(A,F10.6,A,F10.6)') '初期速度: v1_init = ', initial_v1, ' v2_init = ', initial_v2
                write(*,'(A,F10.6,A,F10.6)') '理論最終速度: v1_theo = ', v1_theo, ' v2_theo = ', v2_theo
                write(*,'(A,F10.6,A,F10.6)') '計算最終速度: v1_calc = ', final_v1, ' v2_calc = ', final_v2
                write(*,'(A,F8.4,A,A,F8.4,A)') '相対誤差: v1誤差 = ', err1, '%', '  v2誤差 = ', err2, '%'
                
                ! エネルギー保存チェック
                ke_initial = 0.5d0 * (m1 * initial_v1**2 + m2 * initial_v2**2)
                ke_final = 0.5d0 * (m1 * final_v1**2 + m2 * final_v2**2)
                ke_theo = 0.5d0 * (m1 * v1_theo**2 + m2 * v2_theo**2)
                write(*,'(A,ES12.5)') '初期運動エネルギー: ', ke_initial
                write(*,'(A,ES12.5)') '最終運動エネルギー: ', ke_final
                write(*,'(A,ES12.5)') '理論運動エネルギー: ', ke_theo
                write(*,'(A,F8.4,A)') 'エネルギー保存誤差: ', dabs((ke_final-ke_initial)/ke_initial)*100.0d0, '%'
                write(*,*) '================================='
                ! 衝突検証完了後も計算を続行 (goto 200を削除)
                goto 200
            end if
        end if

        ! === ▼ 粒子-壁間検証ロジック ▼ ===
        if (wall_validation_mode .and. wall_validation_type == 1 .and. .not. wall_collision_finished) then
            ! 自由落下反発検証
            if (particle_falling .and. z_coord(1) <= radius(1) + 1.0d-6) then
                ! 床に衝突
                wall_collision_started = .true.
                pre_collision_velocity = abs(z_vel(1))
                write(*,*) '床衝突検出: 時刻 = ', current_time
                write(*,*) '衝突直前速度: ', pre_collision_velocity
            else if (wall_collision_started .and. particle_falling .and. z_coord(1) > radius(1) + 1.0d-3) then
                ! 床から離れた（反発開始）
                particle_falling = .false.
                particle_rising = .true.
                write(*,*) '床から離脱: 時刻 = ', current_time, ', 位置 = ', z_coord(1)
            else if (particle_rising .and. z_vel(1) < 0.0d0) then
                ! 上昇から落下に転じた（最高点到達）
                particle_rising = .false.
                particle_falling = .true.
                max_rebound_height = z_coord(1)
                wall_collision_finished = .true.
                
                ! 理論値との比較
                height_error = abs(max_rebound_height - theoretical_rebound_height) / theoretical_rebound_height * 100.0d0
                
                write(*,*) '最高点検出: 時刻 = ', current_time, ', 位置 = ', z_coord(1), ', 速度 = ', z_vel(1)
                write(*,*) '================================='
                write(*,*) '自由落下反発検証結果'
                write(*,*) '================================='
                write(*,*) '初期落下高さ: ', initial_drop_height
                write(*,*) '反発係数: ', wall_validation_restitution_coeff
                write(*,*) '理論反発高さ: ', theoretical_rebound_height
                write(*,*) '計算反発高さ: ', max_rebound_height
                write(*,*) '相対誤差: ', height_error, '%'
                                    write(*,*) '================================='
                    
                    exit  ! シミュレーションループを抜ける
            end if
        else if (wall_validation_mode .and. wall_validation_type == 2 .and. .not. wall_collision_finished) then
            ! 摩擦斜面検証
            if (it_step == 1) then
                ! 初期状態の記録
                initial_drop_height = z_coord(1)
                write(*,*) '摩擦斜面検証パラメータ:'
                write(*,*) '  初期位置: x=', x_coord(1), ', z=', z_coord(1)
                write(*,*) '  斜面角度: ', wall_validation_slope_angle, ' ラジアン'
                write(*,*) '  摩擦係数: ', wall_validation_slope_friction
                write(*,*) '  理論加速度: ', GRAVITY_ACCEL * (sin(wall_validation_slope_angle) - wall_validation_slope_friction * cos(wall_validation_slope_angle))
            end if
            
            if (it_step == 5000) then
                ! 一定時間後の理論値との比較
                theoretical_velocity = GRAVITY_ACCEL * (sin(wall_validation_slope_angle) - wall_validation_slope_friction * cos(wall_validation_slope_angle)) * current_time
                theoretical_position = wall_validation_particle_x + 0.5d0 * GRAVITY_ACCEL * (sin(wall_validation_slope_angle) - wall_validation_slope_friction * cos(wall_validation_slope_angle)) * current_time**2
                
                actual_velocity = sqrt(x_vel(1)**2 + z_vel(1)**2)
                actual_position = x_coord(1)
                
                velocity_error = abs(actual_velocity - theoretical_velocity) / theoretical_velocity * 100.0d0
                position_error = abs(actual_position - theoretical_position) / theoretical_position * 100.0d0
                
                write(*,*) '================================='
                write(*,*) '摩擦斜面検証結果'
                write(*,*) '================================='
                write(*,*) '時刻: ', current_time
                write(*,*) '理論速度: ', theoretical_velocity
                write(*,*) '計算速度: ', actual_velocity
                write(*,*) '速度誤差: ', velocity_error, '%'
                write(*,*) '理論位置: ', theoretical_position
                write(*,*) '計算位置: ', actual_position
                write(*,*) '位置誤差: ', position_error, '%'
                write(*,*) '================================='
                
                wall_collision_finished = .true.
            end if
        end if

        ! 静止状態の判定
        if (static_judge_flag == 1) then
            write(*,*) '静止状態に到達しました。時刻: ', current_time
            goto 200 ! シミュレーションループを抜ける
        end if

        ! 計算状況の出力
        if (mod(it_step, 100) == 0) then
            if (wall_validation_mode .and. wall_validation_type == 1) then
                write(*, '(A,F10.6,A,F12.6,A,F12.6,A,L1)') 'Time= ', current_time, &
                                                     ' Z0(1)= ', z_coord(1), &
                                                     ' V0(1)= ', z_vel(1), &
                                                     ' Falling= ', particle_falling
            else
                write(*, '(A,F10.6,A,F12.6,A,F12.6)') 'Time= ', current_time, &
                                                     ' Z0(N)= ', z_coord(num_particles), &
                                                     ' V0(N)= ', z_vel(num_particles)
            end if
        end if

        ! グラフィック用データの出力
        if (validation_mode .or. wall_validation_mode) then
            ! 検証モードでは指定間隔で出力
            if (it_step == 1 .or. mod(it_step, output_interval_validation) == 0) then
                call gfout_sub(it_step, current_time, rmax_particle_radius)
            end if
        else
            ! 通常モードでは指定間隔で出力
            if (it_step == 1 .or. mod(it_step, output_interval_normal) == 0) then
                call gfout_sub(it_step, current_time, rmax_particle_radius)
            end if
        end if
    end do

200 continue ! シミュレーションループ脱出用のラベル

    ! バックアップデータの出力
    call bfout_sub

    close(10) ! data/graph11.d (グラフィックデータファイル1)
    close(11) ! data/graph21.d (グラフィックデータファイル2)
    close(13) ! data/backl.d (バックアップファイル、bfout_subで開かれていれば)

    ! 計算時間計測終了
    call system_clock(end_time)
    elapsed_time = real(end_time - start_time) / real(clock_rate)
    
    write(*,*) '================================='
    write(*,*) 'シミュレーション実行結果'
    write(*,*) '================================='
    write(*,*) '粒子数: ', num_particles
    write(*,*) '計算ステップ数: ', it_step
    write(*,*) '実行時間: ', elapsed_time, ' 秒'
    write(*,*) '1ステップあたりの平均時間: ', elapsed_time / real(it_step), ' 秒'
    
    if (disable_cell_algorithm .or. cell_size_override > 0.0d0) then
        write(*,*) 'セル法アルゴリズム: 無効化'
        write(*,*) 'セルサイズ: ', cell_size
    else
        write(*,*) 'セル法アルゴリズム: 有効'
        write(*,*) 'セルサイズ: ', cell_size
        write(*,*) 'セル数 (X方向): ', cells_x_dir
        write(*,*) 'セル数 (Z方向): ', cells_z_dir
        write(*,*) '総セル数: ', cells_x_dir * cells_z_dir
    end if
    
    write(*,*) '================================='

    stop
contains

    !> inputファイルからパラメータを読み込むサブルーチン
    subroutine read_input_file
        implicit none
        character(len=256) :: line, keyword
        character(len=256) :: input_filename
        integer :: ios, unit_num
        real(8) :: value
        
        ! inputファイル名の決定（コマンドライン引数または固定名）
        if (command_argument_count() > 0) then
            call get_command_argument(1, input_filename)
            ! 相対パスの場合、inputフォルダを追加
            if (input_filename(1:1) /= '/' .and. input_filename(1:5) /= 'input') then
                input_filename = 'input/' // trim(input_filename)
            end if
        else
            input_filename = "input/input.dat"
        end if
        
        write(*,*) 'inputファイルを読み込み中: ', trim(input_filename)
        
        unit_num = 20
        open(unit=unit_num, file=input_filename, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(*,*) 'エラー: inputファイルを開けません: ', trim(input_filename)
            stop
        end if
        
        ! デフォルト値の設定
        time_step = 5.0d-7
        friction_coeff_particle = 0.25d0
        friction_coeff_wall = 0.17d0
        young_modulus_particle = 4.9d9
        young_modulus_wall = 3.9d9
        poisson_ratio_particle = 0.23d0
        poisson_ratio_wall = 0.25d0
        particle_density = 2.48d3
        particle_radius_large = 1.0d-2
        particle_radius_small = 5.0d-3
        container_width = 5.0d-1
        particle_gen_layers = 30
        random_seed = 584287
        disable_cell_algorithm = .false.
        cell_size_override = 0.0d0
        validation_mode = .false.
        validation_particle1_x = 3.0d0
        validation_particle1_z = 3.0d0
        validation_particle2_x = 6.0d0
        validation_particle2_z = 3.0d0
        validation_particle1_vx = 20.0d0
        validation_particle2_vx = 0.0d0
        validation_particle_radius = 1.0d0
        output_interval_normal = 50000
        output_interval_validation = 10
        max_calculation_steps = 2000000
        
        ! 粒子-壁間検証モードのデフォルト値
        wall_validation_mode = .false.
        wall_validation_type = 1
        wall_validation_drop_height = 10.0d0
        wall_validation_particle_radius = 1.0d0
        wall_validation_particle_x = 5.0d0
        wall_validation_restitution_coeff = 0.8d0
        wall_validation_friction_coeff = 0.3d0
        wall_validation_slope_angle = 0.52359877559d0  ! 30度をラジアンに変換
        wall_validation_slope_friction = 0.3d0
        
        ! パラメータ掃引モードのデフォルト値
        parameter_sweep_mode = .false.
        parameter_sweep_count = 5
        restitution_coeff_min = 0.2d0
        restitution_coeff_max = 0.9d0
        
        do
            read(unit_num, '(A)', iostat=ios) line
            if (ios /= 0) exit
            
            ! コメント行と空行をスキップ
            if (line(1:1) == '#' .or. line(1:1) == '!' .or. len_trim(line) == 0) cycle
            
            ! キーワードと値を分離
            read(line, *, iostat=ios) keyword, value
            if (ios /= 0) cycle
            
            select case (trim(keyword))
                case ('TIME_STEP')
                    time_step = value
                case ('FRICTION_COEFF_PARTICLE')
                    friction_coeff_particle = value
                case ('FRICTION_COEFF_WALL')
                    friction_coeff_wall = value
                case ('YOUNG_MODULUS_PARTICLE')
                    young_modulus_particle = value
                case ('YOUNG_MODULUS_WALL')
                    young_modulus_wall = value
                case ('POISSON_RATIO_PARTICLE')
                    poisson_ratio_particle = value
                case ('POISSON_RATIO_WALL')
                    poisson_ratio_wall = value
                case ('PARTICLE_DENSITY')
                    particle_density = value
                case ('PARTICLE_RADIUS_LARGE')
                    particle_radius_large = value
                case ('PARTICLE_RADIUS_SMALL')
                    particle_radius_small = value
                case ('CONTAINER_WIDTH')
                    container_width = value
                case ('PARTICLE_GEN_LAYERS')
                    particle_gen_layers = int(value)
                case ('RANDOM_SEED')
                    random_seed = int(value)
                case ('DISABLE_CELL_ALGORITHM')
                    disable_cell_algorithm = (int(value) == 1)
                case ('CELL_SIZE_OVERRIDE')
                    cell_size_override = value
                case ('VALIDATION_MODE')
                    validation_mode = (int(value) == 1)
                case ('VALIDATION_PARTICLE1_X')
                    validation_particle1_x = value
                case ('VALIDATION_PARTICLE1_Z')
                    validation_particle1_z = value
                case ('VALIDATION_PARTICLE2_X')
                    validation_particle2_x = value
                case ('VALIDATION_PARTICLE2_Z')
                    validation_particle2_z = value
                case ('VALIDATION_PARTICLE1_VX')
                    validation_particle1_vx = value
                case ('VALIDATION_PARTICLE2_VX')
                    validation_particle2_vx = value
                case ('VALIDATION_PARTICLE_RADIUS')
                    validation_particle_radius = value
                case ('OUTPUT_INTERVAL_NORMAL')
                    output_interval_normal = int(value)
                case ('OUTPUT_INTERVAL_VALIDATION')
                    output_interval_validation = int(value)
                case ('MAX_CALCULATION_STEPS')
                    max_calculation_steps = int(value)
                case ('WALL_VALIDATION_MODE')
                    wall_validation_mode = (int(value) == 1)
                case ('WALL_VALIDATION_TYPE')
                    wall_validation_type = int(value)
                case ('WALL_VALIDATION_DROP_HEIGHT')
                    wall_validation_drop_height = value
                case ('WALL_VALIDATION_PARTICLE_RADIUS')
                    wall_validation_particle_radius = value
                case ('WALL_VALIDATION_PARTICLE_X')
                    wall_validation_particle_x = value
                case ('WALL_VALIDATION_RESTITUTION_COEFF')
                    wall_validation_restitution_coeff = value
                case ('WALL_VALIDATION_FRICTION_COEFF')
                    wall_validation_friction_coeff = value
                case ('WALL_VALIDATION_SLOPE_ANGLE')
                    wall_validation_slope_angle = value
                case ('WALL_VALIDATION_SLOPE_FRICTION')
                    wall_validation_slope_friction = value
                case ('PARAMETER_SWEEP_MODE')
                    parameter_sweep_mode = (int(value) == 1)
                case ('PARAMETER_SWEEP_COUNT')
                    parameter_sweep_count = int(value)
                case ('RESTITUTION_COEFF_MIN')
                    restitution_coeff_min = value
                case ('RESTITUTION_COEFF_MAX')
                    restitution_coeff_max = value
                case default
                    write(*,*) '警告: 不明なキーワード: ', trim(keyword)
            end select
        end do
        
        close(unit_num)
        
        write(*,*) 'inputファイルの読み込み完了'
        if (validation_mode) then
            write(*,*) '粒子間衝突検証モードで実行します'
        else if (wall_validation_mode) then
            write(*,*) '粒子-壁間検証モードで実行します'
            if (wall_validation_type == 1) then
                write(*,*) '検証タイプ: 自由落下反発'
            else if (wall_validation_type == 2) then
                write(*,*) '検証タイプ: 摩擦斜面'
            end if
        else
            write(*,*) '通常モードで実行します'
        end if
    end subroutine read_input_file

    !> 初期粒子配置と構成を設定するサブルーチン
    subroutine fposit_sub(rmax_out)
        use simulation_constants_mod, only: ni_max, PI_VAL
        use simulation_parameters_mod
        use particle_data_mod, only: radius, x_coord, z_coord, rotation_angle
        use cell_system_mod ! モジュールからセルシステム関連変数を取得
        implicit none

        real(8), intent(out) :: rmax_out ! 出力: 最大粒子半径

        integer :: i_layer, j_particle_in_layer, ipx_calc, current_particle_count
        real(8) :: r1_val, r2_val, rn_val, dx_offset, random_uniform_val
        real(8) :: rmin_val ! 宣言をここに移動
        integer :: particles_this_row ! 宣言をここに移動
        
        ! -------------------------------------------------
        ! 検証モード: 2 粒子のみを手動配置し通常生成をスキップ
        ! -------------------------------------------------
        if (validation_mode) then
            num_particles = 2

            ! 半径
            radius(1) = validation_particle_radius
            radius(2) = validation_particle_radius

            ! 位置
            x_coord(1) = validation_particle1_x
            z_coord(1) = validation_particle1_z
            x_coord(2) = validation_particle2_x
            z_coord(2) = validation_particle2_z

            rotation_angle(1) = 0.0d0
            rotation_angle(2) = 0.0d0

            rmax_out = validation_particle_radius
            rmin_val = validation_particle_radius

            ! 検証モード用のコンテナ設定
            container_width = max(validation_particle1_x, validation_particle2_x) + 2.0d0 * validation_particle_radius
            cell_size = validation_particle_radius * 2.0d0  ! セルサイズを大きくして両方の粒子を同じセルに配置

            cells_x_dir = idint(container_width / cell_size) + 1
            cells_z_dir = 5

            if (cells_x_dir * cells_z_dir > nc_max) then
                write(*,*) 'セル数がnc_maxを超えています (検証モード)'
                stop 'fposit_sub: セル配列が小さすぎます'
            end if

            return
        end if

        ! -------------------------------------------------
        ! 粒子-壁間検証モード: 1粒子を手動配置
        ! -------------------------------------------------
        if (wall_validation_mode) then
            num_particles = 1

            ! 半径
            radius(1) = wall_validation_particle_radius

            ! 位置
            x_coord(1) = wall_validation_particle_x
            if (wall_validation_type == 1) then
                ! 自由落下反発検証: 指定された高さから落下
                z_coord(1) = wall_validation_drop_height
            else if (wall_validation_type == 2) then
                ! 摩擦斜面検証: 斜面上に配置
                z_coord(1) = wall_validation_drop_height
            end if

            rotation_angle(1) = 0.0d0

            rmax_out = wall_validation_particle_radius
            rmin_val = wall_validation_particle_radius

            ! 粒子-壁間検証モード用のコンテナ設定
            container_width = max(wall_validation_particle_x * 2.0d0, 10.0d0)
            cell_size = wall_validation_particle_radius * 2.0d0

            cells_x_dir = idint(container_width / cell_size) + 1
            cells_z_dir = idint((wall_validation_drop_height * 2.0d0) / cell_size) + 10

            if (cells_x_dir * cells_z_dir > nc_max) then
                write(*,*) 'セル数がnc_maxを超えています (粒子-壁間検証モード)'
                stop 'fposit_sub: セル配列が小さすぎます'
            end if

            return
        end if

        ! -------------------------------------------------
        ! 通常モード: ランダム配置による粒子生成
        ! -------------------------------------------------

        ! 粒子半径 (inputファイルから取得)
        r1_val = particle_radius_large
        r2_val = particle_radius_small

        rmax_out = r1_val             ! 最大半径をr1_valとする
        rmin_val = r2_val             ! 最小半径をr2_valとする

        rn_val = rmax_out + 1.0d-5    ! パッキングのための有効半径
        ipx_calc = idint(container_width / (2.0d0 * rn_val)) ! 1行あたりの粒子数 (概算)

        current_particle_count = 0
        do i_layer = 1, particle_gen_layers
            if (mod(i_layer, 2) == 0) then  ! 偶数層
                dx_offset = 2.0d0 * rn_val
                particles_this_row = ipx_calc - 1
            else                            ! 奇数層
                dx_offset = rn_val
                particles_this_row = ipx_calc
            end if

            do j_particle_in_layer = 1, particles_this_row
                call custom_random(random_seed, random_uniform_val)
                if (random_uniform_val < 2.0d-1) cycle ! 一部の位置をスキップ

                current_particle_count = current_particle_count + 1
                if (current_particle_count > ni_max) then
                    write(*,*) '粒子数がni_maxを超えました: ', ni_max
                    stop 'fposit_sub: 粒子が多すぎます'
                end if
                num_particles = current_particle_count ! グローバルな粒子数を更新

                x_coord(num_particles) = 2.0d0 * rn_val * (j_particle_in_layer - 1) + dx_offset
                z_coord(num_particles) = 2.0d0 * rn_val * (i_layer - 1) + rn_val
                rotation_angle(num_particles) = 0.0d0 ! 回転角を初期化

                call custom_random(random_seed, random_uniform_val)
                if (random_uniform_val < 0.5d0) then
                    radius(num_particles) = r1_val
                else
                    radius(num_particles) = r2_val
                end if
            end do
        end do
        write(*,*) '生成された粒子数: ', num_particles

        ! セルサイズ計算 (原文PDF p.35 eq 3.25: C < sqrt(2)*rmin)
        if (disable_cell_algorithm .or. cell_size_override > 0.0d0) then
            ! セル法アルゴリズムを無効化する場合
            if (cell_size_override > 0.0d0) then
                cell_size = cell_size_override
            else
                ! 計算領域全体を1つのセルとして設定
                cell_size = max(container_width, 2.0d0 * rmax_out * particle_gen_layers) * 2.0d0
            end if
            write(*,*) 'セル法アルゴリズムを無効化: cell_size = ', cell_size
        else
            ! 通常のセルサイズ計算
            if (rmin_val > 0.0d0) then
                 cell_size = rmin_val * 1.30d0 ! または入力から。原文ではrmin*1.35d0はコメントアウト
            else
                 cell_size = rmax_out * 1.30d0 ! rminが適切に定義されない場合のフォールバック
            end if
        end if
        
        if (cell_size <= 0.0d0) then
            write(*,*) "エラー: fposit_subでcell_sizeが正ではありません。"
            stop
        endif

        cells_x_dir = idint(container_width / cell_size) + 1
        ! cells_z_dirは粒子が到達しうる最大高さをカバーする必要がある
        ! 元のコード: idz=idint(z0(n)/c)+10。生成時の最上部粒子のz座標を使用。
        if (num_particles > 0 .and. cell_size > 0.0d0) then
            cells_z_dir = idint(z_coord(num_particles) / cell_size) + 10 
        else if (particle_gen_layers > 0 .and. cell_size > 0.0d0 .and. rn_val > 0.0d0) then ! 粒子がない場合でも推定
             if (particle_gen_layers > 0 .and. rn_val > 0 .and. cell_size > 0) then
                cells_z_dir = idint( (2.0d0 * rn_val * (real(particle_gen_layers) -1.0d0) + rn_val) / cell_size) + 10
             else
                cells_z_dir = 20 ! Fallback if values are still problematic
             end if
        else
            cells_z_dir = 20 ! デフォルト値 (粒子も層もない、またはセルサイズが0の場合)
        end if


        if (cells_x_dir * cells_z_dir > nc_max) then
            write(*,*) 'ncl (cell_particle_map)がオーバーフローしました!! 要求セル数: ', cells_x_dir * cells_z_dir
            stop 'fposit_sub: セル配列が小さすぎます'
        end if

    end subroutine fposit_sub

    !> 材料物性値を初期化し、定数を計算するサブルーチン
    subroutine inmat_sub
        use simulation_constants_mod, only: PI_VAL
        use simulation_parameters_mod, only: time_step, validation_mode
        use particle_data_mod, only: radius, mass, moment_inertia 
        use cell_system_mod, only: num_particles 
        implicit none
        integer :: i
        
        ! time_step, particle_densityなどの値はモジュールsimulation_parameters_modで設定されていると仮定
        ! 粒子のポアソン比に基づいてせん断弾性係数と法線方向弾性係数の比(so)を計算
        shear_to_normal_stiffness_ratio = 1.0d0 / (2.0d0 * (1.0d0 + poisson_ratio_particle))

        do i = 1, num_particles
            ! 質量: 3D球体 V = 4/3 pi r^3。2Dディスク (面積 pi r^2)の場合、deが面密度ならば。
            ! 元のコードは2Dシミュレーションの文脈でも3D球体の体積で質量を計算しているように見える。
            mass(i) = (4.0d0 / 3.0d0) * PI_VAL * radius(i)**3 * particle_density
            
            ! 慣性モーメント: 3D球体 I = 2/5 m r^2。
            ! 元のコード: pmi(i)=8.d0/15.d0*de*pi*(rr(i)**5)
            ! これは (2/5) * (4/3 pi r^3 de) * r^2 = (2/5) * mass * r^2。球体として正しい。
            moment_inertia(i) = (8.0d0 / 15.0d0) * particle_density * PI_VAL * (radius(i)**5)
        end do
    end subroutine inmat_sub

    !> 接触力関連の配列を初期化するサブルーチン
    subroutine init_sub
        use simulation_constants_mod, only: nj_max
        use particle_data_mod
        use cell_system_mod, only: num_particles
        implicit none
        integer :: i, j

        ! 実際の粒子数まで繰り返す
        if (num_particles > 0) then
            do i = 1, num_particles 
                do j = 1, nj_max
                    normal_force_contact(i, j) = 0.0d0
                    shear_force_contact(i, j) = 0.0d0
                    contact_partner_idx(i, j) = 0
                end do
            end do
            ! 増分変位を最初にゼロで初期化
            x_disp_incr(1:num_particles) = 0.0d0
            z_disp_incr(1:num_particles) = 0.0d0
            rot_disp_incr(1:num_particles) = 0.0d0
        end if
    end subroutine init_sub

    !> 近傍探索のために粒子をセルに割り当てるサブルーチン
    subroutine ncel_sub
        use simulation_constants_mod, only: nc_max
        use particle_data_mod, only: x_coord, z_coord
        use cell_system_mod
        implicit none
        integer :: i, cell_block_idx
        integer :: ix_cell, iz_cell ! 宣言をここに移動

        ! 連結リストをクリア
        if (nc_max > 0) cell_head(1:nc_max) = 0
        if (nc_max > 0) cell_particle_map(1:nc_max) = 0  ! デバッグ用途
        if (num_particles > 0) particle_cell_next(1:num_particles) = 0

        do i = 1, num_particles
            particle_cell_idx(i) = 0 ! 初期化
            if (cell_size <= 0.0d0) then
                write(*,*) "エラー: ncel_subでcell_sizeが0または負です。"
                stop
            endif
            if (cells_x_dir <= 0) then
                 write(*,*) "エラー: ncel_subでcells_x_dirが0または負です。"
                 stop
            endif

            ! 粒子iを含むセルの1次元インデックスを計算
            ix_cell = idint(x_coord(i) / cell_size) ! 座標が負にならないように注意
            iz_cell = idint(z_coord(i) / cell_size)

            ! インデックスが有効範囲 [0, cells_x_dir-1] および [0, cells_z_dir-1] 内にあることを保証
            ix_cell = max(0, min(ix_cell, cells_x_dir - 1))
            iz_cell = max(0, min(iz_cell, cells_z_dir - 1))
            
            cell_block_idx = iz_cell * cells_x_dir + ix_cell + 1 ! Fortranの1ベースインデックス

            if (cell_block_idx > 0 .and. cell_block_idx <= nc_max) then
                 ! 連結リストの先頭に追加
                 particle_cell_next(i) = cell_head(cell_block_idx)
                 cell_head(cell_block_idx) = i

                 ! 旧 single-map も更新（最後に登録された粒子）
                 cell_particle_map(cell_block_idx) = i

                 particle_cell_idx(i) = cell_block_idx
            else
                write(*,*) 'エラー: ncel_subで粒子', i, 'のcell_block_idxが範囲外です。'
                write(*,*) 'x0, z0, c: ', x_coord(i), z_coord(i), cell_size
                write(*,*) 'ix_cell, iz_cell, idx, computed_block_idx: ', ix_cell, iz_cell, cells_x_dir, cell_block_idx
                stop 'ncel_sub: 粒子が無効なセルインデックスにマッピングされました'
            end if
            
            ! デバッグ出力 (検証モードの最初の呼び出しのみ)
            ! if (validation_mode .and. i <= 2) then
            !     write(*,'(A,I2,A,I5,A,F10.6,A,F10.6)') 'Particle ',i,' -> cell ',cell_block_idx,', pos=(',x_coord(i),',',z_coord(i),')'
            ! end if
        end do
    end subroutine ncel_sub

    !> 粒子iと壁との接触力を計算するサブルーチン
    subroutine wcont_sub(particle_idx)
        use simulation_parameters_mod, only: time_step, validation_mode, container_width
        use particle_data_mod
        use cell_system_mod, only: num_particles
        implicit none
        integer, intent(in) :: particle_idx ! 対象の粒子インデックス

        real(8) :: xi, zi, ri_particle
        real(8) :: wall_angle_sin, wall_angle_cos, overlap_gap
        integer :: wall_contact_slot_idx, wall_partner_id
        real(8) :: slope_distance, slope_normal_x, slope_normal_z ! 斜面壁用変数

        xi = x_coord(particle_idx)
        zi = z_coord(particle_idx)
        ri_particle = radius(particle_idx)

        ! 左壁 (contact_partner_idx = num_particles + 1)
        wall_contact_slot_idx = 11 ! 元のコードでの左壁用の固定スロット
        wall_partner_id = num_particles + 1
        if (xi < ri_particle) then  ! 左壁と接触
            wall_angle_sin = 0.0d0  ! 法線ベクトル成分 sin(alpha_ij) (粒子中心から壁中心へ向かうベクトル)
            wall_angle_cos = -1.0d0 ! 法線ベクトル成分 cos(alpha_ij)
            overlap_gap = ri_particle - xi ! 元のコードでは dabs(xi)、ここでは重なり量を正とする
            contact_partner_idx(particle_idx, wall_contact_slot_idx) = wall_partner_id
            if (validation_mode) then
                call actf_sub(particle_idx, wall_partner_id, wall_contact_slot_idx, wall_angle_sin, wall_angle_cos, overlap_gap, 1.0d0, 1.0d0)
            else if (wall_validation_mode) then
                call actf_sub(particle_idx, wall_partner_id, wall_contact_slot_idx, wall_angle_sin, wall_angle_cos, overlap_gap, wall_validation_restitution_coeff, wall_validation_restitution_coeff)
            else
                call actf_sub(particle_idx, wall_partner_id, wall_contact_slot_idx, wall_angle_sin, wall_angle_cos, overlap_gap, 0.5d0, 0.5d0)
            end if
        else                        ! 接触なし
            normal_force_contact(particle_idx, wall_contact_slot_idx) = 0.0d0
            shear_force_contact(particle_idx, wall_contact_slot_idx) = 0.0d0
            contact_partner_idx(particle_idx, wall_contact_slot_idx) = 0
        end if

        ! 下壁 (contact_partner_idx = num_particles + 2)
        wall_contact_slot_idx = 12 ! 元のコードでの下壁用の固定スロット
        wall_partner_id = num_particles + 2
        if (zi < ri_particle) then  ! 下壁と接触
            wall_angle_sin = -1.0d0
            wall_angle_cos = 0.0d0
            overlap_gap = ri_particle - zi 
            contact_partner_idx(particle_idx, wall_contact_slot_idx) = wall_partner_id
            if (validation_mode) then
                call actf_sub(particle_idx, wall_partner_id, wall_contact_slot_idx, wall_angle_sin, wall_angle_cos, overlap_gap, 1.0d0, 1.0d0)
            else if (wall_validation_mode) then
                ! デバッグ出力
                if (overlap_gap > 1.0d-6) then
                    write(*,*) 'DEBUG: 下壁接触 - 粒子:', particle_idx, ', 重なり:', overlap_gap, ', 反発係数:', wall_validation_restitution_coeff
                end if
                call actf_sub(particle_idx, wall_partner_id, wall_contact_slot_idx, wall_angle_sin, wall_angle_cos, overlap_gap, wall_validation_restitution_coeff, wall_validation_restitution_coeff)
            else
                call actf_sub(particle_idx, wall_partner_id, wall_contact_slot_idx, wall_angle_sin, wall_angle_cos, overlap_gap, 0.5d0, 0.5d0)
            end if
        else                        ! 接触なし
            normal_force_contact(particle_idx, wall_contact_slot_idx) = 0.0d0
            shear_force_contact(particle_idx, wall_contact_slot_idx) = 0.0d0
            contact_partner_idx(particle_idx, wall_contact_slot_idx) = 0
        end if

        ! 右壁 (contact_partner_idx = num_particles + 3)
        wall_contact_slot_idx = 13 ! 元のコードでの右壁用の固定スロット
        wall_partner_id = num_particles + 3
        if (xi + ri_particle > container_width) then ! 右壁と接触
            wall_angle_sin = 0.0d0
            wall_angle_cos = 1.0d0
            overlap_gap = (xi + ri_particle) - container_width 
            contact_partner_idx(particle_idx, wall_contact_slot_idx) = wall_partner_id
            if (validation_mode) then
                call actf_sub(particle_idx, wall_partner_id, wall_contact_slot_idx, wall_angle_sin, wall_angle_cos, overlap_gap, 1.0d0, 1.0d0)
            else if (wall_validation_mode) then
                call actf_sub(particle_idx, wall_partner_id, wall_contact_slot_idx, wall_angle_sin, wall_angle_cos, overlap_gap, wall_validation_restitution_coeff, wall_validation_restitution_coeff)
            else
                call actf_sub(particle_idx, wall_partner_id, wall_contact_slot_idx, wall_angle_sin, wall_angle_cos, overlap_gap, 0.5d0, 0.5d0)
            end if
        else                                        ! 接触なし
            normal_force_contact(particle_idx, wall_contact_slot_idx) = 0.0d0
            shear_force_contact(particle_idx, wall_contact_slot_idx) = 0.0d0
            contact_partner_idx(particle_idx, wall_contact_slot_idx) = 0
        end if
        
        ! 斜面壁 (摩擦斜面検証用) (contact_partner_idx = num_particles + 4)
        if (wall_validation_mode .and. wall_validation_type == 2) then
            wall_contact_slot_idx = 10 ! 斜面壁用の固定スロット
            wall_partner_id = num_particles + 4
            
            ! 斜面方程式: z = tan(angle) * x + 0 (原点を通る斜面)
            ! 粒子中心から斜面への距離
            slope_normal_x = -sin(wall_validation_slope_angle)
            slope_normal_z = cos(wall_validation_slope_angle)
            slope_distance = slope_normal_x * xi + slope_normal_z * zi
            
            if (slope_distance < ri_particle) then ! 斜面と接触
                wall_angle_sin = slope_normal_z
                wall_angle_cos = slope_normal_x
                overlap_gap = ri_particle - slope_distance
                contact_partner_idx(particle_idx, wall_contact_slot_idx) = wall_partner_id
                call actf_sub(particle_idx, wall_partner_id, wall_contact_slot_idx, wall_angle_sin, wall_angle_cos, overlap_gap, wall_validation_restitution_coeff, wall_validation_restitution_coeff)
            else                                    ! 接触なし
                normal_force_contact(particle_idx, wall_contact_slot_idx) = 0.0d0
                shear_force_contact(particle_idx, wall_contact_slot_idx) = 0.0d0
                contact_partner_idx(particle_idx, wall_contact_slot_idx) = 0
            end if
        end if
    end subroutine wcont_sub

    !> 粒子iと他の粒子との接触力を計算するサブルーチン
    subroutine pcont_sub(particle_i_idx, rmax_val)
        use simulation_constants_mod, only: nj_max
        use particle_data_mod
        use cell_system_mod
        implicit none

        integer, intent(in) :: particle_i_idx     ! 対象の粒子iのインデックス
        real(8), intent(in) :: rmax_val           ! 最大粒子半径 (探索範囲に使用)

        real(8) :: xi, zi, ri_particle_i
        real(8) :: xj, zj, rj_particle_j
        real(8) :: center_dist, overlap_gap
        real(8) :: contact_angle_sin, contact_angle_cos ! 粒子iから粒子jへの法線ベクトル成分
        integer :: particle_j_idx                 ! 接触相手の粒子jのインデックス
        integer :: contact_slot_for_i_j, contact_slot_for_j_i ! 接触リストのスロット
        integer :: iz_cell_min, iz_cell_max, ix_cell_min, ix_cell_max ! セル探索範囲
        integer :: current_iz_cell, current_ix_cell, cell_block_idx
        integer :: jj, max_particle_contacts_check
        real(8) :: dx, dz               ! 位置差計算用
        integer :: curr_j_idx
        real(8) :: search_extent        ! 近傍探索半径

        max_particle_contacts_check = 10 ! 元のコードのループ(do 11 jj=1,10)から、粒子間接触は最大10個と仮定

        xi = x_coord(particle_i_idx)
        zi = z_coord(particle_i_idx)
        ri_particle_i = radius(particle_i_idx)

        ! セル格子内での探索範囲を決定
        if (cell_size <= 0.0d0) stop "pcont_sub: cell_sizeが正ではありません。"

        ! 以前は ±(2*rmax) で探索していたが、cell_size が大きい場合に隣接セル2つ分を
        ! またぐ衝突を取りこぼすことがあった。そこで探索半径を (2*rmax + cell_size) に拡張する。
        search_extent = 2.0d0 * rmax_val + cell_size

        iz_cell_max = idint((zi + search_extent) / cell_size)
        iz_cell_min = idint((zi - search_extent) / cell_size)
        ix_cell_min = idint((xi - search_extent) / cell_size)
        ix_cell_max = idint((xi + search_extent) / cell_size)

        ! セルインデックスを有効なグリッド境界内に収める
        iz_cell_max = min(iz_cell_max, cells_z_dir - 1)
        iz_cell_min = max(iz_cell_min, 0)
        ix_cell_min = max(ix_cell_min, 0)
        ix_cell_max = min(ix_cell_max, cells_x_dir - 1)
        
        if (iz_cell_max < iz_cell_min .or. ix_cell_max < ix_cell_min) then
             ! 粒子が通常の領域外にある場合や、rmax_valが小さすぎて探索範囲が無効になる場合に発生しうる
             return
        end if

        ! デバッグ出力 (検証モードで粒子1のみ)
        if (validation_mode .and. particle_i_idx == 1) then
            write(*,'(A,I2,A,I3,A,I3,A,I3,A,I3)') 'Search p',particle_i_idx,': cells x[',ix_cell_min,':',ix_cell_max,'] z[',iz_cell_min,':',iz_cell_max,']'
            write(*,'(A,ES12.5,A,ES12.5,A,ES12.5)') 'Position: x=',xi,', z=',zi,', search_extent=',search_extent
        end if

        do current_iz_cell = iz_cell_min, iz_cell_max      ! z方向のセルループ
            do current_ix_cell = ix_cell_min, ix_cell_max  ! x方向のセルループ
                if (cells_x_dir <=0) stop "pcont_sub: cells_x_dirが正ではありません。"
                cell_block_idx = current_iz_cell * cells_x_dir + current_ix_cell + 1

                if (cell_block_idx <= 0 .or. cell_block_idx > (cells_x_dir * cells_z_dir) ) cycle ! セルが範囲外ならスキップ

                curr_j_idx = cell_head(cell_block_idx) ! セルに属する最初の粒子

                do while (curr_j_idx > 0)
                    particle_j_idx = curr_j_idx

                    if (particle_j_idx == particle_i_idx) then
                        curr_j_idx = particle_cell_next(curr_j_idx)
                        cycle
                    end if

                    xj = x_coord(particle_j_idx)
                    zj = z_coord(particle_j_idx)
                    rj_particle_j = radius(particle_j_idx)

                    dx = xi - xj
                    dz = zi - zj
                    center_dist = sqrt(dx*dx + dz*dz) 
                    overlap_gap = (ri_particle_i + rj_particle_j) - center_dist ! 重なり量 (正なら接触)

                    if (overlap_gap > 0.0d0) then ! 粒子が接触している (重なっている)
                        ! デバッグ出力 (検証モードのみ)
                        if (validation_mode) then
                            write(*,'(A,I2,A,I2,A,ES12.5,A,ES12.5,A,ES12.5)') &
                                'Contact p',particle_i_idx,' <-> p',particle_j_idx,': dist=',center_dist,', overlap=',overlap_gap,', sumr=',ri_particle_i+rj_particle_j
                        end if

                        if (center_dist < 1.0d-12) then ! 粒子中心が完全に一致する場合のゼロ除算を回避
                            contact_angle_cos = 1.0d0   ! 暫定的にx軸方向とする
                            contact_angle_sin = 0.0d0
                        else
                            ! 法線ベクトル (i から j へ向かう方向) の成分
                            contact_angle_cos = (xj - xi) / center_dist 
                            contact_angle_sin = (zj - zi) / center_dist
                        end if
                        
                        ! ---- 連絡スロット確保 (i -> j) ----------------
                        contact_slot_for_i_j = 0
                        do jj = 1, max_particle_contacts_check
                            if (contact_partner_idx(particle_i_idx, jj) == particle_j_idx) then
                                contact_slot_for_i_j = jj
                                exit
                            end if
                        end do
                        if (contact_slot_for_i_j == 0) then
                            do jj = 1, max_particle_contacts_check
                                if (contact_partner_idx(particle_i_idx, jj) == 0) then
                                    contact_slot_for_i_j = jj
                                    contact_partner_idx(particle_i_idx, jj) = particle_j_idx
                                    exit
                                end if
                            end do
                        end if
                        if (contact_slot_for_i_j == 0) then
                            curr_j_idx = particle_cell_next(curr_j_idx)
                            cycle  ! スロット不足
                        end if

                        ! -------------------------------------------------

                        ! 実際の力計算
                        if (validation_mode) then
                            call actf_sub(particle_i_idx, particle_j_idx, contact_slot_for_i_j, &
                                          contact_angle_sin, contact_angle_cos, overlap_gap, 1.0d0, 1.0d0)
                        else
                            call actf_sub(particle_i_idx, particle_j_idx, contact_slot_for_i_j, &
                                          contact_angle_sin, contact_angle_cos, overlap_gap, 0.5d0, 0.5d0)
                        end if
                    else ! 幾何学的な接触なし / 粒子が離れた
                        ! デバッグ出力 (検証モードのみ、最初の数回のみ)
                        if (validation_mode .and. particle_i_idx == 1 .and. particle_j_idx == 2) then
                            write(*,'(A,I2,A,I2,A,ES12.5,A,ES12.5,A,ES12.5)') &
                                'No contact p',particle_i_idx,' <-> p',particle_j_idx,': dist=',center_dist,', overlap=',overlap_gap,', sumr=',ri_particle_i+rj_particle_j
                        end if

                        ! particle_i_idx の particle_j_idx に関する接触情報をクリア
                        do jj = 1, max_particle_contacts_check
                            if (contact_partner_idx(particle_i_idx, jj) == particle_j_idx) then
                                normal_force_contact(particle_i_idx, jj) = 0.0d0
                                shear_force_contact(particle_i_idx, jj) = 0.0d0
                                contact_partner_idx(particle_i_idx, jj) = 0
                                exit
                            end if
                        end do
                    end if

                    curr_j_idx = particle_cell_next(curr_j_idx) ! セル内の次の粒子へ
                end do !! セル内連結リスト走査
            end do ! x方向セルループ
        end do     ! z方向セルループ
    end subroutine pcont_sub

    !> 粒子の位置と速度を更新するサブルーチン
    subroutine nposit_sub(judge_static)
        use simulation_constants_mod, only: GRAVITY_ACCEL
        use simulation_parameters_mod, only: time_step, validation_mode
        use particle_data_mod
        use cell_system_mod, only: num_particles
        implicit none
        integer, intent(out) :: judge_static ! 出力: 静止判定フラグ (1なら静止)

        integer :: i
        real(8) :: sum_abs_disp, avg_abs_disp
        real(8) :: grav

        if (validation_mode) then
            grav = 0.0d0
        else if (wall_validation_mode) then
            grav = GRAVITY_ACCEL
        else
            grav = GRAVITY_ACCEL
        end if

        sum_abs_disp = 0.0d0
        do i = 1, num_particles
            ! 速度の更新 (オイラー積分)
            z_vel(i) = z_vel(i) + (z_force_sum(i)/mass(i) - grav) * time_step
            x_vel(i) = x_vel(i) + (x_force_sum(i)/mass(i)) * time_step
            rotation_vel(i) = rotation_vel(i) + (moment_sum(i)/moment_inertia(i)) * time_step

            ! 現ステップの変位増分を更新 (元のコードの修正オイラー/平均化スキーム)
            ! v(i)=(v0(i)*dt+v(i))/2.d0 ここで右辺のv(i)は前ステップの変位増分。
            z_disp_incr(i) = (z_vel(i) * time_step + z_disp_incr(i)) / 2.0d0
            x_disp_incr(i) = (x_vel(i) * time_step + x_disp_incr(i)) / 2.0d0
            rot_disp_incr(i) = (rotation_vel(i) * time_step + rot_disp_incr(i)) / 2.0d0
            
            ! 位置と回転角の更新
            z_coord(i) = z_coord(i) + z_disp_incr(i)
            x_coord(i) = x_coord(i) + x_disp_incr(i)
            rotation_angle(i) = rotation_angle(i) + rot_disp_incr(i)

            sum_abs_disp = sum_abs_disp + abs(x_disp_incr(i)) + abs(z_disp_incr(i))
        end do

        ! 静止状態の判定
        if (num_particles > 0) then
            avg_abs_disp = sum_abs_disp / real(num_particles, 8) / 2.0d0 ! 元のコードの /2.d0
            if (avg_abs_disp < (time_step * time_step * GRAVITY_ACCEL * 1.0d-1)) then
                judge_static = 1
            else
                judge_static = 0
            end if
        else
            judge_static = 1 ! 粒子がなければ静止とみなす
        end if
    end subroutine nposit_sub

    !> 粒子iと粒子/壁jとの間の実際の接触力（法線方向およびせん断方向）を計算するサブルーチン
    subroutine actf_sub(p_i, p_j, contact_slot_idx_for_pi, angle_sin, angle_cos, initial_overlap, en_coeff, et_coeff)
        ! p_i: 主となる粒子のインデックス
        ! p_j: 他方の粒子インデックス (<= num_particles の場合) または壁ID (> num_particles の場合)
        ! contact_slot_idx_for_pi: p_i の接触配列における p_j のスロット
        ! angle_sin, angle_cos: p_i から p_j の中心/接触点への法線ベクトル成分
        ! initial_overlap: 幾何学的な重なり量、接触していれば正
        ! en_coeff, et_coeff: 反発係数 (normal, tangential)
        use simulation_constants_mod, only: ni_max
        use simulation_parameters_mod, only: time_step, validation_mode
        use particle_data_mod
        use cell_system_mod, only: num_particles
        implicit none

        integer, intent(in) :: p_i, p_j, contact_slot_idx_for_pi
        real(8), intent(in) :: angle_sin, angle_cos, initial_overlap
        real(8), intent(in) :: en_coeff, et_coeff ! 反発係数 (normal, tangential)

        real(8) :: ri_val, rj_val, effective_mass
        real(8) :: kn_normal_stiffness, ks_shear_stiffness ! 法線・せん断バネ定数 Kn, Ks
        real(8) :: damping_coeff_normal, damping_coeff_shear ! 法線・せん断粘性係数 ηn, ηs
        real(8) :: rel_disp_normal_incr, rel_disp_shear_incr ! 法線・せん断方向の相対変位増分 Δun, Δus
        real(8) :: damping_force_normal, damping_force_shear ! 法線・せん断方向の粘性抵抗力 dn, ds
        real(8) :: total_normal_force, total_shear_force     ! 全法線力 fn, 全せん断力 fs
        real(8) :: friction_coeff_current      ! 現在の摩擦係数 μ
        real(8) :: critical_time_step_check    ! 安定性チェック用の時間刻み (ddt)

        real(8) :: x_disp_incr_pi, z_disp_incr_pi, rot_disp_incr_pi ! 粒子iの変位増分
        real(8) :: x_disp_incr_pj, z_disp_incr_pj, rot_disp_incr_pj ! 粒子j(または壁=0)の変位増分
        real(8) :: mass_pi, mass_pj                                 ! 粒子i,jの質量
        real(8) :: r_eff ! 宣言をここに移動

        ri_val = radius(p_i)
        x_disp_incr_pi = x_disp_incr(p_i)
        z_disp_incr_pi = z_disp_incr(p_i)
        rot_disp_incr_pi = rot_disp_incr(p_i)
        mass_pi = mass(p_i)

        if (p_j <= num_particles) then ! 粒子間
            rj_val = radius(p_j)
            x_disp_incr_pj = x_disp_incr(p_j)
            z_disp_incr_pj = z_disp_incr(p_j)
            rot_disp_incr_pj = rot_disp_incr(p_j)
            mass_pj = mass(p_j)
            if (mass_pi + mass_pj > 1.0d-20) then
                 effective_mass = mass_pi * mass_pj / (mass_pi + mass_pj) ! 等価質量 m_eff = m1*m2/(m1+m2)
            else
                 effective_mass = mass_pi ! フォールバック
            end if
            friction_coeff_current = friction_coeff_particle
            r_eff = ri_val * rj_val / (ri_val + rj_val) ! 等価半径 Reff = ri*rj/(ri+rj)
            
            ! Hertz接触理論に基づく法線剛性 (粒子間)
            kn_normal_stiffness = (4.0d0/3.0d0) * sqrt(r_eff) * &
                                 young_modulus_particle / (1.0d0 - poisson_ratio_particle**2) * &
                                 sqrt(initial_overlap)
            if (kn_normal_stiffness < 1.0d6) kn_normal_stiffness = 1.0d6 ! 最小値設定
        else ! 粒子-壁接触
            rj_val = 0.0d0 ! 壁の半径は実質無限大、または変位にrjは使用しない
            x_disp_incr_pj = 0.0d0 ! 壁は動かないと仮定
            z_disp_incr_pj = 0.0d0
            rot_disp_incr_pj = 0.0d0
            effective_mass = mass_pi   ! 等価質量 m_eff = m1
            friction_coeff_current = friction_coeff_wall
            r_eff = ri_val ! 粒子-壁接触では粒子半径を使用
            
            ! Hertz接触理論に基づく法線剛性 (粒子-壁)
            kn_normal_stiffness = (4.0d0/3.0d0) * sqrt(r_eff) * &
                                 young_modulus_particle * young_modulus_wall / &
                                 ((1.0d0-poisson_ratio_particle**2)*young_modulus_wall + &
                                  (1.0d0-poisson_ratio_wall**2)*young_modulus_particle) * &
                                 sqrt(initial_overlap)
            if (kn_normal_stiffness < 1.0d6) kn_normal_stiffness = 1.0d6 ! 最小値設定
        end if
        ks_shear_stiffness = kn_normal_stiffness * shear_to_normal_stiffness_ratio

        ! 完全弾性衝突 (e=1) のための粘性係数設定
        if (validation_mode .and. en_coeff >= 0.99d0) then
            ! 完全弾性衝突の場合、粘性をほぼゼロに設定
            damping_coeff_normal = 0.0d0
            damping_coeff_shear = 0.0d0
        else
            ! 粘性係数 (反発係数に基づいて計算)
            if (effective_mass > 0.0d0 .and. kn_normal_stiffness > 0.0d0 .and. en_coeff > 1.0d-6) then
                damping_coeff_normal = -2.0d0 * log(en_coeff) * sqrt(effective_mass * kn_normal_stiffness / (log(en_coeff)**2 + PI_VAL**2))
            else
                damping_coeff_normal = 0.0d0
            end if
            if (effective_mass > 0.0d0 .and. ks_shear_stiffness > 0.0d0 .and. et_coeff > 1.0d-6) then
                damping_coeff_shear = -2.0d0 * log(et_coeff) * sqrt(effective_mass * ks_shear_stiffness / (log(et_coeff)**2 + PI_VAL**2))
            else
                damping_coeff_shear = 0.0d0
            end if
        end if
        
        ! 安定性基準のチェック (元のddt、レイリー時間刻みに関連)
        if (kn_normal_stiffness > 1.0d-12) then
           critical_time_step_check = 0.1d0 * sqrt(effective_mass / kn_normal_stiffness)
           if (critical_time_step_check < time_step .and. critical_time_step_check > 1.0d-12) then 
                write(*,*) '警告: 安定性基準違反 - 推奨時間刻み:', critical_time_step_check, '現在:', time_step
                write(*,*) '  Kn=', kn_normal_stiffness, ' M_eff=', effective_mass, ' 粒子:', p_i, p_j
           end if
        end if

        ! 相対変位増分 (現時間ステップ time_step における)
        ! angle成分は粒子p_iからp_jへの法線ベクトルを定義
        rel_disp_normal_incr = (x_disp_incr_pi - x_disp_incr_pj) * angle_cos + &
                               (z_disp_incr_pi - z_disp_incr_pj) * angle_sin
        rel_disp_shear_incr = -(x_disp_incr_pi - x_disp_incr_pj) * angle_sin + &
                               (z_disp_incr_pi - z_disp_incr_pj) * angle_cos + &
                               (ri_val * rot_disp_incr_pi + rj_val * rot_disp_incr_pj)


        ! 弾性力成分の更新 (式3.7, 3.10)
        if (abs(normal_force_contact(p_i, contact_slot_idx_for_pi)) < 1.0d-8) then ! 新規接触
            ! 新規接触の場合、弾性力を重なり量に基づいて初期化
            normal_force_contact(p_i, contact_slot_idx_for_pi) = kn_normal_stiffness * initial_overlap
            shear_force_contact(p_i, contact_slot_idx_for_pi) = 0.0d0 ! せん断力は初期化時にゼロ
        else
            ! 既存接触の場合、増分で更新
            normal_force_contact(p_i, contact_slot_idx_for_pi) = normal_force_contact(p_i, contact_slot_idx_for_pi) + &
                                                                 kn_normal_stiffness * rel_disp_normal_incr
            shear_force_contact(p_i, contact_slot_idx_for_pi) = shear_force_contact(p_i, contact_slot_idx_for_pi) + &
                                                                ks_shear_stiffness * rel_disp_shear_incr
        end if
        
        ! 粘性抵抗力成分の計算 (式3.6, 3.9)
        if (time_step > 1.0d-20) then
             damping_force_normal = damping_coeff_normal * rel_disp_normal_incr / time_step
             damping_force_shear = damping_coeff_shear * rel_disp_shear_incr / time_step
        else
             damping_force_normal = 0.0d0
             damping_force_shear  = 0.0d0
        end if

        ! 引張力のチェック (粒子が引き離される場合) - 付着力は考慮しない
        if (normal_force_contact(p_i, contact_slot_idx_for_pi) < 0.0d0) then
            normal_force_contact(p_i, contact_slot_idx_for_pi) = 0.0d0
            shear_force_contact(p_i, contact_slot_idx_for_pi) = 0.0d0
            damping_force_normal = 0.0d0
            damping_force_shear = 0.0d0
            contact_partner_idx(p_i, contact_slot_idx_for_pi) = 0 ! 引張なら接触を切る
            return ! 引き離される場合は力なし
        end if

        ! クーロンの摩擦法則を適用 (式3.11)
        if (abs(shear_force_contact(p_i, contact_slot_idx_for_pi)) > &
            friction_coeff_current * normal_force_contact(p_i, contact_slot_idx_for_pi)) then
            shear_force_contact(p_i, contact_slot_idx_for_pi) = friction_coeff_current * &
                normal_force_contact(p_i, contact_slot_idx_for_pi) * &
                sign(1.0d0, shear_force_contact(p_i, contact_slot_idx_for_pi))
            damping_force_shear = 0.0d0 ! 滑りが発生している場合はせん断粘性なし
        end if

        ! 粘性を含む合計の力 (式3.8, 3.12)
        total_normal_force = normal_force_contact(p_i, contact_slot_idx_for_pi) + damping_force_normal
        total_shear_force = shear_force_contact(p_i, contact_slot_idx_for_pi) + damping_force_shear

        ! 粒子p_iに力を適用 (式3.13)
        ! 法線力は中心を結ぶ線に沿って作用 (angle_cos, angle_sin で定義される iからjへの方向)
        ! せん断力はそれに垂直。
        x_force_sum(p_i) = x_force_sum(p_i) - total_normal_force * angle_cos + total_shear_force * angle_sin
        z_force_sum(p_i) = z_force_sum(p_i) - total_normal_force * angle_sin - total_shear_force * angle_cos
        moment_sum(p_i) = moment_sum(p_i) - ri_val * total_shear_force

        ! 粒子p_jに反作用力を適用 (相手が粒子の場合)
        if (p_j <= num_particles .and. contact_slot_idx_for_pi <= 10) then ! 元の jk < 10 は粒子間接触のチェック
            x_force_sum(p_j) = x_force_sum(p_j) + total_normal_force * angle_cos - total_shear_force * angle_sin
            z_force_sum(p_j) = z_force_sum(p_j) + total_normal_force * angle_sin + total_shear_force * angle_cos
            moment_sum(p_j) = moment_sum(p_j) - rj_val * total_shear_force ! p_jに対するせん断力は同じ大きさ、逆向きの回転効果
            
            ! デバッグ出力 (検証モードのみ)
            if (validation_mode .and. abs(total_normal_force) > 1.0d-6) then
                write(*,'(A,I2,A,I2,A,ES12.5,A,ES12.5)') 'Force p',p_i,' -> p',p_j,': Fn=',total_normal_force,', Fs=',total_shear_force
                write(*,'(A,ES12.5,A,ES12.5)') '  overlap=',initial_overlap,', dist=',sqrt((x_coord(p_i)-x_coord(p_j))**2 + (z_coord(p_i)-z_coord(p_j))**2)
            end if
        end if
    end subroutine actf_sub

    !> グラフィック用データを出力するサブルーチン
    subroutine gfout_sub(iter_step, time_val, rmax_val)
        use simulation_constants_mod, only: nj_max
        use simulation_parameters_mod, only: container_width
        use particle_data_mod
        use cell_system_mod, only: num_particles
        implicit none
        integer, intent(in) :: iter_step    ! 現在のイテレーションステップ
        real(8), intent(in) :: time_val, rmax_val ! 現在時刻、最大粒子半径
        integer :: i,j

        if (iter_step == 1) then
            open(unit=10, file='data/graph11.d', status='replace', action='write')
            open(unit=11, file='data/graph21.d', status='replace', action='write')
        end if
        ! 初回以降は追記モードで開くか、性能が許せば開いたままにする。
        ! 簡単のため、ここでは毎回開く(replace)。元のコードはファイルを開きっぱなしにするか、適切に再オープン。

        write(10,*) num_particles, time_val, container_width, rmax_val
        if (num_particles > 0) then
            write(10,'(1000(ES12.5,1X,ES12.5,1X,ES12.5,2X))') (sngl(x_coord(i)), sngl(z_coord(i)), sngl(radius(i)), i=1,num_particles)
            write(10,'(1000(ES12.5,1X,ES12.5,1X,ES12.5,2X))') (sngl(x_vel(i)), sngl(z_vel(i)), sngl(rotation_vel(i)), i=1,num_particles)
            write(10,'(1000(ES12.5,2X))') (sngl(rotation_angle(i)), i=1,num_particles)
        end if
        
        ! 接触力の出力 (オプション、graph21.dより)
        write(11,*) 'Time: ', time_val 
        if (num_particles > 0) then
            do i = 1, num_particles
                 write(11,'(A,I5,A,I5)') 'Particle: ', i, ' NumContacts: ', count(contact_partner_idx(i,1:nj_max) > 0)
                 write(11,'(2X,A,13(ES10.3,1X))') 'ShearF: ', (shear_force_contact(i,j), j=1,nj_max)
                 write(11,'(2X,A,13(ES10.3,1X))') 'NormalF:', (normal_force_contact(i,j), j=1,nj_max)
                 write(11,'(2X,A,13(I5,2X))')    'Partner:', (contact_partner_idx(i,j), j=1,nj_max)
            end do
        end if

        ! 元のコードではここでファイルを閉じない。プログラム終了時に閉じられることを示唆。
    end subroutine gfout_sub

    !> バックアップデータを出力するサブルーチン
    subroutine bfout_sub
        use simulation_constants_mod
        use simulation_parameters_mod
        use particle_data_mod
        use cell_system_mod
        implicit none
        integer :: i, j
        real(8) :: rmax_dummy_val ! 元のbfoutはrmaxを必要とするが、メインの呼び出しからは渡されない。
                                  ! リストアに不可欠でないか、粒子半径から取得すると仮定。
        
        if (num_particles > 0) then
           rmax_dummy_val = maxval(radius(1:num_particles))
        else
           rmax_dummy_val = 0.0d0
        end if

        open(unit=13, file='data/backl.d', status='replace', action='write')

        write(13,*) num_particles, cells_x_dir, cells_z_dir, particle_gen_layers
        write(13,*) rmax_dummy_val, 0.0d0, container_width, cell_size, time_step ! current_timeではなく初期t=0を保存すると仮定
        write(13,*) particle_density, friction_coeff_particle, friction_coeff_wall, GRAVITY_ACCEL, PI_VAL
        write(13,*) young_modulus_particle, young_modulus_wall, poisson_ratio_particle, poisson_ratio_wall, shear_to_normal_stiffness_ratio
        
        if (num_particles > 0) then
            write(13,*) (mass(i), moment_inertia(i), i=1,num_particles)
            write(13,*) (x_coord(i), z_coord(i), radius(i), i=1,num_particles)
            write(13,*) (x_disp_incr(i), z_disp_incr(i), rot_disp_incr(i), i=1,num_particles) ! u,v,f (dpm)
            write(13,*) (x_vel(i), z_vel(i), rotation_vel(i), i=1,num_particles)            ! u0,v0,f0
            do i = 1, num_particles
                write(13,*) (shear_force_contact(i,j), normal_force_contact(i,j), j=1,nj_max)
                write(13,*) (contact_partner_idx(i,j), j=1,nj_max)
            end do
        end if
        close(13)
    end subroutine bfout_sub

    !> 擬似乱数を生成するサブルーチン
    subroutine custom_random(seed_io, random_val_out)
        implicit none
        integer, intent(inout) :: seed_io          ! ジェネレータの現在の状態を保持
        real(8), intent(out)   :: random_val_out   ! [0,1) の一様乱数

        ! 元のコードの定数とロジックを可能な限り再現
        seed_io = seed_io * 65539 
        if (seed_io < 0) then
             seed_io = (seed_io + 2147483647) + 1 
        end if
        random_val_out = dble(seed_io) * 0.4656613d-9 ! 元の正規化定数 (1.0 / 2147483648.0)

    end subroutine custom_random

    !> パラメータ掃引による自由落下反発検証を実行するサブルーチン
    subroutine parameter_sweep_validation
        implicit none
        integer :: sweep_idx
        real(8) :: current_restitution_coeff, restitution_step
        real(8) :: theoretical_height, actual_height, height_error_percent
        
        write(*,*) '================================='
        write(*,*) 'パラメータ掃引による自由落下反発検証'
        write(*,*) '================================='
        write(*,*) '反発係数範囲: ', restitution_coeff_min, ' - ', restitution_coeff_max
        write(*,*) '掃引回数: ', parameter_sweep_count
        write(*,*) '================================='
        
        ! 結果出力ファイルを開く
        open(unit=30, file='data/parameter_sweep_results.csv', status='replace', action='write')
        write(30,*) 'Restitution_Coefficient,Theoretical_Height,Actual_Height,Error_Percent'
        
        ! 反発係数のステップサイズを計算
        if (parameter_sweep_count > 1) then
            restitution_step = (restitution_coeff_max - restitution_coeff_min) / real(parameter_sweep_count - 1, 8)
        else
            restitution_step = 0.0d0
        end if
        
        do sweep_idx = 1, parameter_sweep_count
            ! 現在の反発係数を設定
            current_restitution_coeff = restitution_coeff_min + real(sweep_idx - 1, 8) * restitution_step
            wall_validation_restitution_coeff = current_restitution_coeff
            
            write(*,*) ''
            write(*,*) '--- 掃引 ', sweep_idx, '/', parameter_sweep_count, ' ---'
                    write(*,*) '反発係数: ', current_restitution_coeff
        
        ! 単一の検証シミュレーションを実行
        call single_validation_run(theoretical_height, actual_height, height_error_percent)
        
        ! デバッグ出力
        write(*,*) '  -> 理論値:', theoretical_height, ', 実際値:', actual_height
            
            ! 結果をファイルに出力
            write(30,'(F6.3,A,F10.6,A,F10.6,A,F8.4)') current_restitution_coeff, ',', &
                theoretical_height, ',', actual_height, ',', height_error_percent
            
            write(*,*) '理論反発高さ: ', theoretical_height
            write(*,*) '計算反発高さ: ', actual_height
            write(*,*) '相対誤差: ', height_error_percent, '%'
        end do
        
        close(30)
        write(*,*) ''
        write(*,*) '================================='
        write(*,*) 'パラメータ掃引完了'
        write(*,*) '結果は data/parameter_sweep_results.csv に保存されました'
        write(*,*) '================================='
    end subroutine parameter_sweep_validation

    !> 単一の自由落下反発検証シミュレーションを実行するサブルーチン
    subroutine single_validation_run(theoretical_height_out, actual_height_out, error_percent_out)
        implicit none
        real(8), intent(out) :: theoretical_height_out, actual_height_out, error_percent_out
        
        integer :: it_step, static_judge_flag
        real(8) :: current_time, rmax_particle_radius
        logical :: wall_collision_started_local, wall_collision_finished_local
        logical :: particle_falling_local, particle_rising_local
        real(8) :: initial_drop_height_local, pre_collision_velocity_local
        real(8) :: max_rebound_height_local, theoretical_rebound_height_local
        real(8) :: height_error  ! 変数定義を追加
        integer :: i
        
        ! 初期化
        call fposit_sub(rmax_particle_radius)
        call inmat_sub
        call init_sub
        
        ! 検証用変数の初期化
        wall_collision_started_local = .false.
        wall_collision_finished_local = .false.
        particle_falling_local = .true.
        particle_rising_local = .false.
        initial_drop_height_local = wall_validation_drop_height
        theoretical_rebound_height_local = wall_validation_restitution_coeff**2 * initial_drop_height_local
        max_rebound_height_local = 0.0d0  ! 最高点の初期化
        
        ! 粒子-壁間検証モードの初期速度設定
        x_vel(1) = 0.0d0
        z_vel(1) = 0.0d0
        rotation_vel(1) = 0.0d0
        friction_coeff_particle = 0.0d0
        friction_coeff_wall = 0.0d0
        
        current_time = 0.0d0
        
        ! シミュレーションループ
        do it_step = 1, max_calculation_steps
            current_time = current_time + time_step

            call ncel_sub

            do i = 1, num_particles
                x_force_sum(i) = 0.0d0
                z_force_sum(i) = 0.0d0
                moment_sum(i) = 0.0d0
            end do
            
            do i = 1, num_particles
                call wcont_sub(i)
                call pcont_sub(i, rmax_particle_radius)
            end do

            call nposit_sub(static_judge_flag)

            ! 自由落下反発検証ロジック
            if (.not. wall_collision_finished_local) then
                if (particle_falling_local .and. z_coord(1) <= radius(1) + 1.0d-6) then
                    ! 床に衝突
                    wall_collision_started_local = .true.
                    pre_collision_velocity_local = abs(z_vel(1))
                    write(*,*) '床衝突検出: 時刻 = ', current_time
                    write(*,*) '衝突直前速度: ', pre_collision_velocity_local
                else if (wall_collision_started_local .and. particle_falling_local .and. z_coord(1) > radius(1) + 1.0d-3) then
                    ! 床から離れた（反発開始）
                    particle_falling_local = .false.
                    particle_rising_local = .true.
                    write(*,*) '床から離脱: 時刻 = ', current_time, ', 位置 = ', z_coord(1)
                else if (particle_rising_local .and. z_vel(1) < 0.0d0) then
                    ! 上昇から落下に転じた（最高点到達）
                    particle_rising_local = .false.
                    particle_falling_local = .true.
                    max_rebound_height_local = z_coord(1)
                    wall_collision_finished_local = .true.
                    
                    ! 理論値との比較
                    height_error = abs(max_rebound_height_local - theoretical_rebound_height_local) / theoretical_rebound_height_local * 100.0d0
                    
                    write(*,*) '最高点検出: 時刻 = ', current_time, ', 位置 = ', z_coord(1), ', 速度 = ', z_vel(1)
                    write(*,*) '================================='
                    write(*,*) '自由落下反発検証結果'
                    write(*,*) '================================='
                    write(*,*) '初期落下高さ: ', initial_drop_height_local
                    write(*,*) '反発係数: ', wall_validation_restitution_coeff
                    write(*,*) '理論反発高さ: ', theoretical_rebound_height_local
                    write(*,*) '計算反発高さ: ', max_rebound_height_local
                    write(*,*) '相対誤差: ', height_error, '%'
                    write(*,*) '================================='
                    
                    exit  ! シミュレーションループを抜ける
                end if
            end if

            if (static_judge_flag == 1) then
                exit
            end if
        end do
        
        ! 最高点検出が失敗した場合のフォールバック処理
        if (.not. wall_collision_finished_local) then
            if (particle_rising_local) then
                ! 上昇中にシミュレーションが終了した場合、現在の高さを最高点とする
                max_rebound_height_local = z_coord(1)
            else
                ! その他の場合、警告を出力
                write(*,*) '警告: 最高点検出が失敗しました'
                max_rebound_height_local = 0.0d0
            end if
        end if
        
        ! 結果を返す
        theoretical_height_out = theoretical_rebound_height_local
        actual_height_out = max_rebound_height_local
        if (theoretical_height_out > 0.0d0) then
            error_percent_out = abs(actual_height_out - theoretical_height_out) / theoretical_height_out * 100.0d0
        else
            error_percent_out = 0.0d0
        end if
    end subroutine single_validation_run

end program two_dimensional_pem