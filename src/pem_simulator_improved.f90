! 改良版メインプログラムおよびモジュール: 2次元粒子要素法 (球体モデル)
! dem_calc.pyを参考にした自由落下シミュレーション改良版

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
    
    ! 自由落下検証モード設定
    logical :: free_fall_mode                 ! 自由落下モードフラグ
    real(8) :: drop_height                    ! 落下開始高さ
    real(8) :: particle_x_pos                 ! 粒子のx座標
    real(8) :: particle_radius               ! 粒子半径
    real(8) :: restitution_coefficient       ! 反発係数
    
    ! 出力制御パラメータ
    integer :: output_interval                 ! 出力間隔
    integer :: max_calculation_steps          ! 最大計算ステップ数
    
    ! 詳細解析用パラメータ
    logical :: detailed_analysis              ! 詳細解析モードフラグ
    logical :: energy_analysis                ! エネルギー解析フラグ
    logical :: contact_analysis               ! 接触解析フラグ

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
    real(8), dimension(ni_max) :: z_coord        ! z0(ni): 粒子中心のz座標
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
    integer, dimension(ni_max, nj_max) :: contact_partner_idx ! je(ni,nj): 接触点番号配列

    ! 現時間ステップにおける増分変位
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
    integer :: cells_z_dir            ! idz: z方向のセル数

    real(8) :: cell_size              ! c: セルの幅/サイズ

    ! 各セルに属する粒子を連結リストで保持する
    integer, dimension(nc_max) :: cell_head        ! そのセルで最初に登録された粒子インデックス
    integer, dimension(ni_max) :: particle_cell_next ! 次の粒子インデックス

    integer, dimension(nc_max) :: cell_particle_map
    integer, dimension(ni_max) :: particle_cell_idx ! 粒子iが格納されているセル番号

    save
end module cell_system_mod

! モジュール: 解析用データ
module analysis_data_mod
    use simulation_constants_mod, only: ni_max
    implicit none

    ! 接触追跡用変数
    logical :: in_contact                         ! 接触状態フラグ
    integer :: contact_start_step                 ! 接触開始ステップ
    integer :: contact_end_step                   ! 接触終了ステップ
    real(8) :: max_overlap                        ! 最大オーバーラップ
    real(8) :: pre_collision_velocity             ! 衝突前速度
    real(8) :: post_collision_velocity            ! 衝突後速度
    
    ! エネルギー追跡用変数
    real(8) :: initial_potential_energy          ! 初期位置エネルギー
    real(8) :: initial_kinetic_energy            ! 初期運動エネルギー
    real(8) :: current_potential_energy          ! 現在の位置エネルギー
    real(8) :: current_kinetic_energy            ! 現在の運動エネルギー
    
    ! 自由落下解析用変数
    logical :: falling_phase                      ! 落下中フラグ
    logical :: rising_phase                       ! 上昇中フラグ
    logical :: collision_detected                 ! 衝突検出フラグ
    real(8) :: max_rebound_height                 ! 最大反発高さ
    real(8) :: theoretical_rebound_height         ! 理論的反発高さ
    
    save
end module analysis_data_mod

! メインプログラム
program free_fall_simulation
    use simulation_constants_mod
    use simulation_parameters_mod
    use particle_data_mod
    use cell_system_mod
    use analysis_data_mod
    implicit none

    integer :: it_step, static_judge_flag
    integer :: i
    real(8) :: current_time
    real(8) :: rmax_particle_radius
    
    ! 計算時間計測用変数
    integer :: start_time, end_time, clock_rate
    real(8) :: elapsed_time
    
    ! 理論解計算用変数
    real(8) :: theoretical_height, theoretical_velocity, theoretical_time
    real(8) :: height_error, velocity_error
    
    write(*,*) '========================================='
    write(*,*) '   改良版自由落下シミュレーション'
    write(*,*) '   (dem_calc.py準拠の詳細解析機能付き)'
    write(*,*) '========================================='
    
    ! 計算時間計測開始
    call system_clock(start_time, clock_rate)
    
    ! パラメータ読み込み
    call read_improved_input_file
    
    ! 出力ディレクトリの作成
    call create_output_directories
    
    ! 初期化
    call initialize_simulation(rmax_particle_radius)
    
    ! 解析データの初期化
    call initialize_analysis_data
    
    ! 初期状態の保存
    call save_state_to_csv(0, 0.0d0)
    
    write(*,*) '計算開始'
    write(*,*) '目標ステップ数: ', max_calculation_steps
    write(*,*) '時間刻み幅: ', time_step, ' s'
    write(*,*) '出力間隔: ', output_interval, ' ステップ'
    write(*,*) '-----------------------------------------'
    
    current_time = 0.0d0
    
    ! メインシミュレーションループ
    do it_step = 1, max_calculation_steps
        current_time = current_time + time_step

        ! セル更新
        call ncel_sub

        ! 力の初期化
        do i = 1, num_particles
            x_force_sum(i) = 0.0d0
            z_force_sum(i) = 0.0d0
            moment_sum(i) = 0.0d0
        end do
        
        ! 接触力計算
        do i = 1, num_particles
            call wcont_sub(i)
            call pcont_sub(i, rmax_particle_radius)
        end do

        ! 重力の追加
        do i = 1, num_particles
            z_force_sum(i) = z_force_sum(i) - mass(i) * GRAVITY_ACCEL
        end do

        ! 位置・速度更新
        call nposit_sub(static_judge_flag)

        ! 詳細解析の実行
        if (detailed_analysis) then
            call analyze_free_fall(it_step, current_time)
        end if

        ! 定期的な状態保存
        if (mod(it_step, output_interval) == 0) then
            call save_state_to_csv(it_step, current_time)
            call print_progress(it_step, max_calculation_steps, current_time)
        end if
        
        ! 終了条件チェック
        if (static_judge_flag == 1) then
            write(*,*) '静止状態に到達しました。'
            exit
        end if
        
        ! 自由落下モードでの終了判定
        if (free_fall_mode .and. collision_detected .and. .not. rising_phase .and. &
            abs(z_vel(1)) < 1.0d-6 .and. z_coord(1) < drop_height * 0.1d0) then
            write(*,*) '自由落下解析完了: 粒子が停止しました。'
            exit
        end if
    end do
    
    ! 計算時間計測終了
    call system_clock(end_time, clock_rate)
    elapsed_time = real(end_time - start_time, 8) / real(clock_rate, 8)
    
    ! 最終結果の出力
    call output_final_results(it_step, current_time, elapsed_time)
    
    ! 理論解との比較
    if (free_fall_mode) then
        call compare_with_theory
    end if
    
    write(*,*) '-----------------------------------------'
    write(*,*) '計算完了'
    write(*,*) '総計算時間: ', elapsed_time, ' 秒'
    write(*,*) '========================================='

contains

    !> 改良版入力ファイル読み込みサブルーチン
    subroutine read_improved_input_file
        implicit none
        character(len=256) :: line, keyword
        real(8) :: value
        integer :: ios
        
        ! デフォルト値の設定
        time_step = 1.0d-5
        particle_density = 2480.0d0
        container_width = 20.0d0
        max_calculation_steps = 100000
        output_interval = 100
        
        ! 自由落下モードのデフォルト設定
        free_fall_mode = .true.
        drop_height = 10.0d0
        particle_x_pos = 10.0d0
        particle_radius = 1.0d0
        restitution_coefficient = 0.8d0
        
        ! 解析モードのデフォルト設定
        detailed_analysis = .true.
        energy_analysis = .true.
        contact_analysis = .true.
        
        ! 材料物性値のデフォルト設定
        young_modulus_particle = 4.9d7
        young_modulus_wall = 4.9d7
        poisson_ratio_particle = 0.23d0
        poisson_ratio_wall = 0.23d0
        friction_coeff_particle = 0.0d0
        friction_coeff_wall = 0.0d0
        
        write(*,*) '入力ファイルを読み込み中...'
        
        open(unit=10, file='input/input_free_fall.dat', status='old', iostat=ios)
        if (ios /= 0) then
            write(*,*) '警告: input_free_fall.datが見つかりません。デフォルト値を使用します。'
            return
        end if
        
        do
            read(10, '(A)', iostat=ios) line
            if (ios /= 0) exit
            
            ! コメント行とブランク行をスキップ
            if (line(1:1) == '#' .or. len_trim(line) == 0) cycle
            
            read(line, *) keyword, value
            
            select case (trim(keyword))
                case ('TIME_STEP')
                    time_step = value
                case ('PARTICLE_DENSITY')
                    particle_density = value
                case ('CONTAINER_WIDTH')
                    container_width = value
                case ('MAX_CALCULATION_STEPS')
                    max_calculation_steps = int(value)
                case ('OUTPUT_INTERVAL')
                    output_interval = int(value)
                case ('DROP_HEIGHT')
                    drop_height = value
                case ('PARTICLE_X_POS')
                    particle_x_pos = value
                case ('PARTICLE_RADIUS')
                    particle_radius = value
                case ('RESTITUTION_COEFFICIENT')
                    restitution_coefficient = value
                case ('YOUNG_MODULUS_PARTICLE')
                    young_modulus_particle = value
                case ('YOUNG_MODULUS_WALL')
                    young_modulus_wall = value
                case ('POISSON_RATIO_PARTICLE')
                    poisson_ratio_particle = value
                case ('POISSON_RATIO_WALL')
                    poisson_ratio_wall = value
                case ('FRICTION_COEFF_PARTICLE')
                    friction_coeff_particle = value
                case ('FRICTION_COEFF_WALL')
                    friction_coeff_wall = value
            end select
        end do
        
        close(10)
        write(*,*) '入力ファイル読み込み完了'
    end subroutine read_improved_input_file

    !> 出力ディレクトリ作成サブルーチン
    subroutine create_output_directories
        implicit none
        
        ! データディレクトリの作成（システムコマンド使用）
        call system('mkdir -p data')
        call system('mkdir -p data/detailed')
        call system('mkdir -p data/analysis')
        
        write(*,*) '出力ディレクトリを作成しました'
    end subroutine create_output_directories

    !> シミュレーション初期化サブルーチン
    subroutine initialize_simulation(rmax_radius)
        implicit none
        real(8), intent(out) :: rmax_radius
        
        ! 粒子数を1に設定（自由落下）
        num_particles = 1
        
        ! 粒子の物理特性設定
        radius(1) = particle_radius
        mass(1) = (4.0d0/3.0d0) * PI_VAL * particle_density * radius(1)**3
        moment_inertia(1) = (2.0d0/5.0d0) * mass(1) * radius(1)**2
        
        ! 初期位置設定
        x_coord(1) = particle_x_pos
        z_coord(1) = drop_height
        rotation_angle(1) = 0.0d0
        
        ! 初期速度設定（静止状態から開始）
        x_vel(1) = 0.0d0
        z_vel(1) = 0.0d0
        rotation_vel(1) = 0.0d0
        
        ! 力とモーメントの初期化
        x_force_sum(1) = 0.0d0
        z_force_sum(1) = 0.0d0
        moment_sum(1) = 0.0d0
        
        ! 接触力の初期化
        do i = 1, nj_max
            normal_force_contact(1, i) = 0.0d0
            shear_force_contact(1, i) = 0.0d0
            contact_partner_idx(1, i) = 0
        end do
        
        rmax_radius = radius(1)
        
        write(*,*) 'シミュレーション初期化完了'
        write(*,*) '粒子質量: ', mass(1), ' kg'
        write(*,*) '慣性モーメント: ', moment_inertia(1), ' kg⋅m²'
    end subroutine initialize_simulation

    !> 解析データ初期化サブルーチン
    subroutine initialize_analysis_data
        implicit none
        
        in_contact = .false.
        contact_start_step = 0
        contact_end_step = 0
        max_overlap = 0.0d0
        pre_collision_velocity = 0.0d0
        post_collision_velocity = 0.0d0
        
        initial_potential_energy = mass(1) * GRAVITY_ACCEL * z_coord(1)
        initial_kinetic_energy = 0.0d0
        
        falling_phase = .true.
        rising_phase = .false.
        collision_detected = .false.
        max_rebound_height = 0.0d0
        theoretical_rebound_height = restitution_coefficient**2 * drop_height
        
        write(*,*) '解析データ初期化完了'
        write(*,*) '初期位置エネルギー: ', initial_potential_energy, ' J'
        write(*,*) '理論反発高さ: ', theoretical_rebound_height, ' m'
    end subroutine initialize_analysis_data

    !> 状態をCSVファイルに保存するサブルーチン
    subroutine save_state_to_csv(step, time)
        implicit none
        integer, intent(in) :: step
        real(8), intent(in) :: time
        character(len=256) :: filename
        
        write(filename, '(A,I0,A)') 'data/state_step_', step, '.csv'
        
        open(unit=20, file=trim(filename), status='replace', action='write')
        
        ! ヘッダー書き込み
        write(20, '(A)') 'Step,Time,Particle,X,Z,Angle,Vx,Vz,Va,Fx,Fz,Fm'
        
        ! データ書き込み
        write(20, '(I0,A,F15.8,A,I0,A,F15.8,A,F15.8,A,F15.8,A,F15.8,A,F15.8,A,F15.8,A,F15.8,A,F15.8,A,F15.8)') &
            step, ',', time, ',', 1, ',', &
            x_coord(1), ',', z_coord(1), ',', rotation_angle(1), ',', &
            x_vel(1), ',', z_vel(1), ',', rotation_vel(1), ',', &
            x_force_sum(1), ',', z_force_sum(1), ',', moment_sum(1)
        
        close(20)
    end subroutine save_state_to_csv

    !> 自由落下解析サブルーチン
    subroutine analyze_free_fall(step, time)
        implicit none
        integer, intent(in) :: step
        real(8), intent(in) :: time
        real(8) :: particle_height, overlap_wall
        logical :: currently_in_contact
        
        particle_height = z_coord(1)
        currently_in_contact = (particle_height <= radius(1) + 1.0d-6)
        
        ! エネルギー計算
        current_potential_energy = mass(1) * GRAVITY_ACCEL * particle_height
        current_kinetic_energy = 0.5d0 * mass(1) * (x_vel(1)**2 + z_vel(1)**2) + &
                                0.5d0 * moment_inertia(1) * rotation_vel(1)**2
        
        ! 接触状態の変化を検出
        if (.not. in_contact .and. currently_in_contact) then
            ! 接触開始
            in_contact = .true.
            contact_start_step = step
            pre_collision_velocity = abs(z_vel(1))
            collision_detected = .true.
            
            write(*,*) '接触開始: ステップ ', step, ', 時刻 ', time
            write(*,*) '衝突前速度: ', pre_collision_velocity, ' m/s'
            
        else if (in_contact .and. .not. currently_in_contact) then
            ! 接触終了
            in_contact = .false.
            contact_end_step = step
            post_collision_velocity = abs(z_vel(1))
            falling_phase = .false.
            rising_phase = .true.
            
            write(*,*) '接触終了: ステップ ', step, ', 時刻 ', time
            write(*,*) '衝突後速度: ', post_collision_velocity, ' m/s'
            
            ! 接触統計の保存
            call save_contact_analysis(step, time)
        end if
        
        ! オーバーラップの計算と記録
        if (currently_in_contact) then
            overlap_wall = radius(1) - particle_height
            if (overlap_wall > max_overlap) then
                max_overlap = overlap_wall
            end if
        end if
        
        ! 最高点の検出
        if (rising_phase .and. z_vel(1) <= 0.0d0) then
            rising_phase = .false.
            falling_phase = .true.
            max_rebound_height = particle_height
            
            write(*,*) '最高点到達: ステップ ', step, ', 時刻 ', time
            write(*,*) '反発高さ: ', max_rebound_height, ' m'
        end if
        
        ! 詳細解析データの保存
        if (mod(step, output_interval) == 0) then
            call save_detailed_analysis(step, time)
        end if
    end subroutine analyze_free_fall

    !> 接触解析結果保存サブルーチン
    subroutine save_contact_analysis(step, time)
        implicit none
        integer, intent(in) :: step
        real(8), intent(in) :: time
        real(8) :: contact_duration, overlap_ratio, restitution_calc
        
        contact_duration = real(contact_end_step - contact_start_step, 8) * time_step
        overlap_ratio = (max_overlap / radius(1)) * 100.0d0
        
        if (abs(pre_collision_velocity) > 1.0d-10) then
            restitution_calc = post_collision_velocity / pre_collision_velocity
        else
            restitution_calc = 0.0d0
        end if
        
        open(unit=21, file='data/analysis/contact_analysis.csv', status='unknown', position='append', action='write')
        
        ! ファイルが空の場合はヘッダーを書き込み
        if (contact_start_step == step) then
            write(21, '(A)') 'Contact_Start_Step,Contact_End_Step,Duration_s,Pre_Velocity_m/s,Post_Velocity_m/s,' // &
                           'Max_Overlap_m,Overlap_Ratio_%,Calculated_Restitution,Set_Restitution'
        end if
        
        write(21, '(I0,A,I0,A,F12.6,A,F12.6,A,F12.6,A,F12.8,A,F8.4,A,F8.6,A,F8.6)') &
            contact_start_step, ',', contact_end_step, ',', contact_duration, ',', &
            pre_collision_velocity, ',', post_collision_velocity, ',', &
            max_overlap, ',', overlap_ratio, ',', restitution_calc, ',', restitution_coefficient
        
        close(21)
        
        ! オーバーラップをリセット
        max_overlap = 0.0d0
    end subroutine save_contact_analysis

    !> 詳細解析データ保存サブルーチン
    subroutine save_detailed_analysis(step, time)
        implicit none
        integer, intent(in) :: step
        real(8), intent(in) :: time
        real(8) :: total_energy, energy_error
        
        total_energy = current_potential_energy + current_kinetic_energy
        energy_error = abs(total_energy - (initial_potential_energy + initial_kinetic_energy))
        
        open(unit=22, file='data/detailed/energy_analysis.csv', status='unknown', position='append', action='write')
        
        ! ヘッダー（初回のみ）
        if (step == output_interval) then
            write(22, '(A)') 'Step,Time_s,Height_m,Velocity_m/s,Potential_Energy_J,' // &
                           'Kinetic_Energy_J,Total_Energy_J,Energy_Error_J'
        end if
        
        write(22, '(I0,A,F12.6,A,F12.6,A,F12.6,A,F15.8,A,F15.8,A,F15.8,A,F15.10)') &
            step, ',', time, ',', z_coord(1), ',', z_vel(1), ',', &
            current_potential_energy, ',', current_kinetic_energy, ',', total_energy, ',', energy_error
        
        close(22)
    end subroutine save_detailed_analysis

    !> 進行状況表示サブルーチン
    subroutine print_progress(step, max_steps, time)
        implicit none
        integer, intent(in) :: step, max_steps
        real(8), intent(in) :: time
        real(8) :: progress_percent
        
        progress_percent = (real(step, 8) / real(max_steps, 8)) * 100.0d0
        
        write(*, '(A,I0,A,I0,A,F5.1,A,A,F8.4,A,A,F8.4,A)') &
            'ステップ: ', step, '/', max_steps, ' (', progress_percent, '%), ', &
            '時刻: ', time, ' s, ', '高さ: ', z_coord(1), ' m'
    end subroutine print_progress

    !> 最終結果出力サブルーチン
    subroutine output_final_results(final_step, final_time, elapsed_time)
        implicit none
        integer, intent(in) :: final_step
        real(8), intent(in) :: final_time, elapsed_time
        real(8) :: avg_step_time
        
        avg_step_time = elapsed_time / real(final_step, 8) * 1000.0d0  ! ミリ秒
        
        open(unit=23, file='data/simulation_summary.csv', status='replace', action='write')
        
        write(23, '(A)') 'Parameter,Value'
        write(23, '(A,I0)') 'Total_Steps,', final_step
        write(23, '(A,F12.6)') 'Total_Time_s,', final_time
        write(23, '(A,F12.6)') 'Computation_Time_s,', elapsed_time
        write(23, '(A,F12.6)') 'Avg_Step_Time_ms,', avg_step_time
        write(23, '(A,F12.6)') 'Drop_Height_m,', drop_height
        write(23, '(A,F12.6)') 'Particle_Radius_m,', particle_radius
        write(23, '(A,F12.6)') 'Particle_Mass_kg,', mass(1)
        write(23, '(A,F12.6)') 'Restitution_Coefficient,', restitution_coefficient
        write(23, '(A,F12.6)') 'Max_Rebound_Height_m,', max_rebound_height
        write(23, '(A,F12.6)') 'Theoretical_Rebound_Height_m,', theoretical_rebound_height
        
        close(23)
        
        write(*,*) ''
        write(*,*) '========================================='
        write(*,*) '          シミュレーション結果'
        write(*,*) '========================================='
        write(*,*) '総ステップ数: ', final_step
        write(*,*) '総計算時間: ', elapsed_time, ' 秒'
        write(*,*) '1ステップ平均時間: ', avg_step_time, ' ミリ秒'
        if (collision_detected) then
            write(*,*) '最大反発高さ: ', max_rebound_height, ' m'
            write(*,*) '理論反発高さ: ', theoretical_rebound_height, ' m'
        end if
    end subroutine output_final_results

    !> 理論解との比較サブルーチン
    subroutine compare_with_theory
        implicit none
        real(8) :: height_error_percent, free_fall_time, theoretical_impact_velocity
        real(8) :: energy_loss_expected, energy_loss_actual, energy_error_percent
        
        if (.not. collision_detected) then
            write(*,*) '警告: 衝突が検出されませんでした。理論解との比較を実行できません。'
            return
        end if
        
        ! 自由落下時間の理論値
        free_fall_time = sqrt(2.0d0 * drop_height / GRAVITY_ACCEL)
        
        ! 衝突前速度の理論値
        theoretical_impact_velocity = sqrt(2.0d0 * GRAVITY_ACCEL * drop_height)
        
        ! 高さの相対誤差
        if (theoretical_rebound_height > 0.0d0) then
            height_error_percent = abs(max_rebound_height - theoretical_rebound_height) / &
                                  theoretical_rebound_height * 100.0d0
        else
            height_error_percent = 0.0d0
        end if
        
        ! エネルギー損失の比較
        energy_loss_expected = initial_potential_energy * (1.0d0 - restitution_coefficient**2)
        energy_loss_actual = initial_potential_energy - (mass(1) * GRAVITY_ACCEL * max_rebound_height)
        
        if (energy_loss_expected > 0.0d0) then
            energy_error_percent = abs(energy_loss_actual - energy_loss_expected) / &
                                  energy_loss_expected * 100.0d0
        else
            energy_error_percent = 0.0d0
        end if
        
        ! 理論解比較結果の保存
        open(unit=24, file='data/theory_comparison.csv', status='replace', action='write')
        
        write(24, '(A)') 'Parameter,Theoretical,Actual,Error_%'
        write(24, '(A,F12.6,A,F12.6,A,F8.4)') 'Free_Fall_Time_s,', free_fall_time, ',', &
            real(contact_start_step, 8) * time_step, ',', &
            abs(real(contact_start_step, 8) * time_step - free_fall_time) / free_fall_time * 100.0d0
        write(24, '(A,F12.6,A,F12.6,A,F8.4)') 'Impact_Velocity_m/s,', theoretical_impact_velocity, ',', &
            pre_collision_velocity, ',', &
            abs(pre_collision_velocity - theoretical_impact_velocity) / theoretical_impact_velocity * 100.0d0
        write(24, '(A,F12.6,A,F12.6,A,F8.4)') 'Rebound_Height_m,', theoretical_rebound_height, ',', &
            max_rebound_height, ',', height_error_percent
        write(24, '(A,F12.6,A,F12.6,A,F8.4)') 'Energy_Loss_J,', energy_loss_expected, ',', &
            energy_loss_actual, ',', energy_error_percent
        
        close(24)
        
        write(*,*) ''
        write(*,*) '========================================='
        write(*,*) '          理論解との比較結果'
        write(*,*) '========================================='
        write(*,*) '自由落下時間:'
        write(*,*) '  理論値: ', free_fall_time, ' s'
        write(*,*) '  計算値: ', real(contact_start_step, 8) * time_step, ' s'
        write(*,*) '衝突前速度:'
        write(*,*) '  理論値: ', theoretical_impact_velocity, ' m/s'
        write(*,*) '  計算値: ', pre_collision_velocity, ' m/s'
        write(*,*) '反発高さ:'
        write(*,*) '  理論値: ', theoretical_rebound_height, ' m'
        write(*,*) '  計算値: ', max_rebound_height, ' m'
        write(*,*) '  相対誤差: ', height_error_percent, ' %'
        write(*,*) 'エネルギー損失:'
        write(*,*) '  理論値: ', energy_loss_expected, ' J'
        write(*,*) '  計算値: ', energy_loss_actual, ' J'
        write(*,*) '  相対誤差: ', energy_error_percent, ' %'
    end subroutine compare_with_theory

    ! 以下、元のpem_simulator.f90から必要なサブルーチンを簡略化して実装
    ! (セル更新、接触判定、位置更新などの核となる計算部分)
    
    subroutine ncel_sub
        ! セル更新の簡略実装（1粒子なので簡単な処理）
        implicit none
        particle_cell_idx(1) = 1
    end subroutine ncel_sub
    
    subroutine wcont_sub(particle_idx)
        ! 壁との接触判定・力計算
        implicit none
        integer, intent(in) :: particle_idx
        real(8) :: overlap, normal_force, contact_stiffness
        
        ! 床との接触判定
        if (z_coord(particle_idx) <= radius(particle_idx)) then
            overlap = radius(particle_idx) - z_coord(particle_idx)
            
            ! 簡単なばね-ダッシュポット模型
            contact_stiffness = young_modulus_wall * radius(particle_idx)
            normal_force = contact_stiffness * overlap
            
            ! 反発係数を考慮した減衰力
            normal_force = normal_force - 2.0d0 * sqrt(contact_stiffness * mass(particle_idx)) * &
                          log(restitution_coefficient) / sqrt(PI_VAL**2 + log(restitution_coefficient)**2) * z_vel(particle_idx)
            
            if (normal_force > 0.0d0) then
                z_force_sum(particle_idx) = z_force_sum(particle_idx) + normal_force
            end if
        end if
    end subroutine wcont_sub
    
    subroutine pcont_sub(particle_idx, rmax_radius)
        ! 粒子間接触判定（1粒子なので実質何もしない）
        implicit none
        integer, intent(in) :: particle_idx
        real(8), intent(in) :: rmax_radius
        ! 1粒子のシミュレーションなので実装不要
    end subroutine pcont_sub
    
    subroutine nposit_sub(static_flag)
        ! 位置・速度更新（オイラー法）
        implicit none
        integer, intent(out) :: static_flag
        real(8) :: acceleration_x, acceleration_z, angular_acceleration
        integer :: i
        
        static_flag = 0
        
        do i = 1, num_particles
            ! 加速度計算
            acceleration_x = x_force_sum(i) / mass(i)
            acceleration_z = z_force_sum(i) / mass(i)
            angular_acceleration = moment_sum(i) / moment_inertia(i)
            
            ! 速度更新
            x_vel(i) = x_vel(i) + acceleration_x * time_step
            z_vel(i) = z_vel(i) + acceleration_z * time_step
            rotation_vel(i) = rotation_vel(i) + angular_acceleration * time_step
            
            ! 変位増分計算
            x_disp_incr(i) = x_vel(i) * time_step
            z_disp_incr(i) = z_vel(i) * time_step
            rot_disp_incr(i) = rotation_vel(i) * time_step
            
            ! 位置更新
            x_coord(i) = x_coord(i) + x_disp_incr(i)
            z_coord(i) = z_coord(i) + z_disp_incr(i)
            rotation_angle(i) = rotation_angle(i) + rot_disp_incr(i)
            
            ! 静止判定（速度が十分小さく、位置が低い場合）
            if (abs(x_vel(i)) < 1.0d-8 .and. abs(z_vel(i)) < 1.0d-8 .and. z_coord(i) < drop_height * 0.01d0) then
                static_flag = 1
            end if
        end do
    end subroutine nposit_sub

end program free_fall_simulation 