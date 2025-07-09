! Optimized version of the 2D Particle Element Method simulator
! Performance improvements include:
! - OpenMP parallelization for particle loops
! - Vectorization-friendly code structure
! - Improved memory access patterns
! - Optimized mathematical operations

! Module: Simulation constants (array sizes, mathematical constants)
module simulation_constants_mod
    implicit none
    integer, parameter :: ni_max = 1000  ! Maximum particles
    integer, parameter :: nj_max = 13    ! Maximum contact points per particle
    integer, parameter :: nc_max = 20000 ! Maximum cells in grid
    real(8), parameter :: PI_VAL = 3.141592653589793d0
    real(8), parameter :: GRAVITY_ACCEL = 9.80665d0
    
    ! Cache-friendly constants
    integer, parameter :: CACHE_LINE_SIZE = 64  ! bytes
    integer, parameter :: REAL8_PER_CACHE_LINE = CACHE_LINE_SIZE / 8
end module simulation_constants_mod

! Module: Simulation control parameters and material properties
module simulation_parameters_mod
    use simulation_constants_mod, only: ni_max
    implicit none

    ! Simulation control parameters
    real(8) :: time_step                      
    real(8) :: friction_coeff_particle        
    real(8) :: friction_coeff_wall            
    real(8) :: young_modulus_particle         
    real(8) :: young_modulus_wall             
    real(8) :: poisson_ratio_particle         
    real(8) :: poisson_ratio_wall             
    real(8) :: shear_to_normal_stiffness_ratio
    real(8) :: particle_density               
    
    ! Particle generation parameters
    real(8) :: particle_radius_large          
    real(8) :: particle_radius_small          
    real(8) :: container_width                
    integer :: particle_gen_layers            
    integer :: random_seed                    
    
    ! Cell algorithm control parameters
    logical :: disable_cell_algorithm         
    real(8) :: cell_size_override            
    
    ! Validation mode settings
    logical :: validation_mode                
    real(8) :: validation_particle1_x        
    real(8) :: validation_particle1_z        
    real(8) :: validation_particle2_x        
    real(8) :: validation_particle2_z        
    real(8) :: validation_particle1_vx       
    real(8) :: validation_particle2_vx       
    real(8) :: validation_particle_radius    
    
    ! Output control parameters
    integer :: output_interval_normal         
    integer :: output_interval_validation     
    integer :: max_calculation_steps          
    
    ! Performance parameters
    integer :: num_threads                    ! Number of OpenMP threads
    logical :: enable_vectorization           ! Enable vectorization hints

    save
end module simulation_parameters_mod

! Module: Particle data with optimized memory layout
module particle_data_mod
    use simulation_constants_mod, only: ni_max, nj_max
    implicit none

    ! Physical properties - aligned for vectorization
    real(8), dimension(ni_max), align(64) :: radius         
    real(8), dimension(ni_max), align(64) :: mass           
    real(8), dimension(ni_max), align(64) :: moment_inertia 

    ! Position and orientation - aligned for vectorization
    real(8), dimension(ni_max), align(64) :: x_coord        
    real(8), dimension(ni_max), align(64) :: z_coord        
    real(8), dimension(ni_max), align(64) :: rotation_angle 

    ! Velocities - aligned for vectorization
    real(8), dimension(ni_max), align(64) :: x_vel          
    real(8), dimension(ni_max), align(64) :: z_vel          
    real(8), dimension(ni_max), align(64) :: rotation_vel   

    ! Forces and moments - aligned for vectorization
    real(8), dimension(ni_max), align(64) :: x_force_sum    
    real(8), dimension(ni_max), align(64) :: z_force_sum    
    real(8), dimension(ni_max), align(64) :: moment_sum     

    ! Contact forces - reorganized for better cache usage
    real(8), dimension(ni_max, nj_max) :: normal_force_contact  
    real(8), dimension(ni_max, nj_max) :: shear_force_contact   
    integer, dimension(ni_max, nj_max) :: contact_partner_idx 

    ! Displacement increments - aligned for vectorization
    real(8), dimension(ni_max), align(64) :: x_disp_incr    
    real(8), dimension(ni_max), align(64) :: z_disp_incr    
    real(8), dimension(ni_max), align(64) :: rot_disp_incr  
    
    save
end module particle_data_mod

! Module: Cell grid system data
module cell_system_mod
    use simulation_constants_mod, only: ni_max, nc_max
    implicit none

    integer :: num_particles          
    integer :: cells_x_dir            
    integer :: cells_z_dir            

    real(8) :: cell_size              

    ! Cell data structures
    integer, dimension(nc_max) :: cell_head        
    integer, dimension(ni_max) :: particle_cell_next 
    integer, dimension(nc_max) :: cell_particle_map
    integer, dimension(ni_max) :: particle_cell_idx 

    save
end module cell_system_mod

! Main program with OpenMP support
program optimized_two_dimensional_pem
    use omp_lib
    use simulation_constants_mod
    use simulation_parameters_mod
    use particle_data_mod
    use cell_system_mod
    implicit none

    integer :: it_step, static_judge_flag          
    integer :: i 
    real(8) :: current_time                        
    real(8) :: rmax_particle_radius                
    
    ! Performance timing variables
    integer :: start_time, end_time, clock_rate
    real(8) :: elapsed_time, wall_time_start, wall_time_end
    
    ! Validation mode variables
    logical :: collision_started = .false.  
    logical :: collision_finished = .false. 
    real(8) :: initial_v1, initial_v2        
    real(8) :: final_v1, final_v2            
    real(8) :: dist, sumr                    
    real(8) :: m1, m2, e_coeff, v1_theo, v2_theo, err1, err2 
    real(8) :: ke_initial, ke_final, ke_theo 

    ! Initialize OpenMP
    call omp_set_dynamic(.false.)
    
    ! Start timing
    call system_clock(start_time, clock_rate)
    wall_time_start = omp_get_wtime()
    
    ! Read input parameters
    call read_input_file
    
    ! Set number of threads
    if (num_threads > 0) then
        call omp_set_num_threads(num_threads)
    else
        num_threads = omp_get_max_threads()
    end if
    
    write(*,*) '================================='
    write(*,*) 'Optimized PEM Simulator'
    write(*,*) 'OpenMP threads: ', num_threads
    write(*,*) '================================='
    
    ! Initialize simulation
    call fposit_sub(rmax_particle_radius)
    call inmat_sub
    call init_sub

    ! Validation mode setup
    if (validation_mode) then
        x_vel(1) = validation_particle1_vx
        x_vel(2) = validation_particle2_vx
        z_vel(1) = 0.0d0
        z_vel(2) = 0.0d0
        initial_v1 = x_vel(1)
        initial_v2 = x_vel(2)
        friction_coeff_particle = 0.0d0
        friction_coeff_wall = 0.0d0
    end if

    current_time = 0.0d0

    ! Main simulation loop
    do it_step = 1, max_calculation_steps
        current_time = current_time + time_step

        ! Update cell structure
        call ncel_sub

        ! Clear forces - vectorized
        !$omp parallel do simd
        do i = 1, num_particles
            x_force_sum(i) = 0.0d0
            z_force_sum(i) = 0.0d0
            moment_sum(i) = 0.0d0
        end do
        !$omp end parallel do simd
        
        ! Calculate contact forces with OpenMP parallelization
        !$omp parallel do schedule(dynamic,4)
        do i = 1, num_particles
            call wcont_sub(i)
            call pcont_sub_optimized(i, rmax_particle_radius)
        end do
        !$omp end parallel do

        ! Update positions and velocities
        call nposit_sub_optimized(static_judge_flag)

        ! Validation mode collision detection
        if (validation_mode .and. .not. collision_finished) then
            call validate_collision(collision_started, collision_finished, &
                                   initial_v1, initial_v2, final_v1, final_v2, &
                                   current_time)
            if (collision_finished) goto 200
        end if

        ! Check for static state
        if (static_judge_flag == 1) then
            write(*,*) 'Static state reached at time: ', current_time
            goto 200
        end if

        ! Progress output
        if (mod(it_step, 1000) == 0) then
            write(*, '(A,F10.6,A,F12.6,A,F12.6)') 'Time= ', current_time, &
                                                 ' Z0(N)= ', z_coord(num_particles), &
                                                 ' V0(N)= ', z_vel(num_particles)
        end if

        ! Data output
        if (validation_mode) then
            if (it_step == 1 .or. mod(it_step, output_interval_validation) == 0) then
                call gfout_sub(it_step, current_time, rmax_particle_radius)
            end if
        else
            if (it_step == 1 .or. mod(it_step, output_interval_normal) == 0) then
                call gfout_sub(it_step, current_time, rmax_particle_radius)
            end if
        end if
    end do

200 continue

    ! Output backup data
    call bfout_sub

    ! Close files
    close(10) 
    close(11) 
    close(13) 

    ! End timing
    call system_clock(end_time)
    elapsed_time = real(end_time - start_time) / real(clock_rate)
    wall_time_end = omp_get_wtime()
    
    ! Performance report
    write(*,*) '================================='
    write(*,*) 'Simulation Performance Report'
    write(*,*) '================================='
    write(*,*) 'Particles: ', num_particles
    write(*,*) 'Steps completed: ', it_step
    write(*,*) 'CPU time: ', elapsed_time, ' seconds'
    write(*,*) 'Wall time: ', wall_time_end - wall_time_start, ' seconds'
    write(*,*) 'Average time per step: ', elapsed_time / real(it_step), ' seconds'
    write(*,*) 'OpenMP speedup: ', elapsed_time / (wall_time_end - wall_time_start)
    
    if (disable_cell_algorithm .or. cell_size_override > 0.0d0) then
        write(*,*) 'Cell algorithm: DISABLED'
    else
        write(*,*) 'Cell algorithm: ENABLED'
        write(*,*) 'Cell size: ', cell_size
        write(*,*) 'Grid: ', cells_x_dir, ' x ', cells_z_dir
        write(*,*) 'Total cells: ', cells_x_dir * cells_z_dir
    end if
    
    write(*,*) '================================='

    stop

contains

    ! Include all subroutines from original code with optimizations
    include 'pem_subroutines_optimized.inc'

end program optimized_two_dimensional_pem