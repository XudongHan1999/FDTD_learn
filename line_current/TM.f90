program TM

    implicit none
    real, parameter :: pi = 3.1415926535                                ! pi
    real, parameter :: c = 3.0e8                                        ! 光速
    real, parameter :: f = 1.0e7                                        ! 频率
    real :: omega = 2*pi*f                                              ! 圆频率
    real :: lambda, T, delta, delta_t                                   ! 波长, 周期

    real :: epsilon = 8.8541878e-12, mu = pi*4e-7, eta                  ! 介电系数, 磁导系数
    real :: sigma = 0., sigma_m = 0.                                    ! 电导率, 磁导率
    real :: CA, CB, CP, CQ                                              ! CA, CB, CP, CQ

    real, dimension(-5*20:5*20, -5*20:5*20, 0:400) :: E_z = 0           ! E_z 数组
    real, dimension(-5*20:5*20, -5*20:5*20-1, 0:400) :: H_x = 0         ! H_x 数组  
    real, dimension(-5*20:5*20-1, -5*20:5*20, 0:400) :: H_y = 0         ! H_y 数组
    real, dimension(-5*20:5*20, -5*20:5*20, 0:400) :: j_l = 0           ! 电流数组
    real, dimension(-5*20:5*20) :: x, y                                 ! (x, y) 坐标数组
    integer :: n, i, j, bound_l, bound_r, bound_b, bound_t              ! n:时间循环变量, i,j:空间循环变量

    character(len=20) :: filename                                       ! 输出文件
    integer :: status                                                   ! 文件打开状态
    character(len=50) :: msg                                            ! 文件打开信息

    lambda = c/f                                                        ! 波长
    T = 2.*pi/omega                                                     ! 周期
    delta = lambda/20.                                                  ! 空间步长
    delta_t = T/30.                                                     ! 时间步长
    bound_l = lbound(E_z, dim=1)
    bound_r = ubound(E_z, dim=1)
    bound_b = lbound(E_z, dim=2)
    bound_t = ubound(E_z, dim=2)
    forall (i = bound_l:bound_r)
        x(i) = i*delta                                                  ! x 坐标
    end forall
    forall (j = bound_b:bound_t)
        y(j) = j*delta                                                  ! y 坐标        
    end forall

    CA = (epsilon/delta_t-sigma/2.)/(epsilon/delta_t+sigma/2.)
    CB = 1./(epsilon/delta_t+sigma/2.)
    CP = (mu/delta_t-sigma_m/2.)/(mu/delta_t+sigma_m/2.)
    CQ = 1./(mu/delta_t+sigma_m/2.)

    ! 线电流源
    do n = 0, 400
        j_l(0, 0, n) = sin(omega*n*delta_t)
    end do

    do n = 1, 400
        ! FDTD 计算
        forall (i = bound_l:bound_r, j = bound_b:bound_t-1)
            H_x(i, j, n) = CP*H_x(i, j, n-1)-CQ/delta*(E_z(i, j+1, n-1)-E_z(i, j, n-1))          
        end forall
        forall (i = bound_l:bound_r-1, j = bound_b:bound_t)
            H_y(i, j, n) = CP*H_y(i, j, n-1)+CQ/delta*(E_z(i+1, j, n-1)-E_z(i, j, n-1))            
        end forall
        forall (i = bound_l+1:bound_r-1, j = bound_b+1:bound_t-1)
            E_z(i, j, n) = CA*E_z(i, j, n-1)+CB*((H_y(i, j, n)-H_y(i-1, j, n))/delta-(H_x(i, j, n)-H_x(i, j-1, n))/delta) &
                & -CB*j_l(i, j, n)/(delta**2)
        end forall

        ! 二维边界: 二阶 Mur 吸收边界
        do j = bound_b+1, bound_t-1
            E_z(bound_l, j, n) = E_z(bound_l+1, j, n-1)+(c*delta_t-delta)/(c*delta_t+delta) &
                & *(E_z(bound_l+1, j, n)-E_z(bound_l, j, n-1))-(c**2*mu*delta_t)/(2*(c*delta_t+delta)) &
                & *(H_x(bound_l, j, n)+H_x(bound_l+1, j, n)-H_x(bound_l, j-1, n)-H_x(bound_l+1, j-1, n))
            E_z(bound_r, j, n) = E_z(bound_r-1, j, n-1)+(c*delta_t-delta)/(c*delta_t+delta) &
                & *(E_z(bound_r-1, j, n)-E_z(bound_r, j, n-1))-(c**2*mu*delta_t)/(2*(c*delta_t+delta)) &
                & *(H_x(bound_r, j, n)+H_x(bound_r-1, j, n)-H_x(bound_r, j-1, n)-H_x(bound_r-1, j-1, n))           
        end do
        do i = bound_l+1, bound_r-1
            E_z(i, bound_b, n) = E_z(i, bound_b+1, n-1)+(c*delta_t-delta)/(c*delta_t+delta) &
                & *(E_z(i, bound_b+1, n)-E_z(i, bound_b, n-1))+(c**2*mu*delta_t)/(2*(c*delta_t+delta)) &
                & *(H_y(i, bound_b, n)+H_y(i, bound_b+1, n)-H_y(i-1, bound_b, n)-H_y(i-1, bound_b+1, n))
            E_z(i, bound_t, n) = E_z(i, bound_t-1, n-1)+(c*delta_t-delta)/(c*delta_t+delta) &
                & *(E_z(i, bound_t-1, n)-E_z(i, bound_t, n-1))+(c**2*mu*delta_t)/(2*(c*delta_t+delta)) &
                & *(H_y(i, bound_t, n)+H_y(i, bound_t-1, n)-H_y(i-1, bound_t, n)-H_y(i-1, bound_t-1, n))
        end do

        ! 角点处理: 一阶 Mur 吸收边界
        E_z(bound_l, bound_b, n) = E_z(bound_l+1, bound_b+1, n-1)+(c*delta_t-sqrt(2.)*delta)/(c*delta_t+sqrt(2.)*delta) &
            & *(E_z(bound_l+1, bound_b+1, n)-E_z(bound_l, bound_b, n-1))
        E_z(bound_r, bound_b, n) = E_z(bound_r-1, bound_b+1, n-1)+(c*delta_t-sqrt(2.)*delta)/(c*delta_t+sqrt(2.)*delta) &
            & *(E_z(bound_r-1, bound_b+1, n)-E_z(bound_r, bound_b, n-1))
        E_z(bound_l, bound_t, n) = E_z(bound_l+1, bound_t-1, n-1)+(c*delta_t-sqrt(2.)*delta)/(c*delta_t+sqrt(2.)*delta) &
            & *(E_z(bound_l+1, bound_t-1, n)-E_z(bound_l, bound_t, n-1))
        E_z(bound_r, bound_t, n) = E_z(bound_r-1, bound_t-1, n-1)+(c*delta_t-sqrt(2.)*delta)/(c*delta_t+sqrt(2.)*delta) &
            & *(E_z(bound_r-1, bound_t-1, n)-E_z(bound_r, bound_t, n-1))
    end do

    ! 结果输出
    filename = "coordinate.txt"
    open(unit = 20, file = filename, status = "replace", action = "write", iostat=status, iomsg = msg)
    if ( status /= 0 ) stop "Error opening file name"
    write(20, '(201F15.5)') (x(i), i = bound_l, bound_r)
    write(20, '(201F15.5)') (y(j), j = bound_b, bound_t)
    
    filename = "Ez.txt"
    open(unit = 20, file = filename, status = "replace", action = "write", iostat=status, iomsg = msg)
    if ( status /= 0 ) stop "Error opening file name"
    do n = 0, 400
        do i = bound_b, bound_t
            write(20, '(201F15.5)') (E_z(i, j, n), j = bound_l, bound_r)
        end do
    end do






end program TM
