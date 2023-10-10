program plane_wave
    implicit none
    real, parameter :: pi = 3.1415926535                                                    ! pi
    real, parameter :: c = 3.0e8                                                            ! 光速
    real, parameter :: f = 1.0e7                                                            ! 频率

    real :: omega = 2*pi*f                                                                  ! 圆频率
    real :: lambda = c/f                                                                    ! 波长
    real :: dlt, dt                                                                         ! dlt: 空间步长, dt: 时间步长
    integer :: n, i, j, k                                                                   ! n: 时间迭代变量, i,j,k: 空间循环变量
    integer, parameter :: t_max = 600                                                       ! t_max: 时间范围指标
    integer, parameter :: x_min = -4*15, x_max = 4*15, &                                    ! 空间范围指标
        y_min = -4*15, y_max = 4*15, z_min = -4*15, z_max = 4*15
    integer, parameter :: xc_min = -2*15, xc_max = 2*15, &                                  ! 总场区范围
        yc_min = -2*15, yc_max = 2*15

    real, dimension(x_min:x_max, y_min:y_max) :: Ezn = 0.                                   ! n+1 步 Ez 数组
    real, dimension(x_min:x_max, y_min:y_max) :: Ezp = 0.                                   ! n 步 Ez 数组
    real, dimension(x_min:x_max, y_min:y_max-1) :: Hx = 0.                                  ! Hx 数组
    real, dimension(x_min:x_max-1, y_min:y_max) :: Hy = 0.                                  ! Hy 数组
    real, dimension(x_min:x_max, y_min:y_max) :: Ezn_in = 0.                                ! 入射场 n+1 步 Ez 数组
    real, dimension(x_min:x_max, y_min:y_max) :: Ezp_in = 0.                                ! 入射场 n 步 Ez 数组
    real, dimension(x_min:x_max, y_min:y_max-1) :: Hx_in = 0.                               ! 入射场 Hx 数组
    real, dimension(x_min:x_max-1, y_min:y_max) :: Hy_in = 0.                               ! 入射场 Hy 数组

    real, parameter :: epsilon = 8.8541878e-12, mu = pi*4e-7                                ! 真空介电系数, 真空磁导系数
    real :: sigma = 0., sigma_m = 0.                                                        ! 电导率, 磁导率
    real :: CA, CB, CP, CQ                                                                  ! FDTD 迭代常数
    real :: C1, C2, C3                                                                      ! 吸收边界常数

    character(len=20) :: filename                                                           ! 输出文件
    integer :: status                                                                       ! 文件打开状态
    character(len=50) :: msg                                                                ! 文件打开信息


    dlt = lambda/15.                                                                        ! 空间步长
    dt = dlt/(2*c)                                                                          ! Courant 稳定性条件

    CA = (epsilon/dt-sigma/2.)/(epsilon/dt+sigma/2.)                                        ! CA
    CB = 1./(epsilon/dt+sigma/2.)                                                           ! CB
    CP = (mu/dt-sigma_m/2.)/(mu/dt+sigma_m/2.)                                              ! CP
    CQ = 1./(mu/dt+sigma_m/2.)                                                              ! CQ

    C1 = (c*dt-dlt)/(c*dt+dlt)
    C2 = (c**2*mu*dt)/(2*(c*dt+dlt))
    C3 = (c*dt-sqrt(2.)*dlt)/(c*dt+sqrt(2.)*dlt)

    filename = "Ez.txt"
    open(unit = 20, file = filename, status = "replace", action = "write", iostat=status, iomsg = msg)
    if ( status /= 0 ) stop "Error opening file name"


    do n = 1, t_max
        ! 平面入射波
        do j = y_min, y_max
            Hy_in(x_min, j) = sin(omega*n*dt)
            do i = x_min+1, x_max-1
                Hy_in(i, j) = CP*Hy_in(i, j)-CQ/dlt*(Ezp_in(i+1, j)-Ezp_in(i, j))                    
            end do
        end do
        do j = y_min, y_max
            do i = x_min+1, x_max-1
                Ezn_in(i, j) = CA*Ezp_in(i, j)-CB/dlt*(Hy_in(i, j)-Hy_in(i-1, j))
            end do
            Ezn_in(x_max, j) = Ezp_in(x_max-1, j)+C1*(Ezn_in(x_max-1, j)-Ezp_in(x_max, j))
            Ezn_in(x_min, j) = Ezp_in(x_min+1, j)+C1*(Ezn_in(x_min+1, j)-Ezp_in(x_min, j))
        end do

        ! FDTD 迭代计算
        forall (i = x_min:x_max, j = y_min:y_max-1)                                         ! Hx 迭代计算
            Hx(i, j) = CP*Hx(i, j)-CQ/dlt*(Ezp(i, j+1)-Ezp(i, j))                           
        end forall
        do i = xc_min, xc_max                                                               ! Hx 总场边界处理
            Hx(i, yc_min-1) = Hx(i, yc_min-1)-dt/mu/dlt*Ezp_in(i, yc_min)                                      
            Hx(i, yc_max) = Hx(i, yc_max)+dt/mu/dlt*Ezp_in(i, yc_max)
        end do

        forall (i = x_min:x_max-1, j = y_min:y_max)                                         ! Hy 迭代计算
            Hy(i, j) = CP*Hy(i, j)+CQ/dlt*(Ezp(i+1, j)-Ezp(i, j))
        end forall
        do j = yc_min, yc_max                                                               ! Hy 总场边界处理
            Hy(xc_min-1, j) = Hy(xc_min-1, j)+dt/mu/dlt*Ezp_in(xc_min, j)
            Hy(xc_max, j) = Hy(xc_max, j)-dt/mu/dlt*Ezp_in(xc_max, j)
        end do

        forall (i = x_min+1:x_max-1, j = y_min+1:y_max-1)                                   ! Ez 迭代计算
            Ezn(i, j) = CA*Ezp(i, j)+CB/dlt*((Hy(i, j)-Hy(i-1, j))-(Hx(i, j)-Hx(i, j-1)))
        end forall
        do i = xc_min+1, xc_max-1                                                           ! Ez 总场边界处理
            Ezn(i, yc_min) = Ezn(i, yc_min)+dt/epsilon/dlt*Hx_in(i, yc_min-1)
            Ezn(i, yc_max) = Ezn(i, yc_max)-dt/epsilon/dlt*Hx_in(i, yc_max)
        end do
        do j = yc_min+1, yc_max-1                                                           ! Ez 总场边界处理
            Ezn(xc_min, j) = Ezn(xc_min, j)-dt/epsilon/dlt*Hy_in(xc_min-1, j)
            Ezn(xc_max, j) = Ezn(xc_max, j)+dt/epsilon/dlt*Hy_in(xc_max, j)            
        end do
        ! 总场区角点处理
        Ezn(xc_min, yc_min) = Ezn(xc_min, yc_min)-dt/epsilon/dlt*(Hy_in(xc_min-1, yc_min)-Hx_in(xc_min, yc_min-1))
        Ezn(xc_min, yc_max) = Ezn(xc_min, yc_max)-dt/epsilon/dlt*(Hy_in(xc_min-1, yc_max)+Hx_in(xc_min, yc_max))
        Ezn(xc_max, yc_min) = Ezn(xc_max, yc_min)+dt/epsilon/dlt*(Hy_in(xc_max, yc_min)+Hx_in(xc_max, yc_min-1))
        Ezn(xc_max, yc_max) = Ezn(xc_max, yc_max)+dt/epsilon/dlt*(Hy_in(xc_max, yc_max)-Hx_in(xc_max, yc_max))


        ! 二维边界: 二阶 Mur 吸收边界
        do j = y_min+1, y_max-1
            Ezn(x_min, j) = Ezp(x_min+1, j)+C1*(Ezn(x_min+1, j)-Ezp(x_min, j)) &
                & -C2*(Hx(x_min, j)+Hx(x_min+1, j)-Hx(x_min, j-1)-Hx(x_min+1, j-1))
            Ezn(x_max, j) = Ezp(x_max-1, j)+C1*(Ezn(x_max-1, j)-Ezp(x_max, j)) &
                & -C2*(Hx(x_max, j)+Hx(x_max-1, j)-Hx(x_max, j-1)-Hx(x_max-1, j-1))
        end do
        do i = x_min+1, x_max-1
            Ezn(i, y_min) = Ezp(i, y_min+1)+C1*(Ezn(i, y_min+1)-Ezp(i, y_min)) &
                & +C2*(Hy(i, y_min)+Hy(i, y_min+1)-Hy(i-1, y_min)-Hy(i-1, y_min+1))
            Ezn(i, y_max) = Ezp(i, y_max-1)+C1*(Ezn(i, y_max-1)-Ezp(i, y_max)) &
                & +C2*(Hy(i, y_max)+Hy(i, y_max-1)-Hy(i-1, y_max)-Hy(i-1, y_max-1))
        end do

        ! 角点处理: 一阶 Mur 吸收边界
        Ezn(x_min, y_min) = Ezp(x_min+1, y_min+1)+C3*(Ezn(x_min+1, y_min+1)-Ezp(x_min, y_min))
        Ezn(x_max, y_min) = Ezp(x_max-1, y_min+1)+C3*(Ezn(x_max-1, y_min+1)-Ezp(x_max, y_min))
        Ezn(x_min, y_max) = Ezp(x_min+1, y_max-1)+C3*(Ezn(x_min+1, y_max-1)-Ezp(x_min, y_max))
        Ezn(x_max, y_max) = Ezp(x_max-1, y_max-1)+C3*(Ezn(x_max-1, y_max-1)-Ezp(x_max, y_max))
        
        Ezp_in = Ezn_in
        Ezp = Ezn
        ! 计算结果输出
        do j = y_min, y_max
            write(20, '(121F15.5)') (Ezp(i, j), i = x_min, x_max)
        end do
    end do

    ! 计算结果输出
    ! do j = y_min, y_max
    !     write(20, '(121F15.5)') (Ezn(i, j), i = x_min, x_max)
    ! end do

    ! filename = "coordinate.txt"
    ! open(unit = 20, file = filename, status = "replace", action = "write", iostat=status, iomsg = msg)
    ! if ( status /= 0 ) stop "Error opening file name"
    ! write(20, '(121F15.5)') (i*dlt, i = x_min, x_max)
    ! write(20, '(121F15.5)') (j*dlt, j = y_min, y_max)


end program plane_wave