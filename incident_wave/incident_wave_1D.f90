program plane_wave

    implicit none
    real, parameter :: pi = 3.1415926535                                                    ! pi
    real, parameter :: c = 3.0e8                                                            ! 光速
    real, parameter :: f = 1.0e7                                                            ! 频率

    real :: omega = 2*pi*f                                                                  ! 圆频率
    real :: lambda = c/f                                                                    ! 波长
    real :: dlt, dt                                                                         ! dlt: 空间步长, dt: 时间步长

    integer :: n, k                                                                         ! n: 时间迭代变量, k: 空间循环变量
    integer,parameter :: t_max = 1000                                                       ! t_max: 时间范围指标   
    integer,parameter :: Z_min = -5*20, Z_max = 5*20                                        ! 吸收边界指标
    integer,parameter :: Zt_min = -3*20, Zt_max = 3*20                                      ! 总场边界指标
    real, dimension(Z_min:Z_max) :: Ex = 0., Ex_n = 0.                                      ! Ex 数组
    real, dimension(Z_min:Z_max-1) :: Hy = 0.                                               ! Hy 数组
    real, dimension(Z_min:Z_max-1) :: Hy_in = 0.                                            ! 入射波 Hy
    real, dimension(Z_min:Z_max) :: Ex_in = 0., Exn_in = 0.                                 ! 入射波 Ex

    real :: epsilon = 8.8541878e-12, mu = pi*4e-7, eta                                      ! 介电系数, 磁导系数, 
    real :: sigma = 0., sigma_m = 0.                                                        ! 电导率, 磁导率
    real :: CA, CB, CP, CQ                                                                  ! FDTD 迭代常数
    real :: C1                                                                              ! 吸收边界计算常数

    character(len=20) :: filename                                                           ! 输出文件
    integer :: status                                                                       ! 文件打开状态
    character(len=50) :: msg                                                                ! 文件打开信息

    dlt = lambda/30.                                                                        ! Yee 元胞大小
    dt = dlt/(2*c)                                                                          ! 时间步长(Courant 稳定性条件)

    eta = sqrt(mu/epsilon)                                                                  ! 波阻抗

    CA = (epsilon/dt-sigma/2.)/(epsilon/dt+sigma/2.)                                        ! CA
    CB = 1./(epsilon/dt+sigma/2.)                                                           ! CB
    CP = (mu/dt-sigma_m/2.)/(mu/dt+sigma_m/2.)                                              ! CP
    CQ = 1./(mu/dt+sigma_m/2.)                                                              ! CQ

    C1 = (c*dt-dlt)/(c*dt+dlt)                                                              ! C1

    filename = "Hy.txt"
    open(unit = 20, file = filename, status = "replace", action = "write", iostat=status, iomsg = msg)
    if ( status /= 0 ) stop "Error opening file name"


    do n = 1, t_max
        ! 入射波
        do k = Z_min, Z_max
            if ( k == Z_min ) then
                Hy_in(k) = sin(omega*n*dt)
            else
                Hy_in(k) = CP*Hy_in(k)-CQ/dlt*(Ex_in(k+1)-Ex_in(k))
            end if
        end do
        do k = Z_min, Z_max
            Exn_in(k) = CA*Ex_in(k)-CB/dlt*(Hy_in(k)-Hy_in(k-1))
        end do
        Exn_in(Z_max) = Ex_in(Z_max-1)+C1*(Exn_in(Z_max-1)-Ex_in(Z_max))
        Exn_in(Z_min) = Ex_in(Z_min+1)+C1*(Exn_in(Z_min+1)-Ex_in(Z_min))


        ! FDTD 迭代计算
        do k = Z_min, Z_max-1
            Hy(k) = CP*Hy(k)-CQ/dlt*(Ex(k+1)-Ex(k))                                         ! Hy 迭代计算
        end do
        Hy(Zt_min-1) = Hy(Zt_min-1)+CQ*Ex_in(Zt_min)                                        ! Hy 左端散射场-总场边界
        Hy(Zt_max) = Hy(Zt_max)-dt/(mu*dlt)*Ex_in(Zt_max)                                   ! Hy 右端散射场-总场边界

        do k = Z_min+1, Z_max-1
            Ex_n(k) = CA*Ex(k)-CB/dlt*(Hy(k)-Hy(k-1))                                       ! Ex 迭代计算
        end do
        Ex_n(Zt_min) = Ex_n(Zt_min)+CB*Hy_in(Zt_min-1)                                      ! Ex 左端散射场-总场边界
        Ex_n(Zt_max) = Ex_n(Zt_max)-dt/(epsilon*dlt)*Hy_in(Zt_max)                          ! Ex 右端散射场-总场边界

        ! 一阶 Mur 吸收边界
        Ex_n(Z_max) = Ex(Z_max-1)+C1*(Ex_n(Z_max-1)-Ex(Z_max))                              ! Ex 右端吸收边界
        Ex_n(Z_min) = Ex(Z_min+1)+C1*(Ex_n(Z_min+1)-Ex(Z_min))                              ! Ex 左端吸收边界
        Ex_in = Exn_in
        Ex = Ex_n
        write(20, '(201F15.5)') (Hy(k), k = Z_min, Z_max-1)                                 ! 输出 Hy
    end do

    ! Hy 结果输出
    ! write(20, '(201F15.5)') (Hy(k), k = Z_min, Z_max-1)

    ! Z 坐标输出
    filename = "coordinate.txt"
    open(unit = 20, file = filename, status = "replace", action = "write", iostat=status, iomsg = msg)
    if ( status /= 0 ) stop "Error opening file name"
    write(20, '(201F15.5)') (k*dlt, k = Z_min, Z_max-1)


end program plane_wave
