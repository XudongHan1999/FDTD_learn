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
    integer,parameter :: Za_min = -5*20, Za_max = 5*20                                      ! 吸收边界指标
    integer,parameter :: Zt_min = -3*20, Zt_max = 3*20                                      ! 总场边界指标
    real, dimension(Za_min:Za_max) :: Ex_n = 0., Ex_p = 0.                                  ! Ex 数组
    real, dimension(Za_min:Za_max-1) :: Hy = 0.                                             ! Hy 数组
    ! real :: Exi_min = 0., Exi_max = 0., Hyi_min = 0., Hyi_max = 0.                          ! 入射电磁波
    real, dimension(Za_min:Za_max) :: Exi_min = 0., Hyi_min = 0.
    real, dimension(Za_min:Za_max) :: Exi_max = 0., Hyi_max = 0.


    real :: epsilon = 8.8541878e-12, mu = pi*4e-7, eta                                      ! 介电系数, 磁导系数, 
    real :: sigma = 0., sigma_m = 0.                                                        ! 电导率, 磁导率
    real :: CA, CB, CP, CQ                                                                  ! FDTD 迭代常数
    real :: C1                                                                              ! 吸收边界计算常数


    character(len=20) :: filename                                                           ! 输出文件
    integer :: status                                                                       ! 文件打开状态
    character(len=50) :: msg                                                                ! 文件打开信息

    dlt = lambda/30.
    dt = dlt/(2*c)

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
        ! if ( n*dt > (Zt_min-Za_min)*dlt/c ) then
        !     Hyi_min = sin(omega*(n*dt-(Zt_min-1-Za_min)*dlt/c))
        !     Exi_min = eta*sin(omega*((n+0.5)*dt-(Zt_min-Za_min)*dlt/c))
        ! end if
        ! if ( n*dt > (Zt_max-Za_min)*dlt/c ) then
        !     Hyi_max = sin(omega*((n)*dt-(Zt_max+0.25-Za_min)*dlt/c))
        !     Exi_max = eta*sin(omega*((n+0.5)*dt-(Zt_max+0.25-Za_min)*dlt/c))
        ! end if

        do k = Za_min, Za_max
            if ( n*dt > (k-Za_min)*dlt/c ) then
                Hyi_min(k) = sin(omega*((n)*dt-(k-Za_min)*dlt/c))
                Exi_min(k) = eta*sin(omega*((n+0.5)*dt-(k-Za_min)*dlt/c))
            end if
        end do
        do k = Za_min, Za_max
            if ( n*dt > (k-Za_min)*dlt/c ) then
                Hyi_max(k) = sin(omega*((n)*dt-(k-Za_min)*dlt/c))
                Exi_max(k) = eta*sin(omega*((n+0.5)*dt-(k-Za_min)*dlt/c))
            end if
        end do

        ! FDTD 迭代计算
        do k = Za_min, Za_max-1
            Hy(k) = CP*Hy(k)-CQ/dlt*(Ex_p(k+1)-Ex_p(k))                                     ! Hy 迭代计算
        end do
        ! Hy(Zt_min-1) = Hy(Zt_min-1)+dt/(mu*dlt)*Exi_min                                     ! Hy 左端散射场-总场边界
        ! Hy(Zt_max) = Hy(Zt_max)-dt/(mu*dlt)*Exi_max                                         ! Hy 右端散射场-总场边界
        Hy(Zt_min-1) = Hy(Zt_min-1)+dt/(mu*dlt)*Exi_min(Zt_min)
        Hy(Zt_max) = Hy(Zt_max)-dt/(mu*dlt)*Exi_min(Zt_max)


        do k = Za_min+1, Za_max-1
            Ex_n(k) = CA*Ex_p(k)-CB/dlt*(Hy(k)-Hy(k-1))                                     ! Ex 迭代计算
        end do
        ! Ex_n(Zt_min) = Ex_n(Zt_min)+dt/(epsilon*dlt)*Hyi_min                                ! Ex 左端散射场-总场边界
        ! Ex_n(Zt_max) = Ex_n(Zt_max)-dt/(epsilon*dlt)*Hyi_max                                ! Ex 右端散射场-总场边界
        Ex_n(Zt_min) = Ex_n(Zt_min)+dt/(epsilon*dlt)*Hyi_min(Zt_min-1)
        Ex_n(Zt_max) = Ex_n(Zt_max)-dt/(epsilon*dlt)*Hyi_min(Zt_max)

        ! 一阶 Mur 吸收边界
        Ex_n(Za_max) = Ex_p(Za_max-1)+C1*(Ex_n(Za_max-1)-Ex_p(Za_max))                      ! Ex 右端吸收边界
        Ex_n(Za_min) = Ex_p(Za_min+1)+C1*(Ex_n(Za_min+1)-Ex_p(Za_min))                      ! Ex 左端吸收边界
        Ex_p = Ex_n

        write(20, '(201F15.5)') (Hy(k), k = Za_min, Za_max-1)                               ! 输出 Hy
    end do

    ! Hy 结果输出
    write(20, '(201F15.5)') (Hy(k), k = Za_min, Za_max-1)

    ! Z 坐标输出
    filename = "coordinate.txt"
    open(unit = 20, file = filename, status = "replace", action = "write", iostat=status, iomsg = msg)
    if ( status /= 0 ) stop "Error opening file name"
    write(20, '(201F15.5)') (k*dlt, k = Za_min, Za_max-1)


end program plane_wave
