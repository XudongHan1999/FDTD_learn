program TE

    implicit none
    real, parameter :: pi = 3.1415926535                                ! pi
    real, parameter :: c = 3.0e8                                        ! 光速
    real, parameter :: f = 1.0e7                                        ! 频率
    real :: omega = 2*pi*f                                              ! 圆频率
    real :: lambda = c/f                                                ! 波长
    real :: dlt, dt                                                     ! dlt: 空间步长, dt: 时间步长

    integer :: n, i, j                                                  ! n: 时间迭代变量, i,j: 空间循环变量
    integer,parameter :: t_max = 400                                    ! t_max: 时间范围指标   
    integer,parameter :: x_min = -5*20, x_max = 5*20, &                 ! 空间范围指标
        & y_min = -5*20, y_max = 5*20


    real, dimension(x_min:x_max, y_min:y_max-1, 0:t_max) :: E_x = 0.    ! E_x 数组
    real, dimension(x_min:x_max-1, y_min:y_max, 0:t_max) :: E_y = 0.    ! E_y 数组  
    real, dimension(x_min:x_max, y_min:y_max, 0:t_max) :: H_z = 0.      ! H_x 数组
    real, dimension(x_min:x_max, y_min:y_max, 0:t_max) :: M_z = 0.      ! 磁流数组

    real :: epsilon = 8.8541878e-12, mu = pi*4e-7                       ! 介电系数, 磁导系数
    real :: sigma = 0, sigma_m = 0.                                     ! 电导率, 磁导率
    real :: CA, CB, CP, CQ                                              ! FDTD 迭代常数
    real :: C1, C2, C3                                                  ! 吸收边界计算常数

    character(len=20) :: filename                                       ! 输出文件
    integer :: status                                                   ! 文件打开状态
    character(len=50) :: msg                                            ! 文件打开信息


    dlt = lambda/20.
    dt = dlt/(2*c)

    CA = (epsilon/dt-sigma/2.)/(epsilon/dt+sigma/2.)
    CB = 1./(epsilon/dt+sigma/2.)
    CP = (mu/dt-sigma_m/2.)/(mu/dt+sigma_m/2.)
    CQ = 1./(mu/dt+sigma_m/2.)

    C1 = (c*dt-dlt)/(c*dt+dlt)
    C2 = (c**2*epsilon*dt)/(2*(c*dt+dlt))
    C3 = (c*dt-sqrt(2.)*dlt)/(c*dt+sqrt(2.)*dlt)

    ! 线电流源
    do n = 0, t_max
        M_z(0, 0, n) = sqrt(mu/epsilon)*sin(omega*n*dt)
    end do

    do n = 1, t_max
        ! FDTD 计算
        forall (i = x_min:x_max, j = y_min:y_max-1)
            E_x(i, j, n) = CA*E_x(i, j, n-1)+CB/dlt*(H_z(i, j+1, n-1)-H_z(i, j, n-1))          
        end forall
        forall (i = x_min:x_max-1, j = y_min:y_max)
            E_y(i, j, n) = CA*E_y(i, j, n-1)-CB/dlt*(H_z(i+1, j, n-1)-H_z(i, j, n-1))            
        end forall
        forall (i = x_min+1:x_max-1, j = y_min+1:y_max-1)
            H_z(i, j, n) = CP*H_z(i, j, n-1)-CQ/dlt*((E_y(i, j, n)-E_y(i-1, j, n))-(E_x(i, j, n)-E_x(i, j-1, n))) &
                & +CQ*M_z(i, j, n)/(dlt**2)
        end forall

        ! 二维边界: 二阶 Mur 吸收边界
        do j = y_min+1, y_max-1
            H_z(x_min, j, n) = H_z(x_min+1, j, n-1)+C1*(H_z(x_min+1, j, n)-H_z(x_min, j, n-1)) &
                & +C2*(E_x(x_min, j, n)+E_x(x_min+1, j, n)-E_x(x_min, j-1, n)-E_x(x_min+1, j-1, n))
            H_z(x_max, j, n) = H_z(x_max-1, j, n-1)+C1*(H_z(x_max-1, j, n)-H_z(x_max, j, n-1)) &
                & +C2*(E_x(x_max, j, n)+E_x(x_max-1, j, n)-E_x(x_max, j-1, n)-E_x(x_max-1, j-1, n))           
        end do
        do i = x_min+1, x_max-1
            H_z(i, y_min, n) = H_z(i, y_min+1, n-1)+C1*(H_z(i, y_min+1, n)-H_z(i, y_min, n-1)) &
                & -C2*(E_y(i, y_min, n)+E_y(i, y_min+1, n)-E_y(i-1, y_min, n)-E_y(i-1, y_min+1, n))
            H_z(i, y_max, n) = H_z(i, y_max-1, n-1)+C1*(H_z(i, y_max-1, n)-H_z(i, y_max, n-1)) &
                & -C2*(E_y(i, y_max, n)+E_y(i, y_max-1, n)-E_y(i-1, y_max, n)-E_y(i-1, y_max-1, n))
        end do

        ! 角点处理: 一阶 Mur 吸收边界
        H_z(x_min, y_min, n) = H_z(x_min+1, y_min+1, n-1)+C3*(H_z(x_min+1, y_min+1, n)-H_z(x_min, y_min, n-1))
        H_z(x_max, y_min, n) = H_z(x_max-1, y_min+1, n-1)+C3*(H_z(x_max-1, y_min+1, n)-H_z(x_max, y_min, n-1))
        H_z(x_min, y_max, n) = H_z(x_min+1, y_max-1, n-1)+C3*(H_z(x_min+1, y_max-1, n)-H_z(x_min, y_max, n-1))
        H_z(x_max, y_max, n) = H_z(x_max-1, y_max-1, n-1)+C3*(H_z(x_max-1, y_max-1, n)-H_z(x_max, y_max, n-1))
    end do

    ! 结果输出
    filename = "Hz.txt"
    open(unit = 20, file = filename, status = "replace", action = "write", iostat=status, iomsg = msg)
    if ( status /= 0 ) stop "Error opening file name"
    do j = y_min, y_max
        write(20, '(201F15.5)') (H_z(i, j, t_max), i = x_min, x_max)
    end do

end program TE
