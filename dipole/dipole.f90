program dipole

    implicit none
    real, parameter :: pi = 3.1415926535                                                    ! pi
    real, parameter :: c = 3.0e8                                                            ! 光速
    real, parameter :: f = 1.0e7                                                            ! 频率

    real :: omega = 2*pi*f                                                                  ! 圆频率
    real :: lambda = c/f                                                                    ! 波长
    real :: dlt, dt                                                                         ! dlt: 空间步长, dt: 时间步长

    integer :: n, i, j, k                                                                   ! n: 时间迭代变量, i,j,k: 空间循环变量
    integer,parameter :: t_max = 100                                                        ! t_max: 时间范围指标   
    integer,parameter :: x_min = -2*15, x_max = 2*15, &                                     ! 空间范围指标
        & y_min = -2*15, y_max = 2*15, z_min = -2*15, z_max = 2*15
    real, dimension(x_min:x_max-1, y_min:y_max, z_min:z_max, -1:t_max) :: E_x = 0.          ! E_x 数组
    real, dimension(x_min:x_max, y_min:y_max-1, z_min:z_max, -1:t_max) :: E_y = 0.          ! E_y 数组
    real, dimension(x_min:x_max, y_min:y_max, z_min:z_max-1, -1:t_max) :: E_z = 0.          ! E_z 数组
    real, dimension(x_min:x_max, y_min:y_max-1, z_min:z_max-1, 0:t_max) :: H_x = 0.         ! H_x 数组
    real, dimension(x_min:x_max-1, y_min:y_max, z_min:z_max-1, 0:t_max) :: H_y = 0.         ! H_y 数组
    real, dimension(x_min:x_max-1, y_min:y_max-1, z_min:z_max, 0:t_max) :: H_z = 0.         ! H_z 数组
    real, dimension(x_min:x_max, y_min:y_max, z_min:z_max-1, -1:t_max) :: J_z = 0.          ! J_z 数组

    real :: epsilon = 8.8541878e-12, mu = pi*4e-7                                           ! 介电系数, 磁导系数
    real :: sigma = 0., sigma_m = 0.                                                        ! 电导率, 磁导率
    real :: CA, CB, CP, CQ                                                                  ! FDTD 迭代常数
    real :: C1, C2, C3, C4                                                                  ! 吸收边界常数

    character(len=20) :: filename                                                           ! 输出文件名
    integer :: status                                                                       ! 文件打开状态
    character(len=50) :: msg                                                                ! 文件错误信息

    dlt = lambda/15.
    dt = dlt/(2*c)

    CA = (epsilon/dt-sigma/2.)/(epsilon/dt+sigma/2.)
    CB = 1./(epsilon/dt+sigma/2.)
    CP = (mu/dt-sigma_m/2.)/(mu/dt+sigma_m/2.)
    CQ = 1./(mu/dt+sigma_m/2.)

    C1 = (c*dt-dlt)/(c*dt+dlt)
    C2 = 2.*(c*dt-dlt)/dlt
    C3 = (c*dt)**2./(2.*dlt*(c*dt+dlt))
    C4 = (c*dt-sqrt(2.)*dlt)/(c*dt+sqrt(2.)*dlt)

    do n = -1, t_max
        J_z(0, 0, 0, n) = 0.01*sin(omega*n*dt)/(dlt)**3                                     ! J_z 电流



    end do

    do n = 1, t_max
        ! FDTD 迭代计算
        forall (i = x_min:x_max, j = y_min:y_max-1, k = z_min:z_max-1)
            H_x(i, j, k, n) = CP*H_x(i, j, k, n-1) &
                & -CQ/dlt*((E_z(i, j+1, k, n-1)-E_z(i, j, k, n-1))-(E_y(i, j, k+1, n-1)-E_y(i, j, k, n-1)))
        end forall
        forall (i = x_min:x_max-1, j = y_min:y_max, k = z_min:z_max-1)
            H_y(i, j, k, n) = CP*H_y(i, j, k, n-1) &
                & -CQ/dlt*((E_x(i, j, k+1, n-1)-E_x(i, j, k, n-1))-(E_z(i+1, j, k, n-1)-E_z(i, j, k, n-1)))
        end forall
        forall (i = x_min:x_max-1, j = y_min:y_max-1, k = z_min:z_max)
            H_z(i, j, k, n) = CP*H_z(i, j, k, n-1) &
                & -CQ/dlt*((E_y(i+1, j, k, n-1)-E_y(i, j, k, n-1))-(E_x(i, j+1, k, n-1)-E_x(i, j, k, n-1)))
        end forall
        forall (i = x_min:x_max-1, j = y_min+1:y_max-1, k = z_min+1:z_max-1)
            E_x(i, j, k, n) = CA*E_x(i, j, k, n-1) &
                & +CB/dlt*((H_z(i, j, k, n)-H_z(i, j-1, k, n))-(H_y(i, j, k, n)-H_y(i, j, k-1, n)))
        end forall
        forall (i = x_min+1:x_max-1, j = y_min:y_max-1, k = z_min+1:z_max-1)
            E_y(i, j, k, n) = CA*E_y(i, j, k, n-1) &
                & +CB/dlt*((H_x(i, j, k, n)-H_x(i, j, k-1, n))-(H_z(i, j, k, n)-H_z(i-1, j, k, n)))
        end forall
        forall (i = x_min+1:x_max-1, j = y_min+1:y_max-1, k = z_min:z_max-1)
            E_z(i, j, k, n) = CA*E_z(i, j, k, n-1) &
                & +CB/dlt*((H_y(i, j, k, n)-H_y(i-1, j, k, n))-(H_x(i, j, k, n)-H_x(i, j-1, k, n)))
        end forall


        ! 加入电流源
        forall (i = x_min+1:x_max-1, j = y_min+1:y_max-1, k = z_min:z_max-1)
            E_z(i, j, k, n) = E_z(i, j, k, n) - J_z(i, j, k, n)/(epsilon/dt+sigma/2.)
        end forall


        ! 面处理: 二阶 Mur 吸收边界
        forall (i = x_min+1:x_max-2, j = y_min+1:y_max-1)
            E_x(i, j, z_max, n) = -E_x(i, j, z_max-1, n-2)+C1*(E_x(i, j, z_max-1, n)+E_x(i, j, z_max, n-2)) &
                & -C2*(E_x(i, j, z_max-1, n-1)+E_x(i, j, z_max, n-1)) &
                & +C3*(E_x(i-1, j, z_max, n-1)+E_x(i+1, j, z_max, n-1)+E_x(i, j-1, z_max, n-1)+E_x(i, j+1, z_max, n-1) &
                & +E_x(i-1, j, z_max-1, n-1)+E_x(i+1, j, z_max-1, n-1)+E_x(i, j-1, z_max-1, n-1)+E_x(i, j+1, z_max-1, n-1))
            E_x(i, j, z_min, n) = -E_x(i, j, z_min+1, n-2)+C1*(E_x(i, j, z_min+1, n)+E_x(i, j, z_min, n-2)) &
                & -C2*(E_x(i, j, z_min+1, n-1)+E_x(i, j, z_min, n-1)) &
                & +C3*(E_x(i-1, j, z_min, n-1)+E_x(i+1, j, z_min, n-1)+E_x(i, j-1, z_min, n-1)+E_x(i, j+1, z_min, n-1) &
                & +E_x(i-1, j, z_min+1, n-1)+E_x(i+1, j, z_min+1, n-1)+E_x(i, j-1, z_min+1, n-1)+E_x(i, j+1, z_min+1, n-1))            
        end forall
        forall (i = x_min+1:x_max-2, k = z_min+1:z_max-1)
            E_x(i, y_max, k, n) = -E_x(i, y_max-1, k, n-2)+C1*(E_x(i, y_max-1, k, n)+E_x(i, y_max, k, n-2)) &
                & -C2*(E_x(i, y_max-1, k, n-1)+E_x(i, y_max, k, n-1)) &
                & +C3*(E_x(i-1, y_max, k, n-1)+E_x(i+1, y_max, k, n-1)+E_x(i, y_max, k-1, n-1)+E_x(i, y_max, k+1, n-1) &
                & +E_x(i-1, y_max-1, k, n-1)+E_x(i+1, y_max-1, k, n-1)+E_x(i, y_max-1, k-1, n-1)+E_x(i, y_max-1, k+1, n-1))
            E_x(i, y_min, k, n) = -E_x(i, y_min+1, k, n-2)+C1*(E_x(i, y_min+1, k, n)+E_x(i, y_min, k, n-2)) &
                & -C2*(E_x(i, y_min+1, k, n-1)+E_x(i, y_min, k, n-1)) &
                & +C3*(E_x(i-1, y_min, k, n-1)+E_x(i+1, y_min, k, n-1)+E_x(i, y_min, k-1, n-1)+E_x(i, y_min, k+1, n-1) &
                & +E_x(i-1, y_min+1, k, n-1)+E_x(i+1, y_min+1, k, n-1)+E_x(i, y_min+1, k-1, n-1)+E_x(i, y_min+1, k+1, n-1))
        end forall
        forall (i = x_min+1:x_max-1, j = y_min+1:y_max-2)
            E_y(i, j, z_max, n) = -E_y(i, j, z_max-1, n-2)+C1*(E_y(i, j, z_max-1, n)+E_y(i, j, z_max, n-2)) &
                & -C2*(E_y(i, j, z_max-1, n-1)+E_y(i, j, z_max, n-1)) &
                & +C3*(E_y(i-1, j, z_max, n-1)+E_y(i+1, j, z_max, n-1)+E_y(i, j-1, z_max, n-1)+E_y(i, j+1, z_max, n-1) &
                & +E_y(i-1, j, z_max-1, n-1)+E_y(i+1, j, z_max-1, n-1)+E_y(i, j-1, z_max-1, n-1)+E_y(i, j+1, z_max-1, n-1))
            E_y(i, j, z_min, n) = -E_y(i, j, z_min+1, n-2)+C1*(E_y(i, j, z_min+1, n)+E_y(i, j, z_min, n-2)) &
                & -C2*(E_y(i, j, z_min+1, n-1)+E_y(i, j, z_min, n-1)) &
                & +C3*(E_y(i-1, j, z_min, n-1)+E_y(i+1, j, z_min, n-1)+E_y(i, j-1, z_min, n-1)+E_y(i, j+1, z_min, n-1) &
                & +E_y(i-1, j, z_min+1, n-1)+E_y(i+1, j, z_min+1, n-1)+E_y(i, j-1, z_min+1, n-1)+E_y(i, j+1, z_min+1, n-1))
        end forall
        forall (j = y_min+1:y_max-2, k = z_min+1:z_max-1)
            E_y(x_max, j, k, n) = -E_y(x_max-1, j, k, n-2)+C1*(E_y(x_max-1, j, k, n)+E_y(x_max, j, k, n-2)) &
                & -C2*(E_y(x_max-1, j, k, n-1)+E_y(x_max, j, k, n-1)) &
                & +C3*(E_y(x_max, j-1, k, n-1)+E_y(x_max, j+1, k, n-1)+E_y(x_max, j, k-1, n-1)+E_y(x_max, j, k+1, n-1) &
                & +E_y(x_max-1, j-1, k, n-1)+E_y(x_max-1, j+1, k, n-1)+E_y(x_max-1, j, k-1, n-1)+E_y(x_max-1, j, k+1, n-1))
            E_y(x_min, j, k, n) = -E_y(x_min+1, j, k, n-2)+C1*(E_y(x_min+1, j, k, n)+E_y(x_min, j, k, n-2)) &
                & -C2*(E_y(x_min+1, j, k, n-1)+E_y(x_min, j, k, n-1)) &
                & +C3*(E_y(x_min, j-1, k, n-1)+E_y(x_min, j+1, k, n-1)+E_y(x_min, j, k-1, n-1)+E_y(x_min, j, k+1, n-1) &
                & +E_y(x_min+1, j-1, k, n-1)+E_y(x_min+1, j+1, k, n-1)+E_y(x_min+1, j, k-1, n-1)+E_y(x_min+1, j, k+1, n-1))
        end forall
        forall (j = y_min+1:y_max-1, k = z_min+1:z_max-2)
            E_z(x_max, j, k, n) = -E_z(x_max-1, j, k, n-2)+C1*(E_z(x_max-1, j, k, n)+E_z(x_max, j, k, n-2)) &
                & -C2*(E_z(x_max-1, j, k, n-1)+E_z(x_max, j, k, n-1)) &
                & +C3*(E_z(x_max, j-1, k, n-1)+E_z(x_max, j+1, k, n-1)+E_z(x_max, j, k-1, n-1)+E_z(x_max, j, k+1, n-1) &
                & +E_z(x_max-1, j-1, k, n-1)+E_z(x_max-1, j+1, k, n-1)+E_z(x_max-1, j, k-1, n-1)+E_z(x_max-1, j, k+1, n-1))
            E_z(x_min, j, k, n) = -E_z(x_min+1, j, k, n-2)+C1*(E_z(x_min+1, j, k, n)+E_z(x_min, j, k, n-2)) &
                & -C2*(E_z(x_min+1, j, k, n-1)+E_z(x_min, j, k, n-1)) &
                & +C3*(E_z(x_min, j-1, k, n-1)+E_z(x_min, j+1, k, n-1)+E_z(x_min, j, k-1, n-1)+E_z(x_min, j, k+1, n-1) &
                & +E_z(x_min+1, j-1, k, n-1)+E_z(x_min+1, j+1, k, n-1)+E_z(x_min+1, j, k-1, n-1)+E_z(x_min+1, j, k+1, n-1))
        end forall
        forall (i = x_min+1:x_max-1, k = z_min+1:z_max-2)
            E_z(i, y_max, k, n) = -E_z(i, y_max-1, k, n-2)+C1*(E_z(i, y_max-1, k, n)+E_z(i, y_max, k, n-2)) &
                & -C2*(E_z(i, y_max-1, k, n-1)+E_z(i, y_max, k, n-1)) &
                & +C3*(E_z(i-1, y_max, k, n-1)+E_z(i+1, y_max, k, n-1)+E_z(i, y_max, k-1, n-1)+E_z(i, y_max, k+1, n-1) &
                & +E_z(i-1, y_max-1, k, n-1)+E_z(i+1, y_max-1, k, n-1)+E_z(i, y_max-1, k-1, n-1)+E_z(i, y_max-1, k+1, n-1))
            E_z(i, y_min, k, n) = -E_z(i, y_min+1, k, n-2)+C1*(E_z(i, y_min+1, k, n)+E_z(i, y_min, k, n-2)) &
                & -C2*(E_z(i, y_min+1, k, n-1)+E_z(i, y_min, k, n-1)) &
                & +C3*(E_z(i-1, y_min, k, n-1)+E_z(i+1, y_min, k, n-1)+E_z(i, y_min, k-1, n-1)+E_z(i, y_min, k+1, n-1) &
                & +E_z(i-1, y_min+1, k, n-1)+E_z(i+1, y_min+1, k, n-1)+E_z(i, y_min+1, k-1, n-1)+E_z(i, y_min+1, k+1, n-1))
        end forall

        ! 棱边处理: 一阶 Mur 吸收边界
        do i = x_min, x_max-1
            E_x(i, y_max, z_max, n) = E_x(i, y_max-1, z_max-1, n-1)+C4*(E_x(i, y_max-1, z_max-1, n)-E_x(i, y_max, z_max, n-1))
            E_x(i, y_max, z_min, n) = E_x(i, y_max-1, z_min+1, n-1)+C4*(E_x(i, y_max-1, z_min+1, n)-E_x(i, y_max, z_min, n-1))
            E_x(i, y_min, z_max, n) = E_x(i, y_min+1, z_max-1, n-1)+C4*(E_x(i, y_min+1, z_max-1, n)-E_x(i, y_min, z_max, n-1))
            E_x(i, y_min, z_min, n) = E_x(i, y_min+1, z_min+1, n-1)+C4*(E_x(i, y_min+1, z_min+1, n)-E_x(i, y_min, z_min, n-1))
        end do
        do j = y_min, y_max-1
            E_y(x_max, i, z_max, n) = E_y(x_max-1, i, z_max-1, n-1)+C4*(E_y(x_max-1, i, z_max-1, n)-E_y(x_max, i, z_max, n-1))
            E_y(x_max, i, z_min, n) = E_y(x_max-1, i, z_min+1, n-1)+C4*(E_y(x_max-1, i, z_min+1, n)-E_y(x_max, i, z_min, n-1))
            E_y(x_min, i, z_max, n) = E_y(x_min+1, i, z_max-1, n-1)+C4*(E_y(x_min+1, i, z_max-1, n)-E_y(x_min, i, z_max, n-1))
            E_y(x_min, i, z_min, n) = E_y(x_min+1, i, z_min+1, n-1)+C4*(E_y(x_min+1, i, z_min+1, n)-E_y(x_min, i, z_min, n-1))
        end do
        do k = z_min, z_max-1
            E_z(x_max, y_max, k, n) = E_z(x_max-1, y_max-1, k, n-1)+C4*(E_z(x_max-1, y_max-1, k, n)-E_z(x_max, y_max, k, n-1))
            E_z(x_max, y_min, k, n) = E_z(x_max-1, y_min+1, k, n-1)+C4*(E_z(x_max-1, y_min+1, k, n)-E_z(x_max, y_min, k, n-1))
            E_z(x_min, y_max, k, n) = E_z(x_min+1, y_max-1, k, n-1)+C4*(E_z(x_min+1, y_max-1, k, n)-E_z(x_min, y_max, k, n-1))
            E_z(x_min, y_min, k, n) = E_z(x_min+1, y_min+1, k, n-1)+C4*(E_z(x_min+1, y_min+1, k, n)-E_z(x_min, y_min, k, n-1))
        end do
    end do

    ! 结果输出
    filename = "coordinate.txt"
    open(unit = 20, file = filename, status = "replace", action = "write", iostat=status, iomsg = msg)
    if ( status /= 0 ) stop "Error opening file name"
    write(20, '(60F15.5)') (i*dlt, i = x_min, x_max-1)
    write(20, '(60F15.5)') (j*dlt, j = y_min, y_max-1)
    write(20, '(60F15.5)') (k*dlt, k = z_min, z_max-1)

    filename = "E_y.txt"
    open(unit = 20, file = filename, status = "replace", action = "write", iostat=status, iomsg = msg)
    if ( status /= 0 ) stop "Error opening file name"
    do k = z_min, z_max-1
        write(20, '(60F15.5)') (E_y(0, j, k, 70), j = y_min, y_max-1)
    end do

    filename = "E_z.txt"
    open(unit = 20, file = filename, status = "replace", action = "write", iostat=status, iomsg = msg)
    if ( status /= 0 ) stop "Error opening file name"
    do k = z_min, z_max-1
        write(20, '(60F15.5)') (E_z(0, j, k, 70), j = y_min, y_max-1)
    end do

    filename = "E.txt"
    open(unit = 20, file = filename, status = "replace", action = "write", iostat=status, iomsg = msg)
    if ( status /= 0 ) stop "Error opening file name"
    do k = z_min, z_max-1
        do j = y_min, y_max-1
            do i = x_min, x_max-1
                write(20, '(3F15.5)') E_x(i, j, k, t_max), E_y(i, j, k, t_max), E_z(i, j, k, t_max)
            end do
        end do
    end do

end program dipole
