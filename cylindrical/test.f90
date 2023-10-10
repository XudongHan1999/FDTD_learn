program cylindrical_2D

    implicit none
    real, parameter :: pi = 3.1415926535                                ! pi
    real, parameter :: c = 3.0e8                                        ! 光速
    real, parameter :: f = 1.0e7                                        ! 频率
    real :: omega = 2*pi*f                                              ! 圆频率
    real :: lambda, T                                                   ! 波长, 周期
    real :: dr, dz, dt                                                  ! 空间步长, 时间步长

    real :: epsilon = 8.8541878e-12, mu = pi*4e-7                       ! 介电系数, 磁导系数

    real, dimension(0:10*40-1, -5*40:5*40, 0:500) :: E_r = 0.           ! 电场径向分量数组 E_r
    real, dimension(0:10*40, -5*40:5*40-1, 0:500) :: E_z = 0.           ! 电场轴向分量数组 E_z
    real, dimension(0:10*40-1, -5*40:5*40-1, 0:500) :: H_p = 0.         ! 磁场角向分量数组 H_p

    real, dimension(0:10*40, -5*40:5*40, 0:500) :: sigma = 0.           ! 电导率数组 sigma
    real, dimension(0:10*40-1, -5*40:5*40, 0:500) :: j_r = 0.           ! 径向电流数组 J_r
    real, dimension(0:10*40, -5*40:5*40-1, 0:500) :: j_z = 0.           ! 轴向电流数组 J_z

    real, dimension(0:10*40-1) :: rho                                   ! 径向坐标数组 rho
    real, dimension(-5*40:5*40-1) :: z                                  ! 轴向坐标数组 z
    integer :: n, i, k                                                  ! 时间循环变量 n, 空间循环变量 i, k
    integer :: bound_r, bound_t, bound_b                                ! 边界坐标 

    character(len=20) :: filename                                       ! 输出文件
    integer :: status                                                   ! 文件打开状态
    character(len=50) :: msg                                            ! 文件打开信息

    lambda = c/f                                                        ! 波长
    T = 2*pi/omega                                                      ! 周期
    dr = lambda/20.                                                     ! 径向空间步长
    dz = lambda/20.                                                     ! 轴向空间步长
    dt = T/30.                                                          ! 时间步长
    bound_r = ubound(E_z, dim=1)                                        ! 径向边界坐标
    bound_t = ubound(E_r, dim=2)                                        ! 轴向上边界坐标
    bound_b = lbound(E_r, dim=2)                                        ! 轴向下边界坐标

    do n = 0, 500
        do k = bound_b, bound_t-1
            J_z(2, k, n) = sin(omega*n*dt)
        end do
    end do

    do n = 1, 500
        ! FDTD 计算
        forall (i = 0:bound_r-1, k = bound_b:bound_t-1)
            H_p(i, k, n) = H_p(i, k, n-1)+(dt/(mu*dr))*(E_z(i+1, k, n-1)-E_z(i, k, n-1)) &
                & -(dt/(mu*dz))*(E_r(i, k+1, n-1)-E_r(i, k, n-1))
        end forall
        forall (i = 0:bound_r-1, k = bound_b+1:bound_t-1)
            E_r(i, k, n) = (2*epsilon-sigma(i, k, n)*dt)/(2*epsilon+sigma(i, k, n)*dt)*E_r(i, k, n-1) &
                & -(2*dt)/(2*epsilon+sigma(i, k, n)*dt)*J_r(i, k, n)/(dr*dz) &
                & -(2*dt)/(2*epsilon+sigma(i, k, n)*dt)*(H_p(i, k, n)-H_p(i, k-1, n))/dz
        end forall        
        forall (i = 1:bound_r-1, k = bound_b:bound_t-1)
            E_z(i, k, n) = (2*epsilon-sigma(i, k, n)*dt)/(2*epsilon+sigma(i, k, n)*dt)*E_z(i, k, n-1) &
                & -(2*dt)/(2*epsilon+sigma(i, k, n)*dt)*J_z(i, k, n)/(dr*dz) &
                & +(2*dt)/(2*epsilon+sigma(i, k, n)*dt)*(H_p(i, k, n)-H_p(i-1, k, n))/dr &
                & +dt/(2*epsilon+sigma(i, k, n)*dt)/i*(H_p(i, k, n)+H_p(i-1, k, n))/dr
        end forall

        ! 一阶 Mur 吸收边界
        do k = bound_b+1, bound_t-1
            E_z(0, k, n) = 2*E_z(1, k, n)-E_z(2, k, n)
            E_z(bound_r, k, n) = (4*(bound_r-0.5)*dr*c*dt-dr*c*dt+4*(bound_r-0.5)*dr*dr)/ &
                & (4*(bound_r-0.5)*dr*c*dt+dr*c*dt+4*(bound_r-0.5)*dr*dr)*E_z(bound_r-1, k, n-1)+ &
                & (4*(bound_r-0.5)*dr*c*dt-dr*c*dt-4*(bound_r-0.5)*dr*dr)/ &
                & (4*(bound_r-0.5)*dr*c*dt+dr*c*dt+4*(bound_r-0.5)*dr*dr)*E_z(bound_r-1, k, n)+ &
                & (4*(bound_r-0.5)*dr*c*dt+dr*c*dt-4*(bound_r-0.5)*dr*dr)/ &
                & (4*(bound_r-0.5)*dr*c*dt+dr*c*dt+4*(bound_r-0.5)*dr*dr)*E_z(bound_r, k, n-1)
        end do
        do i = 1, bound_r-1
            E_r(i, bound_t, n) = E_r(i, bound_t-1, n-1)+(c*dt-dz)/(c*dt+dz)*(E_r(i, bound_t-1, n)-E_r(i, bound_t, n-1))
            E_r(i, bound_b, n) = E_r(i, bound_b+1, n-1)+(c*dt-dz)/(c*dt+dz)*(E_r(i, bound_b+1, n)-E_r(i, bound_b, n-1))
        end do

        ! 二维角点处理
        E_r(0, bound_b, n) = E_r(1, bound_b+1, n-1) &
            & +(c*dt-sqrt(dr**2+dz**2))/(c*dt+sqrt(dr**2+dz**2))*(E_r(1, bound_b+1, n)-E_r(0, bound_b, n-1))
        E_r(0, bound_t, n) = E_r(1, bound_t-1, n-1) &
            & +(c*dt-sqrt(dr**2+dz**2))/(c*dt+sqrt(dr**2+dz**2))*(E_r(1, bound_t-1, n)-E_r(0, bound_t, n-1))
        E_z(0, bound_b, n) = E_z(1, bound_b+1, n-1) &
            & +(c*dt-sqrt(dr**2+dz**2))/(c*dt+sqrt(dr**2+dz**2))*(E_z(1, bound_b+1, n)-E_z(0, bound_b, n-1))
        E_z(bound_r, bound_b, n) = E_z(bound_r-1, bound_b+1, n-1) &
            & +(c*dt-sqrt(dr**2+dz**2))/(c*dt+sqrt(dr**2+dz**2))*(E_z(bound_r-1, bound_b+1, n)-E_z(bound_r, bound_b, n-1))        
    end do

    ! rho 方向坐标计算
    do i = 0, bound_r
        rho(i) = i*dr
    end do

    ! 结果输出
    filename = "coordinate.txt"
    open(unit = 20, file = filename, status = "replace", action = "write", iostat=status, iomsg = msg)
    if ( status /= 0 ) stop "Error opening file name"
    write(20, '(401F15.5, /)') (rho(i), i = 0, bound_r)
    
    filename = "result.txt"
    open(unit = 20, file = filename, status = "replace", action = "write", iostat=status, iomsg = msg)
    if ( status /= 0 ) stop "Error opening file name"

    do n = 0, 500
        do k = bound_b, bound_t-1
            write(20, '(401F15.5, /)') (E_z(i, k, n-1), i = 0, bound_r)
        end do
    end do

end program cylindrical_2D