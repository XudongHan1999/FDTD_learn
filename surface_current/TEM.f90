program TEM

    implicit none
    real, parameter :: pi = 3.1415926535                            ! pi
    real, parameter :: c = 3.0e8                                    ! 光速
    real,parameter :: f = 1.0e6                                     ! 频率
    real :: omega = 2*pi*f                                          ! 圆频率
    real :: lambda, T, delta, delta_t                               ! 波长, 周期
    real :: epsilon = 8.8541878e-12, mu = pi*4e-7, eta              ! 介电系数, 磁导系数
    real :: sigma = 0., sigma_m = 0.                                ! 电导率, 磁导率
    real :: CA, CB, CP, CQ                                          ! CA, CB, CP, CQ
    
    real, dimension(-5*40:5*40, 0:400) :: E_x = 0, H_y = 0          ! E_x, H_y 数组
    real, dimension(-5*40:5*40, 0:400) :: j = 0, M = 0              ! 电流, 磁流数组
    real, dimension(-5*40:5*40) :: z                                ! z坐标数组

    integer :: n, k, bound_l, bound_r                               ! n:时间循环变量, k:空间循环变量

    character(len=20) :: filename                                   ! 输出文件
    integer :: status                                               ! 文件打开状态
    character(len=50) :: msg                                        ! 文件打开信息

    lambda = c/f                                                    ! 波长
    T = 2*pi/omega                                                  ! 周期
    delta = lambda/40.                                              ! 空间步长
    delta_t = T/40.                                                 ! 时间步长
    bound_r = ubound(E_x, dim=1)                                    ! 右边界坐标
    bound_l = lbound(E_x, dim=1)                                    ! 左边界坐标
    forall (k = bound_l:bound_r)
        z(k) = k*delta                                              ! z 坐标
    end forall

    CA = (epsilon/delta_t-sigma/2.)/(epsilon/delta_t+sigma/2.)
    CB = 1/(epsilon/delta_t+sigma/2.)
    CP = (mu/delta_t-sigma_m/2)/(mu/delta_t+sigma_m/2.)
    CQ = 1/(mu/delta_t+sigma_m/2.)
    eta = sqrt(mu/epsilon)                                          ! 波阻抗

    ! 设置电流, 磁流
    forall (n = 0:400)
        j(0, n) = sin(omega*n*delta_t)
        ! M(0, n) = 0
        M(0, n) = -eta*sin(omega*n*delta_t)
    end forall

    ! FDTD 计算
    do n = 1, 400
        forall (k = bound_l:bound_r, (k /= bound_l) .or. (k /= bound_r))
            E_x(k, n) = CA*E_x(k, n-1)-CB*(H_y(k, n-1)-H_y(k-1, n-1))/delta-CB/delta*j(k, n)
            H_y(k, n) = CP*H_y(k, n-1)-CQ*(E_x(k+1, n)-E_x(k, n))/delta-CQ/delta*M(k,n)
        end forall
        ! 吸收边界
        E_x(bound_l, n) = E_x(bound_l+1, n-1) + (c*delta_t-delta)/(c*delta_t+delta)*(E_x(bound_l+1, n)-E_x(bound_l, n-1))
        E_x(bound_r, n) = E_x(bound_r-1, n-1) + (c*delta_t-delta)/(c*delta_t+delta)*(E_x(bound_r-1, n)-E_x(bound_r, n-1))
        H_y(bound_l, n) = H_y(bound_l+1, n-1) + (c*delta_t-delta)/(c*delta_t+delta)*(H_y(bound_l+1, n)-H_y(bound_l, n-1))
        H_y(bound_r, n) = H_y(bound_r-1, n-1) + (c*delta_t-delta)/(c*delta_t+delta)*(H_y(bound_r-1, n)-H_y(bound_r, n-1)) 
    end do

    ! 结果输出
    filename = "result.txt"
    open(unit = 20, file = filename, status = "replace", action = "write", iostat=status, iomsg = msg)
    write(*, *) msg
    if ( status /= 0 ) stop "Error opening file name"
    write(20, '(401F15.5, /)') (z(k), k = bound_l, bound_r)
    do n = 0, 400
        write(20, '(401F15.5, /)') (H_y(k, n), k = bound_l, bound_r)
    end do
    
end program TEM
