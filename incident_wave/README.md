# 散射场-总场边界

## 一维平面波入射

*一开始，我尝试直接写一个平面波的解析表达式加在总场边界上，后来发现总会有相位不匹配的问题。只能老老实实写一个一维 FDTD 入射波，然后通过总场边界条件实现入射波的引入*:smile:

<center class="half">
    <img src="./一维入射波_解析表达式.gif" width="400">
    <img src="./一维入射波_数值引入.gif" width="400">
</center>


## 二维平面波入射(TM波)
<center>
    <img src="./incident_2D.gif" width=400>
</center>

*有了一维平面波入射的编程基础，二维的平面波入射就很容易实现，但是在实现计算 TM 入射波时发现课本上的总场边界条件可能有一些小错误*:smile:

```Fortran
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
...

! 吸收边界处理
...
```
