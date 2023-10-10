# 柱坐标系二维 FDTD

[TOC]

柱坐标系 Maxwell 方程组

$$
\left\{\begin{aligned}
    &z_0\frac{\partial E_r}{\partial t}=-J_r-\sigma E_r-\frac{\partial H_\varphi}{\partial z}\\
    &\varepsilon_0\frac{\partial E_z}{\partial t}=-J_z-\sigma E_z+\frac{1}{r}H_\varphi+\frac{\partial H_\varphi}{\partial r}\\
    &\mu_0\frac{\partial H_\varphi}{\partial t}=\frac{\partial E_z}{\partial r}-\frac{\partial E_r}{\partial z}
\end{aligned}\right.
$$

二维 Yee 元胞离散

| 场分量 | $r$ 坐标 | $z$ 坐标 | 时间 $t$ 取样 |
| :--: | :--: | :--: | :--: |
| $E_r$ | $i+1/2$ | $k$ | $n$ |
| $E_z$ | $i$ | $k+1/2$ | $n$ |
| $H_\varphi$ | $i+1/2$ | $k+1/2$ | $n+1/2$ |

FDTD 递推公式

$$
\begin{aligned}
    E_r^{n+1}\left(i+\frac{1}{2},k\right)=&\frac{2\varepsilon-\sigma^{n+1/2}(i+1/2,k)\Delta t}{2\varepsilon+\sigma^{n+1/2}(i+1/2,k)\Delta t}E_r^n\left(i+\frac{1}{2},k\right)-\frac{2\Delta t}{2\varepsilon+\sigma^{n+1/2}(i+1/2,k)\Delta t}J_r^{n+1/2}\left(i+\frac{1}{2},k\right)\\
    &-\frac{2\Delta t}{2\varepsilon+\sigma^{n+1/2}(i+1/2,k)\Delta t}\cdot\frac{H^{n+1/2}_\varphi\left(i+\frac{1}{2},k+\frac{1}{2}\right)-H^{n+1/2}_\varphi\left(i+\frac{1}{2},k-\frac{1}{2}\right)}{\Delta z}\\\\
    E_z^{n+1}\left(i,k+\frac{1}{2}\right)=&\frac{2\varepsilon-\sigma^{n+1/2}(i,k+1/2)\Delta t}{2\varepsilon+\sigma^{n+1/2}(i,k+1/2)\Delta t}E_z^n\left(i,k+\frac{1}{2}\right)-\frac{2\Delta t}{2\varepsilon+\sigma^{n+1/2}(i,k+1/2)\Delta t}J_z^{n+1/2}\left(i,k+\frac{1}{2}\right)\\
    &+\frac{2\Delta t}{2\varepsilon+\sigma^{n+1/2}(i,k+1/2)\Delta t}\cdot\frac{\left(i+\frac{1}{2}\right)H_\varphi^{n+1/2}\left(i+\frac{1}{2},k+\frac{1}{2}\right)-\left(i-\frac{1}{2}\right)H_\varphi^{n+1/2}\left(i-\frac{1}{2},k+\frac{1}{2}\right)}{i\Delta r}\\\\
    H_\varphi^{n+1/2}\left(i+\frac{1}{2},k+\frac{1}{2}\right)=&H_\varphi^{n-1/2}\left(i+\frac{1}{2},k+\frac{1}{2}\right)+\frac{\Delta t}{\mu \Delta r}\left[E_z^n\left(i+1,k+\frac{1}{2}\right)-E_z^n\left(i,k+\frac{1}{2}\right)\right]\\
    &-\frac{\Delta t}{\mu\Delta z}\left[E_r^n\left(i+\frac{1}{2},k+1\right)-E_r^n\left(i+\frac{1}{2},k\right)\right]
\end{aligned}
$$

$r=0$ 处，电场轴向分量 $E_z$ 为奇异值，线性插值处理为

$$
E_z(0,k)=2E_z(1,k)-E_z(2,k)
$$

## 一阶 Mur 吸收边界

Mur 吸收边界的一阶近似

$$
\begin{aligned}
    \left(\frac{\partial}{\partial r}+\frac{1}{2r}+\frac{1}{v}\frac{\partial}{\partial t}\right)\psi&=0\hspace{5ex}r=r_{\max}\\
    \left(\frac{\partial}{\partial z}+\frac{1}{v}\frac{\partial}{\partial t}\right)\psi&=0\hspace{5ex}z=z_{\max}\\
    \left(\frac{\partial}{\partial z}-\frac{1}{v}\frac{\partial}{\partial t}\right)\psi&=0\hspace{5ex}z=z_{\min}
\end{aligned}
$$

离散形式

$$
\begin{aligned}
    E_z^{n+1}(i_{\min},k)=&\frac{4rv\Delta t-\Delta rv\Delta t+4r\Delta r}{4rv\Delta t+\Delta rv\Delta t+4r\Delta r}\cdot E_z^n(i_{\min}-1,k)+\frac{4rv\Delta t-\Delta rv\Delta t-4r\Delta r}{4rv\Delta t+\Delta rv\Delta t+4r\Delta r}\cdot E_z^{n+1}(i_{\min-1},k)\\
    &-\frac{4rv\Delta t+\Delta rv\Delta t-4r\Delta r}{4rv\Delta t+\Delta rv\Delta t+4r\Delta r}\cdot E_z^n(i_{\min},k)\\\\
    E_r^{n+1}(i,k_{\min})=&E_r^n(i,k_{\min}+1)+\frac{v\Delta t-\Delta z}{v\Delta t+\Delta z}\left[E_r^{n+1}(i,k_{\min}+1)-E_r^n(i,k_{\min})\right]\\\\
    E_r^{n+1}(i,k_{\max})=&E_r^n(i,k_{\max}-1)+\frac{v\Delta t-\Delta z}{v\Delta t+\Delta z}\left[E_r^{n+1}(i,k_{\max}-1)-E_r^n(i,k_{\max})\right]
\end{aligned}
$$

Mur 吸收边界的二阶近似有些复杂

$$
\begin{aligned}
    \left(\frac{\partial^2}{\partial r\partial t}+\frac{1}{2r}\frac{\partial}{\partial t}+\frac{1}{v}\frac{\partial^2}{\partial t^2}+\frac{v}{8r^2}-\frac{v}{2}\frac{\partial^2}{\partial z^2}\right)\psi&=0\hspace{5ex}r=r_{\max}\\
    \left(\frac{\partial^2}{\partial z\partial t}-\frac{1}{v}\frac{\partial^2}{\partial t^2}+\frac{v}{2}\frac{\partial^2}{\partial r^2}+\frac{v}{2r}\frac{\partial}{\partial r}\right)\psi&=0\hspace{5ex}z=z_{\min}\\
    \left(\frac{\partial^2}{\partial z\partial t}+\frac{1}{v}\frac{\partial^2}{\partial t^2}-\frac{v}{2}\frac{\partial^2}{\partial r^2}-\frac{v}{2r}\frac{\partial}{\partial r}\right)\psi&=0\hspace{5ex}z=z_{\max}
\end{aligned}
$$

## 初步结果

在 $r=0$ 附近加入正弦变化的激励电流（沿轴向），将产生柱面电波

<center>
    <img src="./line current.gif" width="350">
</center>

## 下一步计划

* 二阶 Mur 吸收边界
* 把 FDTD 程序写得更加规范
* 加入 Compton 电流源
* 考虑随时间、空间变化的电导率