# GSIS 中 Maxwell 漫反射边界通量构造

> 适用范围补充（2026-10-08）：本文原有推导针对漫反射壁面。新 FV 正激波 demo 使用左右远场和上下镜面对称，没有漫反射壁温或反射密度条件。两者共享“半空间 Maxwell 通量加冻结介观差值”的思路，但物理边界不能混用。对应源码与完整说明见 [NSW 边界指南](D:/2025work/dg/gsis_nsw_matlab_demo/maxwell_gsis_boundary_condition.md) 和 [CIS、NS、GSIS 子程序指南](D:/2025work/dg/gsis_nsw_matlab_demo/SOLVER_GUIDE.md)。

## NSW demo 的边界与调用流程补充

以下约定以计算域外法向为准：v·n<0 进入计算域，v·n>0 离开计算域。原壁面推导中的“朝向壁面的入射”需按法向符号理解，不要直接套用叫法。

- 左右远场：dvm_boundary 填上、下游给定 Maxwell 虚拟分布，dvm_convection 仅在进入计算域的半空间使用它；离开域的分布取内场。上下边界按 vy 反号镜面反射，不使用漫反射密度公式。
- 纯 NS：nsw_closure 的 half 计算当前内场 Maxwell 出射通量与给定远场 Maxwell 入射通量之和；ns_boundary 在上下对称面使用压力通量。
- GSIS：dvm_hot 定义 Dfar=Fkin-FMaxwell(qkin)。宏观内迭代冻结 Dfar，使用 FMaxwell(qNS当前)+Dfar。恢复到 Fkin 的恒等式只在 qNS=qkin 处成立，且这里只应用于左右远场。
- 边界差值按单位面长保存，在 ns_residual 中统一乘面长。内部面的格心 HoT 由 ns_viscous 加入，不在远场边界重复添加。
- 本版前10步为 CIS，第10步只初始化一次宏观状态，第11步起每轮调用100步 monoNS；宏观状态持续保留，下一轮反馈碰撞目标，不进行 VDF lift。

上下对称面目前没有显式偏应力修正，验证范围是一维正交条带。边界恢复检查通过不等于整体离散固定点一致；当前固定 HoT 测试仍有约1.57%的宏观场变化。具体公式、数组和验证界限以链接的 demo 文档为准。以下保留原漫反射壁面说明。

本小节整理当前 GSIS 程序中壁面边界通量的构造逻辑，重点说明四件事：

1. 介观漫反射边界中，入射分布函数和反射分布函数如何构造；
2. 介观分布函数如何通过半速度空间积分得到宏观边界通量；
3. Maxwell 平衡态边界通量如何用半空间解析积分得到；
4. GSIS 中边界高阶项如何进入宏观迭代通量。

## 0. 符号约定

记边界面中心、相邻内场单元中心、壁面外法向、分子速度分别为

$$
\boldsymbol{x}_b,\qquad
\boldsymbol{x}_i,\qquad
\boldsymbol{n},\qquad
\boldsymbol{v}_k .
$$

法向速度定义为

$$
v_{n,k}
=
\boldsymbol{v}_k\cdot\boldsymbol{n}.
$$

这里的法向量从内场单元指向边界外侧。因此速度空间可以分为入射壁面的半空间和壁面反射回来的半空间：

$$
\Xi^+
=
\left\{
\boldsymbol{v}:
\boldsymbol{v}\cdot\boldsymbol{n}>0
\right\},
\qquad
\Xi^-
=
\left\{
\boldsymbol{v}:
\boldsymbol{v}\cdot\boldsymbol{n}<0
\right\}.
$$

其中正半空间表示气体内部射向壁面的分子，负半空间表示由壁面反射回计算域的分子。边界宏观通量写成半空间分裂形式：

$$
\boldsymbol{F}_b
=
\boldsymbol{F}_b^+
+
\boldsymbol{F}_b^- .
$$

## 1. 介观漫反射的入射与反射边界条件

在 DVM 介观层面，边界面上的分布函数按照分子速度方向分裂。入射壁面的分子由内场分布函数外推得到，反射回来的分子由壁面漫反射 Maxwell 分布给出：

$$
f_{b,k}
=
\begin{cases}
f_{i,k}
+
\left(\boldsymbol{x}_b-\boldsymbol{x}_i\right)\cdot\nabla f_{i,k},
&
v_{n,k}\ge 0,
\\
\rho_w g_{w,k},
&
v_{n,k}<0.
\end{cases}
$$

壁面单位密度 Maxwell 分布记为

$$
g_{w,k}
=
g_w
\left(
\boldsymbol{v}_k;
\boldsymbol{U}_w,
T_w
\right).
$$

若采用无量纲气体常数为 1 的写法，壁面 Maxwell 分布为

$$
g_w
\left(
\boldsymbol{v};
\boldsymbol{U}_w,
T_w
\right)
=
\frac{1}{\left(2\pi T_w\right)^{3/2}}
\exp
\left[
-
\frac{
\left|
\boldsymbol{v}-\boldsymbol{U}_w
\right|^2
}{
2T_w
}
\right].
$$

以三维单原子气体为例，因此这里已经取平动速度维数为 3。

壁面漫反射密度系数由质量无穿透条件确定：入射壁面的质量通量与壁面反射回来的质量通量相加必须为零，

$$
\sum_{v_{n,k}>0}
v_{n,k}
f_{b,k}^{+}
\Delta\Xi_k
+
\sum_{v_{n,k}<0}
v_{n,k}
f_{b,k}^{-}
\Delta\Xi_k
=
0.
$$

将反射分布函数写成

$$
f_{b,k}^{-}
=
\rho_w g_{w,k},
$$

则壁面反射密度为

$$
\rho_w
=
-
\frac{
\sum_{v_{n,k}>0}
v_{n,k}
f_{b,k}^{+}
\Delta\Xi_k
}{
\sum_{v_{n,k}<0}
v_{n,k}
g_{w,k}
\Delta\Xi_k
}.
$$

因此，介观漫反射边界的计算步骤是：

1. 根据内场分布函数外推得到入射分布；
2. 对入射分布计算入射质量通量；
3. 对单位密度壁面 Maxwell 分布计算反射方向的单位质量通量；
4. 用质量无穿透条件求壁面反射密度；
5. 得到完整的边界分布函数。

## 2. 由介观分布函数积分得到宏观半侧通量

介观边界分布函数确定以后，宏观边界通量由分布函数的速度矩给出。设碰撞不变量为

$$
\boldsymbol{\psi}
=
\left[
1,\,
v_x,\,
v_y,\,
v_z,\,
\frac{1}{2}
\left(
v_x^2+v_y^2+v_z^2
\right)
\right]^T .
$$

则入射半空间宏观通量为

$$
\boldsymbol{F}_b^+
=
\int_{\Xi^+}
v_n
f_b^+
\boldsymbol{\psi}
d\Xi .
$$

反射半空间宏观通量为

$$
\boldsymbol{F}_b^-
=
\int_{\Xi^-}
v_n
f_b^-
\boldsymbol{\psi}
d\Xi .
$$

总边界通量为

$$
\boldsymbol{F}_b
=
\boldsymbol{F}_b^+
+
\boldsymbol{F}_b^- .
$$

在离散速度空间中，上述积分变成求和：

$$
\boldsymbol{F}_b^+
=
\sum_{v_{n,k}>0}
v_{n,k}
f_{b,k}^{+}
\boldsymbol{\psi}_k
\Delta\Xi_k,
$$

$$
\boldsymbol{F}_b^-
=
\sum_{v_{n,k}<0}
v_{n,k}
f_{b,k}^{-}
\boldsymbol{\psi}_k
\Delta\Xi_k.
$$

对于三维单原子气体，边界宏观通量直接写成质量通量、动量通量和总能量通量：

$$
F_{\rho}
=
\sum_k
v_{n,k}
f_{b,k}
\Delta\Xi_k,
$$

$$
\boldsymbol{F}_{\rho\boldsymbol{u}}
=
\sum_k
v_{n,k}
\boldsymbol{v}_k
f_{b,k}
\Delta\Xi_k.
$$

$$
F_E
=
\sum_k
v_{n,k}
\frac{1}{2}
\left|
\boldsymbol{v}_k
\right|^2
f_{b,k}
\Delta\Xi_k .
$$

如果需要显式区分半空间，则只需把上式的求和范围分别限制在入射半空间和反射半空间，再把两部分相加：

$$
F_{\rho}
=
F_{\rho}^{+}
+
F_{\rho}^{-},
\qquad
\boldsymbol{F}_{\rho\boldsymbol{u}}
=
\boldsymbol{F}_{\rho\boldsymbol{u}}^{+}
+
\boldsymbol{F}_{\rho\boldsymbol{u}}^{-},
\qquad
F_E
=
F_E^{+}
+
F_E^{-}.
$$

因此，由 CIS 算出的边界通量可以记为

$$
\boldsymbol{F}_{DVM}^{0}
=
\boldsymbol{F}_{DVM}^{+,0}
+
\boldsymbol{F}_{DVM}^{-,0}.
$$

这个通量包含两类信息：第一类是局部 Maxwell 平衡态对应的主通量，第二类是 DVM 分布函数中非平衡部分贡献的高阶通量。如果在宏观 NS 内迭代中始终固定使用这个 DVM 通量，那么边界通量不会随着宏观变量的快速更新而变化，近壁区域的信息传播会滞后。

## 3. Maxwell 边界条件对应的解析半空间通量

为了让宏观边界通量在 NS 内迭代中实时更新，需要把平衡态部分写成 Maxwell 分布的半空间解析积分。

在边界面上建立三维局部正交坐标系：

$$
\boldsymbol{n},
\qquad
\boldsymbol{t}_1,
\qquad
\boldsymbol{t}_2.
$$

其中

$$
\boldsymbol{n}\cdot\boldsymbol{t}_1
=
0,
\qquad
\boldsymbol{n}\cdot\boldsymbol{t}_2
=
0,
\qquad
\boldsymbol{t}_1\cdot\boldsymbol{t}_2
=
0.
$$

宏观速度在局部坐标中的分量为

$$
u_n
=
\boldsymbol{u}\cdot\boldsymbol{n},
\qquad
u_{t1}
=
\boldsymbol{u}\cdot\boldsymbol{t}_1,
\qquad
u_{t2}
=
\boldsymbol{u}\cdot\boldsymbol{t}_2.
$$

分子速度在局部坐标中的分量为

$$
v_n
=
\boldsymbol{v}\cdot\boldsymbol{n},
\qquad
v_{t1}
=
\boldsymbol{v}\cdot\boldsymbol{t}_1,
\qquad
v_{t2}
=
\boldsymbol{v}\cdot\boldsymbol{t}_2.
$$

局部 Maxwell 分布可以分解成一个法向速度分布和两个切向速度分布的乘积：

$$
M
\left(
\boldsymbol{v};
W
\right)
=
\rho
M_n
\left(
v_n;
u_n,
T
\right)
M_{t1}
\left(
v_{t1};
u_{t1},
T
\right)
M_{t2}
\left(
v_{t2};
u_{t2},
T
\right).
$$

因此，半空间通量的关键只剩下法向一维半矩。定义正半空间法向半矩：

$$
A_p^+
\left(
u_n,T
\right)
=
\int_0^{+\infty}
v_n^p
M_n
\left(
v_n;
u_n,
T
\right)
dv_n .
$$

定义负半空间法向半矩：

$$
A_p^-
\left(
u_n,T
\right)
=
\int_{-\infty}^{0}
v_n^p
M_n
\left(
v_n;
u_n,
T
\right)
dv_n .
$$

令

$$
\lambda
=
\frac{1}{2T},
\qquad
\chi
=
\sqrt{\lambda}\,u_n,
\qquad
E
=
\exp
\left(
-\lambda u_n^2
\right).
$$

零阶半矩为

$$
A_0^+
=
\frac{1}{2}
\left[
1+\operatorname{erf}
\left(
\chi
\right)
\right],
$$

$$
A_0^-
=
\frac{1}{2}
\operatorname{erfc}
\left(
\chi
\right).
$$

一阶半矩为

$$
A_1^+
=
u_n A_0^+
+
\frac{E}{2\sqrt{\pi\lambda}},
$$

$$
A_1^-
=
u_n A_0^-
-
\frac{E}{2\sqrt{\pi\lambda}}.
$$

更高阶半矩由递推公式得到：

$$
A_{p+2}^{\pm}
=
u_n A_{p+1}^{\pm}
+
\left(
p+1
\right)
T A_p^{\pm},
\qquad
p\ge 0.
$$

因此，宏观通量需要的几个半矩

$$
A_1^{\pm},
\qquad
A_2^{\pm},
\qquad
A_3^{\pm}
$$

不需要离散速度积分，而是通过误差函数和递推公式解析得到。

对于三维单原子模型，Maxwell 半空间通量在局部坐标中为

$$
F_{\rho}^{\pm}
=
\rho A_1^{\pm},
$$

$$
F_{\rho u_n}^{\pm}
=
\rho A_2^{\pm},
$$

$$
F_{\rho u_{t1}}^{\pm}
=
\rho u_{t1} A_1^{\pm},
$$

$$
F_{\rho u_{t2}}^{\pm}
=
\rho u_{t2} A_1^{\pm},
$$

$$
F_E^{\pm}
=
\frac{1}{2}
\rho
\left[
A_3^{\pm}
+
\left(
u_{t1}^2
+
u_{t2}^2
+
2T
\right)
A_1^{\pm}
\right].
$$

最后需要把局部坐标中的动量通量投影回全局坐标。三维情况下：

$$
\boldsymbol{F}_{\rho\boldsymbol{u}}^{\pm}
=
F_{\rho u_n}^{\pm}
\boldsymbol{n}
+
F_{\rho u_{t1}}^{\pm}
\boldsymbol{t}_1
+
F_{\rho u_{t2}}^{\pm}
\boldsymbol{t}_2.
$$

因此，Maxwell 平衡态的半空间解析通量可以写成

$$
\boldsymbol{F}_{M}^{\pm}
\left(
W
\right)
=
\left[
F_{\rho}^{\pm},\,
\left(
\boldsymbol{F}_{\rho\boldsymbol{u}}^{\pm}
\right)^T,\,
F_E^{\pm}
\right]^T .
$$

对于完全漫反射壁面，入射壁面的 Maxwell 通量由内场宏观状态决定，反射回来的 Maxwell 通量由壁面状态决定：

$$
\boldsymbol{F}_{M,b}
=
\boldsymbol{F}_{M}^{+}
\left(
W_i
\right)
+
\boldsymbol{F}_{M}^{-}
\left(
\rho_w,
\boldsymbol{U}_w,
T_w
\right).
$$

壁面反射密度仍然由质量无穿透条件决定：

$$
F_{\rho,M}^{+}
\left(
W_i
\right)
+
F_{\rho,M}^{-}
\left(
\rho_w,
\boldsymbol{U}_w,
T_w
\right)
=
0.
$$

若壁面 Maxwell 分布写成单位密度形式，则

$$
\rho_w
=
-
\frac{
F_{\rho,M}^{+}
\left(
W_i
\right)
}{
A_{1,w}^{-}
}.
$$

其中

$$
A_{1,w}^{-}
=
A_1^-
\left(
u_{n,w},
T_w
\right).
$$

对于静止壁面，

$$
u_{n,w}=0,
$$

因此

$$
A_{1,w}^{-}
=
-
\sqrt{
\frac{T_w}{2\pi}
}.
$$

这一部分就是宏观 Maxwell 边界条件的核心。它不需要访问介观分布函数，也不需要在每次 NS 内迭代中做离散速度积分，只需要由当前宏观边界状态计算 Maxwell 半空间解析矩。

## 4. 边界高阶项与迭代法解通量

GSIS 的宏观边界通量不是简单固定为 DVM 通量，也不是完全丢掉介观非平衡信息，而是采用如下结构：

$$
\text{当前 Maxwell 平衡态通量}
+
\text{冻结的 DVM 高阶修正}.
$$

外层 DVM 给出的初始边界通量记为

$$
\boldsymbol{F}_{DVM}^{0}.
$$

用同一个初始宏观边界状态构造 Maxwell 解析通量：

$$
\boldsymbol{F}_{M}^{0}
=
\boldsymbol{F}_{M}
\left(
W^0
\right).
$$

边界高阶项定义为

$$
\boldsymbol{F}_{HoT}
=
\boldsymbol{F}_{DVM}^{0}
-
\boldsymbol{F}_{M}^{0}.
$$

在宏观第

$$
m
$$

次内迭代中，边界宏观状态更新为

$$
W^m.
$$

此时边界通量取为

$$
\boldsymbol{F}_b^m
=
\boldsymbol{F}_{M}
\left(
W^m
\right)
+
\boldsymbol{F}_{HoT}.
$$

代入高阶项定义，可以得到

$$
\boldsymbol{F}_b^m
=
\boldsymbol{F}_{M}
\left(
W^m
\right)
+
\left[
\boldsymbol{F}_{DVM}^{0}
-
\boldsymbol{F}_{M}
\left(
W^0
\right)
\right].
$$

等价地，也可以写成增量修正形式：

$$
\boldsymbol{F}_b^m
=
\boldsymbol{F}_{DVM}^{0}
+
\left[
\boldsymbol{F}_{M}
\left(
W^m
\right)
-
\boldsymbol{F}_{M}
\left(
W^0
\right)
\right].
$$

如果显式分开入射和反射两部分，则有

$$
\boldsymbol{F}_b^{+,m}
=
\boldsymbol{F}_{DVM}^{+,0}
+
\left[
\boldsymbol{F}_{M}^{+}
\left(
W_i^m
\right)
-
\boldsymbol{F}_{M}^{+}
\left(
W_i^0
\right)
\right],
$$

$$
\boldsymbol{F}_b^{-,m}
=
\boldsymbol{F}_{DVM}^{-,0}
+
\left[
\boldsymbol{F}_{M}^{-}
\left(
W_w^m
\right)
-
\boldsymbol{F}_{M}^{-}
\left(
W_w^0
\right)
\right].
$$

总通量为

$$
\boldsymbol{F}_b^m
=
\boldsymbol{F}_b^{+,m}
+
\boldsymbol{F}_b^{-,m}.
$$

对于漫反射壁面，反射通量不能随意给定。每一次宏观内迭代都应先用当前入射 Maxwell 质量通量确定壁面反射密度：

$$
\rho_w^m
=
-
\frac{
F_{\rho,M}^{+}
\left(
W_i^m
\right)
}{
A_{1,w}^{-}
\left(
T_w
\right)
}.
$$

然后再计算当前反射 Maxwell 通量：

$$
\boldsymbol{F}_{M}^{-}
\left(
W_w^m
\right)
=
\boldsymbol{F}_{M}^{-}
\left(
\rho_w^m,
\boldsymbol{U}_w,
T_w
\right).
$$

因此，漫反射壁面的 GSIS 边界通量迭代过程可以写成以下步骤：

1. 外层 DVM 根据介观漫反射边界条件计算边界分布函数；
2. 对边界分布函数进行离散速度积分，得到初始 DVM 边界通量；
3. 用初始宏观边界状态计算 Maxwell 解析通量；
4. 计算并冻结边界高阶项；
5. 在宏观 NS 内迭代中，由当前宏观解重构边界入射状态；
6. 用 Maxwell 半空间解析公式计算当前入射通量；
7. 由质量无穿透条件计算当前壁面反射密度；
8. 用当前壁面反射密度计算当前反射通量；
9. 将当前 Maxwell 通量与冻结高阶项相加，得到最终宏观边界通量；
10. 用该通量进入宏观有限体积残差。

用公式概括，上述步骤就是

$$
\boldsymbol{F}_{HoT}
=
\boldsymbol{F}_{DVM}^{0}
-
\boldsymbol{F}_{M}
\left(
W^0
\right),
$$

$$
\boldsymbol{F}_b^m
=
\boldsymbol{F}_{M}
\left(
W^m
\right)
+
\boldsymbol{F}_{HoT}.
$$

迭代法解边界通量，并不是每次宏观内迭代都重新进行 DVM 速度空间积分，而是每次宏观内迭代只更新 Maxwell 平衡态的解析半空间通量。DVM 与 Maxwell 之间的差值作为高阶项在当前宏观内迭代中冻结。等外层 DVM 分布函数再次更新后，边界高阶项才重新计算。

在线性化版本直接装配 Ax=b 时，可以先只用宏观 Maxwell 边界条件迭代右端项 b；随后再加入 DVM 提供的边界高阶项和内场高阶项，作为对比测试。
