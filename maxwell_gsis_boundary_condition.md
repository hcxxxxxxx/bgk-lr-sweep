# NSW 的介观、宏观与 GSIS 边界条件

本文配合 [求解流程与子程序指南](SOLVER_GUIDE.md)，对应当前 monoCIS.m、monoNS.m、nsw_init.m 和 nsw_closure.m。原参考目录的同名文档讨论 Maxwell 漫反射壁面；本例没有漫反射壁面，不能照搬壁温和反射密度公式。

## 1. 先统一入射、出射和法向

这里的入射指进入计算域，出射指离开计算域。边界法向 n 指向域外，vn=v·n：vn<0 为入射，vn>0 为出射。某些壁面文献把朝向壁面称为入射，其叫法相反，判断时应看 vn 的符号。

| 边界 | 外法向 | 从外部进入计算域的速度 | 当前条件 |
|---|---|---|---|
| 左侧 | (-1,0) | vx>0 | 给定上游 Maxwell 分布 |
| 右侧 | (1,0) | vx<0 | 给定下游 Maxwell 分布 |
| 上侧 | (0,1) | vy<0 | 内部分布按 vy 反号镜面反射 |
| 下侧 | (0,-1) | vy>0 | 内部分布按 vy 反号镜面反射 |

例如左侧正向运动的分子 vx>0，此时 vn=-vx<0，应取上游给定值，而不是左侧内部出射值。

nsw_init 按几何位置设置 faces_tag：0 内部，1 上下对称，2 左远场，3 右远场。标签不是直接照搬 Fluent 区域编号。当前反射索引只针对水平上下边界。

## 2. 远场状态从哪里来

上游 rho_u=T_u=1，u_u=Ma sqrt(gamma)，gamma=5/3。令 M=Ma：

```text
rho_d/rho_u = (gamma+1)M² / [(gamma-1)M²+2]
p_d/p_u = [2 gamma M²-(gamma-1)]/(gamma+1)
T_d/T_u = (p_d/p_u)/(rho_d/rho_u)
u_d/u_u = rho_u/rho_d
```

Ma=3 时，下游 rho=3、T=11/3、ux=sqrt(15)/3，上游 ux=sqrt(15)，两端 uy=0。这些状态保存在 gas.left 和 gas.right。它们给定远场入射信息，不强行覆盖边界相邻物理单元的全部变量。

## 3. CIS 的远场通量

monoCIS 初始化 fBoundary、hBoundary，采用 fM=rho/(2 pi T) exp[-|v-u|²/(2T)]、hM=T fM。dvm_boundary 拼接内部与虚拟分布；dvm_convection 按速度符号选边：

```text
面分布：vn≥0 取内部 f；vn<0 取远场 fM
面输运：max(vn,0) f_internal + min(vn,0) fM_far
```

h 使用相同规则。面长在 dvm_convection 的 pos、neg 中已经包含。虽然虚拟行的所有速度都填了 Maxwell 值，只有入射半空间进入实际通量；这不是完整分布 Dirichlet 边界。

## 4. 上下对称与漫反射的区别

dvm_boundary 对上下边界采用 fghost(vx,vy)=finside(vx,-vy)，h 同理。reflect 通过翻转速度网格的 vy 索引构造，因此不需要插值。边界没有给定壁温，也没有求解反射 Maxwell 密度。

宏观 ns_boundary 在这些面使用 F=[0,p nx,p ny,0]，即无质量穿透和无能量通量的压力贡献。当前没有显式加入对称面偏应力和边界 HoT；一维条带解的上下贡献成对抵消，不能据此宣称这是通用黏性曲壁条件。

漫反射则需要壁面温度、壁面速度，并通常用零法向质量流确定反射密度。它描述另一类物理边界，仅在真正有壁面时使用。

## 5. NS 的半 Maxwell 远场通量

nsw_closure 调用局部 half 两次：内部当前 q 的出射半空间，加给定远场状态的入射半空间。宏观内迭代改变内部 q 时，出射通量必须重新计算。

令 un=u·n，ut=-ux ny+uy nx，s=+1 表示出射、s=-1 表示入射。half 的解析矩为：

```text
m0 = 0.5 erfc[-s un/sqrt(2T)]
m1 = un m0 + s sqrt(T/(2 pi)) exp[-un²/(2T)]
m2 = un m1 + T m0
m3 = un m2 + 2T m1
Fhalf = rho [m1,
             m2 nx-ut m1 ny,
             m2 ny+ut m1 nx,
             0.5(m3+(ut²+2T)m1)]
FM(qinside,qfar) = Fhalf(qinside,+1)+Fhalf(qfar,-1)
```

四列依次是质量、x 动量、y 动量和总能量通量。m1 的入射值可以为负，这是按外法向计的有符号通量，不是错误。能量中的 2T 来自两个切向热运动自由度。

这里的 F 是每单位面长通量；ns_residual 最后才乘面长、按 scatter 求和并除以单元面积。纯 NS 令 hot.boundary=0，直接使用 FM。

## 6. GSIS 的冻结边界差值

介观步完成后，dvm_hot 从真实迎风分布求矩得到 Fkin。使用相同介观矩 qkin 计算解析 FM，定义：

```text
Dfar = Fkin - FM(qkin,qfar)
内迭代边界通量 = FM(qNS_current,qfar) + Dfar
```

本轮全部 NS 内迭代中 Dfar 不变，但 FM 随 qNS_current 更新。下一次介观步完成后重新计算 Dfar。这样既保留宏观边界对当前宏观状态的响应，也保留介观边界的非平衡及离散差值信息。

仅在参考状态 qNS=qkin 时，有 FM(qkin,qfar)+Dfar=Fkin；不能说整个内迭代始终严格等于旧介观通量。本版宏观状态持续保留，因此每次调用 NS 的初始状态也不一定等于当前 qkin。

Dfar 已是完整边界通量差，不应再在同一个远场面叠加一次格心应力/热流 HoT。hot.cell 在 ns_viscous 的内部面使用，hot.boundary 在 ns_boundary 的远场面使用。对称面不在 blks.far 中，不参与这个恢复恒等式。

## 7. 按代码核对的顺序

1. nsw_init：检查法向向外、标签和 gas.left/right。
2. monoCIS 初始化：检查 fBoundary、hBoundary 与速度反射索引。
3. dvm_boundary → dvm_convection：检查虚拟值及入射/出射选择，面长只乘一次。
4. nsw_closure → half：检查解析边界通量与法向符号。
5. dvm_hot：检查 Fkin 和 FM 使用相同的介观参考状态与相同通量单位。
6. ns_boundary：检查冻结差值只加一次，当前 Maxwell 部分持续更新。

已有固定 HoT 测试中远场通量恢复差为 0，但从收敛 CIS 推进 100 步 NS 后宏观场仍变化约 1.57%。这说明边界差值的代数恢复成立，不代表内部 NS+HoT 离散与 CIS 完全一致；不能只凭边界检查通过就判定整个耦合固定点正确。
