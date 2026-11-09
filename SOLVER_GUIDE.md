# 从子程序读懂 CIS、NS 与 GSIS

本文对应当前目录的 MATLAB 源码，不是其他 DG 或混合气体版本的算法说明。算例是单组分单原子 Ma=3 正激波，400 个一阶有限体积单元，采用 Fluent 非结构拓扑读取；目前验证的几何是正交四边形条带。

## 1. 先运行哪一个文件

在 MATLAB 中进入本目录：

```matlab
monoCIS                 % 纯介观求解，默认模式
monoNS                  % 纯宏观求解，无输入参数
demo_method = 'GSIS';
monoCIS                 % 介观主程序调用同一个 monoNS
```

建议先分别读懂并验证纯 CIS 和纯 NS，再看耦合分支。直接在 monoCIS 中按 F5 默认运行 CIS，不会自动选择 GSIS。demo_method 在读入后被清除，因此下次直接运行又是 CIS。

主程序负责安排顺序；文件下方的局部函数负责一个具体步骤。局部函数不需要逐个按 F5。MATLAB 函数有独立工作区：输入是括号中的参数，输出是等号左边的变量；主程序接收输出后才更新自己的状态。

## 2. 变量与数组：先分清物理空间和速度空间

| 变量 | 大小或内容 | 含义 |
|---|---|---|
| N | 400 | 物理单元数量，不是速度点数 |
| vx、vy | 1 × 7680 | 96 × 80 速度网格展开后的坐标 |
| dv | 标量 | 速度求积权重 Δvx Δvy |
| f、h | N × 7680 | 每一行对应一个物理单元，每一列对应一个离散速度 |
| qf、macro.qf | N × 4 | [rho, ux, uy, T]；前者来自介观矩，后者是宏观求解器状态 |
| sigmaq | N × 5 | [sigma_xx, sigma_xy, sigma_yy, qx, qy]，由分布函数求矩 |
| U | N × 4 | [rho, rho ux, rho uy, rho E]，宏观守恒变量 |
| hot.cell | N × 5 | 介观应力/热流减去同一介观状态的 NS 闭合 |
| hot.boundary | 左右远场面数 × 4 | 远场介观通量减去同状态半 Maxwell 通量 |
| cells.A、cells.xy | 面积、格心坐标 | 有限体积几何量 |
| edges.enormal、elength | 面法向、面长 | 面法向从左单元指向右单元，边界处指向域外 |

f 是三维分布对第三速度分量积分后的约化分布；h 是第三速度平方乘分布再积分后的约化分布。h 不是总能量。平衡时 h=T f，因此计算温度和热流时不能漏掉 h。

blks.faces 的第 3、4 列给出面左右单元，边界右单元为 0。faces_ghost 将边界的 0 换成虚拟单元行号。scatter 是单元到面的有向求和矩阵：一个内部面的通量对左单元记正、对右单元记负，所以内部通量在全域相加时抵消。

代码中的 .*、./、.^ 是逐元素运算，* 是矩阵乘法。例如 f.*vx 在每一行乘速度，而 scatter*flux 是对每个单元的各个面求和。sum(f,2) 沿速度方向求矩；sum(q.^2,1) 沿物理单元方向计算每个宏观变量的范数。

## 3. 初始化：nsw_init

输入：无。输出：blks、cells、edges、gas、macro。

1. 调用 importFluentMesh2 读取 mesh/nsw400.cas，保留面与单元的连接关系。
2. 把坐标缩放到 x=[-60,60]、y=[0,1]，计算格心、面积、面长、法向和法向格心距离。
3. 建立边界虚拟单元、scatter、单元邻接表和 Jacobian 着色分组。
4. 设置 gamma=5/3、Pr=2/3、omega=0.74、delta=sqrt(pi/2)，由 Rankine-Hugoniot 关系给出上、下游状态。
5. 左半域取上游，右半域取下游，构成不连续初态；宏观 CFL 为 5。

这里既没有读取收敛解，也没有从其他精度的解初始化。着色仅用于减少 Jacobian 计算次数，不改变物理网格。

## 4. 纯 CIS：一个外迭代做什么

```text
保存旧宏观量
  → dvm_collision：构造冻结碰撞目标和频率
  → dvm_boundary：填充边界虚拟分布
  → dvm_convection：计算每个面的迎风通量
  → 主程序：装配碰撞减输运的残差
  → dvm_implicit：更新 f、h
  → dvm_moments：重新计算宏观量、应力和热流
  → 更新误差、显示、判断停止
```

### dvm_collision：碰撞目标

输入 q、heat、速度网格和 gas；输出 G、H、nu。它不直接修改 f、h。令 c=v-u、c²=cx²+cy²：

```text
E = rho/(2 pi T) exp[-c²/(2T)]
A = (1-Pr)/(5 rho T²) (qx cx + qy cy)
G = E [1 + A(c²/T-4)]
H = T E [1 + A(c²/T-2)]
nu = delta rho T^(1-omega)
```

E 是约化 Maxwell 分布，G、H 是 Shakhov 碰撞目标。纯 CIS 使用当前介观矩提供 q 和 heat。Shakhov 目标不保证速度空间逐点非负，不能把目标函数的负值与负温度混为一谈。

### dvm_boundary：边界数据

输入内部 f、h，远场平衡分布 fb、hb 和反射索引；输出拼接后的 fext、hext。左右边界的虚拟行采用给定远场 Maxwell 分布；上下边界复制内部状态中 vy 反号对应的分布。

填充虚拟行不等于把边界所有速度都固定。真正选择内部或虚拟值的是下一步迎风通量。详细符号见 [边界条件说明](maxwell_gsis_boundary_condition.md)。

### dvm_convection：迎风输运

输入左右单元索引以及 pos=max(v·n 面长,0)、neg=min(v·n 面长,0)。输出 fluxf、fluxh：

```text
fluxf = pos f_left + neg f_right
fluxh = pos h_left + neg h_right
```

这些通量已经乘了面长，之后不能再乘一次。主程序装配 rf=面积 nu (G-f)-scatter fluxf，rh 同理。这是单元积分残差。

### dvm_implicit：冻结碰撞目标的隐式更新

输入当前分布、残差、碰撞频率、出流系数与邻接关系；输出更新后的 f、h。对角项为 D=面积(nu+1e-5)+outflow，其中 1e-5 是伪时间步 1e5 的倒数，不是误差阈值。

先按单元编号正向扫描，再反向扫描；每次用邻居已有的 df、dh 修正本单元增量，最后 f←f+df、h←h+dh。这是当前介观的前后向 LU-SGS 更新，不是精确求解整个非线性系统。j>0 表示邻居是物理单元，不能删成无条件访问。

### dvm_moments：从分布求宏观矩

输入 f、h 和速度求积；输出 q 与 sq。求积符号 Σ 包含 dv：

```text
rho = Σ f                 u = Σ v f / rho
T = Σ[(cx²+cy²)f+h] / (3 rho)
sigma_xx = Σ cx² f-rho T   sigma_xy = Σ cx cy f
sigma_yy = Σ cy² f-rho T
qx = 0.5 Σ cx[(cx²+cy²)f+h]，qy 同理
```

这里的 sigma 是压力张量偏离平衡压力的部分，NS 闭合采用负速度梯度的符号。宏观黏性通量必须沿用该约定。

### dvm_hot：供 GSIS 使用的高阶项

纯 CIS 不调用此函数。输入更新后的介观分布及其矩；先调用 nsw_closure，在同一个介观 q 上计算 NS 应力/热流与解析边界通量，然后构造 hot.cell=sigmaq-closure、hot.boundary=Fkin-Fns。

hot.cell 是格心信息，内部面采用两侧平均。hot.boundary 只对应左右远场面，是每单位面长的通量差，不等同于格心 HoT。二者在一次宏观内迭代期间冻结。

### dvm_postplot：显示，不参与求解

输入窗口、格心坐标、宏观量和历史误差；按 x 排序画 rho、ux、T 和双对数更新误差。排序仅影响显示，不重排求解数组。剖面显示范围 [-15,15] 不是计算域范围。关闭窗口只停止该窗口的更新，不终止计算。

## 5. 纯 NS：同一个宏观求解器独立运行

monoNS 无参数时执行初始化，令两个 hot 数组为零，最多迭代 4000 步，更新误差小于 1e-8 停止。由 GSIS 带参数调用时，不重新初始化；从传入的 macro.qf 连续推进指定的 nstep 步。

```text
ns_conservative：q 转 U
循环：ns_residual → ns_timestep → ns_implicit → ns_primitive → 更新误差
结束：nsw_closure → 返回 macro 和 info
```

### ns_conservative 与 ns_primitive：变量转换

前者把 q=[rho,ux,uy,T] 转成 U，使用 rho E=rho[1.5T+(ux²+uy²)/2]；后者反算速度和温度。宏观守恒方程更新 U，物性与边界计算使用 q。二者只是代数转换，不是求矩，也不做插值。

### ns_residual：组织面通量与单元残差

输入 U、几何和冻结 hot；依次调用 ns_primitive、ns_convection、ns_viscous、ns_boundary，最后返回 R=-(scatter*(F 面长))/面积。

与介观不同，这里的 F 在装配前没有乘面长；R 已除以面积，是单位体积残差。比较不同模块残差时必须先统一量纲与归一化。

### ns_convection：内部面 Rusanov 通量

输入 U、q、面几何和 gas；计算左右 Euler 通量，再取 F=0.5[FL+FR-a(UR-UL)]，a 为两侧 |u·n|+sqrt(gamma T) 的最大值。输出面通量数组，边界面暂留给 ns_boundary。此 FV 版本不是 Roe 通量。

### ns_viscous：NS 闭合与 HoT 的面贡献

输入对流通量 F、q 和 hot；用两侧平均作为面状态，用 (qR-qL)/法向格心距离求法向导数。mu=T^omega/delta，kappa=2.5mu/Pr；黏性应力为负对称偏梯度，热流为 -kappa grad T。

应力和热流加上两侧 hot.cell 的平均值，再把 [0,sigma·n,u·(sigma·n)+qheat·n] 加到面通量。内迭代中，梯度、mu、速度随宏观状态更新，只有 HoT 保持不变。正交两点梯度没有一般非正交修正。

### ns_boundary：覆盖边界面通量

左右远场：调用 nsw_closure 得到当前宏观状态的半 Maxwell 通量，加上 hot.boundary。上下对称面：使用 [0,p nx,p ny,0]。它不是漫反射壁面，也不是给全部守恒量强制 Dirichlet 值。

### ns_timestep：宏观伪时间步

输入 q、CFL 和几何；返回每个单元的 dt=CFL 面积/各面波速乘面长之和。当前实现的面波速取面左单元状态，与 Rusanov 的两侧最大波速不是同一段公式。CFL 控制推进幅度，不是停止误差。

### ns_implicit：数值 Jacobian 与稀疏线性求解

输入 U、R、dt、几何和冻结 hot；用扰动 epsU=1e-7 max(|U|,1) 计算 J=∂R/∂U。邻接图着色将残差影响范围不相交的变量一起扰动，以减少 ns_residual 调用次数，再通过邻接关系提取各列。

求解 [diag(1/dt)-J]ΔU=R，并返回 U+ΔU。MATLAB 反斜杠完成稀疏线性求解；这一步不是宏观 LU-SGS，也不是 GMRES。计算 Jacobian 时 hot 冻结，但 Maxwell 边界和 NS 物性跟随扰动状态变化。

### ns_postplot：纯 NS 的监视窗口

功能与 dvm_postplot 相同，仅独立运行时调用。GSIS 调用 monoNS 的每一次内迭代不会新建窗口，避免大量窗口和绘图开销。

## 6. 共用闭合 nsw_closure 与局部 half

nsw_closure 输入 q、几何和 gas，输出格心 NS 应力/热流，以及左右远场半 Maxwell 通量。格心梯度用面平均值的 Gauss 求和；远场梯度虚拟值取给定状态，对称面虚拟值复制内部状态。该格心闭合供 dvm_hot 做差和记录 macro.sigmaq。

它与 ns_viscous 的面两点闭合用途不同，不能因都计算应力就随意合并。两种离散并不自动保证收敛 CIS 是 NS+HoT 的精确固定点。

局部 half 输入一个 Maxwell 状态、法向及半空间方向，输出质量、两方向动量和总能量的半空间通量。它利用 erfc 的解析半高斯矩，不遍历离散速度点。公式与边界差值含义见边界文档。

## 7. GSIS：外迭代与内迭代的交接

| 时刻 | 介观步骤 | 宏观步骤 |
|---|---|---|
| 第 1–9 步 | 与纯 CIS 相同 | 不推进 |
| 第 10 步 | 仍完成一整步 CIS，形成 HoT | 用介观 qf 初始化一次 macro.qf，开启耦合；尚不执行 100 步 NS |
| 第 11 步起 | 用上一轮 macro.qf 和总热流构造碰撞目标，输运并求矩 | 以新介观矩构造并冻结 HoT，从保留的宏观状态推进 100 步 |
| 本轮结束 | 保留新 f、h 与介观 qf | 保存 macro.qf，total=macro.sigmaq+hot.cell，供下一轮碰撞 |

所谓“总热流”是 total 的第 4、5 列。它包括宏观 NS 热流和介观高阶修正，而不是单独的 Fourier 热流。

本版有两个重要选择：宏观状态只在第 10 步初始化一次，不每轮重置为介观矩；不对分布执行 Maxwell 差值 VDF lift。宏观加速信息通过下一轮的碰撞目标反馈给介观。不能把它解释成“NS 直接覆盖 f”。

## 8. 误差、输出与验证界限

每个原始变量的更新误差为 sqrt[Σ(qnew-qold)²/max(Σqold²,N)]，取四个变量的最大值。纯 CIS 用介观更新误差，GSIS 用介观与宏观外层更新误差的最大值；NS 内迭代误差不是 GSIS 外层误差。分母下限避免零横向速度导致除零，并非密度/温度保正处理。

当前纯 CIS 阈值 1e-8，GSIS 阈值 1e-6、默认 NS100。已有 NS100 测试首次达到 1e-6 为总第 39 步，其中 10 步 CIS 启动；这是历史测试证据，不是所有机器运行时间的保证。

cis/gsis_history.dat 为步数、介观误差、宏观外层误差、最大误差、累计秒数。ns_history.dat 为步数、更新误差、该步更新前方程残差最大绝对值；info.equation 则在结束状态重新计算。更新很小不等于方程残差一定很小。

主要输出剖面来自介观 qf，macro.qf 另存。小更新误差不证明 CIS 与 GSIS 完全相同。已有固定 HoT 测试的远场通量恢复差为 0，但推进 100 步 NS 后宏观场仍变化约 1.57%，因此不能声称已证明严格离散固定点一致。

运行中会打印进度并显示监视图。纯 CIS 图在首步、每 100 步和结束更新，GSIS 另在每 10 步更新；NS 独立运行每 100 步更新。图不动不一定是程序停了，应同时看命令窗口。绘图会有额外开销，带监视功能的耗时不能直接视作之前无监视版本的性能基准。

## 9. 网格读取中的子程序

importFluentMesh2 负责把 Fluent 文本块转成拓扑数据，不负责物理边界方程。

| 子程序 | 阅读时抓住的作用 |
|---|---|
| readBlocks | 建立待填充的网格块结构 |
| nxtBlk | 按括号结构取出下一个文件块 |
| determineType | 判断块类型，决定按节点、面、单元等处理 |
| readMeshGenComment | 读取网格生成信息 |
| readComment | 注释读取辅助分支，不是求解步骤 |
| readHeader | 解析块头，包括十六进制编号 |
| readNodes | 读取节点坐标 |
| readCells | 读取单元块信息 |
| readFaces | 读取面节点及左右单元连接 |
| readZonesInfo | 当前保留的空兼容入口，不参与物理计算 |

初学时先读 nsw_init 如何使用读取结果，再回看解析器。不要把文件中的区域编号直接等同于求解器标签；本例的边界标签在初始化中按几何位置设置。
