# MATLAB 有限体积正激波 GSIS 算例

另有独立的 **1D2V Couette / Poiseuille** 文件 [`mono1D2V.m`](mono1D2V.m)：直接运行默认分别计算两个算例；用法见本文末尾。该文件不读取 Fluent 网格，也不调用原有激波求解器。

另有按用户提供 PDF 实现的独立线性化 BGK / Q1-DG 求解器 `monoBGKDG.m`，使用同一个 `mesh/nsw400.cas`。运行方式与模型说明见本文末尾的 **PDF BGK / Q1-DG solver**。原有 CIS、NS、GSIS 求解器保持原样。

## 阅读入口

- [测试算例设置](CASE_SETTINGS.md)：当前四组 1D2V 主算例的物理参数、网格、速度求积、参考解和复现命令，以及验证脚本和其他求解器的配置区别；已保存配置快照见 [CASE_SETTINGS.json](CASE_SETTINGS.json)。
- [求解流程与子程序指南](SOLVER_GUIDE.md)：变量和数组、纯 CIS、纯 NS、GSIS 调用顺序，以及每个局部子程序的输入输出、公式和作用。
- [NSW 边界条件](maxwell_gsis_boundary_condition.md)：远场入射/出射、上下镜面反射、解析半 Maxwell 通量和冻结边界差值；与原漫反射壁面文档区分。

建议先读变量表，再按纯 CIS → 纯 NS → GSIS 的顺序读代码。文档对应当前一阶 FV、保留宏观状态且不做 VDF lift 的版本。

## 代码结构

本版沿用 `D:/2025work/rot_gsis/nsw_release_rot` 的 Fluent 网格读取和 `blks / cells / edges` 组织方式，不再由一维数组生成物理网格。两个主求解文件为：

- `monoCIS.m`：初始化约化分布，计算碰撞、面迎风通量、LU-SGS 隐式更新以及介观矩；GSIS 模式中调用宏观求解器。
- `monoNS.m`：不带参数独立求纯 NS；带参数从给定宏观状态继续迭代，支持冻结 HoT。

三个辅助文件分别为 `importFluentMesh2.m`（沿用参考代码的网格读取，补充关闭文件句柄）、`nsw_init.m`（几何、邻接、初态）、`nsw_closure.m`（共用物理闭合与半空间边界通量）。没有引入类、通用求解器框架或 DG 数据结构。

两个主文件按阶段调用文件末尾的局部子程序，没有增加求解器文件：

- `monoCIS.m`：`dvm_collision`（碰撞）、`dvm_boundary`（边界分布）、`dvm_convection`（迎风通量）、`dvm_implicit`（LU-SGS）、`dvm_moments`（求矩）、`dvm_hot`（冻结高阶项和边界差值）。
- `monoNS.m`：`ns_residual`（残差装配）、`ns_convection`（对流）、`ns_viscous`（黏性与HoT）、`ns_boundary`（边界）、`ns_timestep`（时间步）、`ns_implicit`（隐式更新）、`ns_conservative/ns_primitive`（变量转换）。

移除了两个主求解循环中的逐步 `isreal/isfinite/正密度正温度` 断言；有限性与正宏观量由debug验证检查。网格合法性检查、实际边界分支、零横向速度所需的误差归一化，以及Jacobian扰动尺度保留，不作为冗余防御代码删除。数值公式、运算顺序和耦合顺序不变。

2026-10-08局部子程序重构回归：纯CIS前10步、纯NS完整2046步、GSIS完整39步的误差历史与重构前逐值一致（最大差为0）。GSIS最终误差9.793615772265e-7，终态分布有限、密度温度为正。未重新进行纯CIS的完整1691步计算。证据为 `debug/local_functions_regression.log`，隔离测试结果为 `debug/local_functions_regression/results/unstructured`，未覆盖原始历史。

介观与宏观均通过 `.cas` 中的左右单元、面法向、面长和单元面积装配残差，未将读入网格重新转换为 `i-1/i+1` 一维差分。宏观 Jacobian 着色也由邻接图生成。

**适用范围：** 当前验证的是正交四边形条带上的一维正激波，保留二维速度和二维宏观动量。粘性面梯度采用正交两点形式；尚未实现一般非正交网格修正，也未验证任意曲壁或三角网格。采用非结构网格数据格式不等于已经验证任意非结构几何。

## 运行顺序

在 MATLAB 中切换至本目录，先独立运行：

```matlab
monoCIS
monoNS
run('debug/check_unstructured.m')
```

确认纯求解器结果后再运行：

```matlab
demo_method = 'GSIS';
monoCIS
run('debug/check_unstructured.m')
```

`demo_method` 用后清除，下次直接执行 `monoCIS` 又是纯 CIS。宏观子求解是同一个 `monoNS`，不另设耦合专用宏观代码。验证脚本 `debug/check_unstructured.m` 用于检查输出和绘图。

直接在 monoCIS 按 F5 默认运行纯 CIS。当前版本在初始化前打印模式，随后显示 rho、ux、T 和双对数更新误差：纯 CIS 首步、每100步和结束刷新图，GSIS 另在每10步刷新；命令窗口每5步打印介观进度。纯 NS 首步、每100步和结束刷新。关掉图窗不会停止求解。监视功能的两步 CIS/NS 冒烟测试通过，见 `debug/monitor_smoke.log`；它不代替完整收敛验证，绘图也有额外耗时。

## 物理与数值设置

- 单组分单原子 Shakhov 模型：Ma=3、gamma=5/3、Pr=2/3、omega=0.74、delta=sqrt(pi/2)。上游 rho=T=1，下游按 Rankine-Hugoniot 关系设置。
- 网格来自参考目录 `nsw401.cas`，另存为 `mesh/nsw400.cas`。实际为 **400 个物理单元、1201 个面**；在初始化中将坐标缩放至 x=[-60,60]、y=[0,1]，保留原拓扑与原文件。
- 一阶 FV；vx=[-18,18] 的96点，vy=[-15,15] 的80点，等权矩形速度求积。
- 左右远场仅指定入射 Maxwell 分布，出射由内部解给定；上下采用镜面反射。宏观左右使用解析半 Maxwell 通量，上下使用无穿透压力通量。不是将全部边界分布强制设成 Dirichlet 值。
- 两种方法从域中心的上、下游不连续平衡态开始，不读取其他阶数或收敛解作为初值。
- 介观采用冻结碰撞目标的前后向 LU-SGS，伪时间步1e5；宏观采用邻接图着色数值 Jacobian 和稀疏隐式求解，CFL=5。
- GSIS 保留前10步 CIS；每次宏观内迭代100步，停止阈值1e-6。纯CIS和纯NS仍使用1e-8。

## 变量与耦合

`qf=[rho,ux,uy,T]`，`sigmaq=[sigma_xx,sigma_xy,sigma_yy,qx,qy]`。`f`、`h` 是第三速度方向积分后的约化分布；平衡态满足 h=T*f。

1. 前10步纯 CIS，第10步后用介观矩初始化一次 NS 状态。
2. 后续碰撞的密度、速度、温度、碰撞频率全部使用持续保存的 NS 状态。Shakhov 热流使用 NS 闭合加 HoT。
3. 推进一步介观，取新介观矩；在同一介观状态计算物理 NS 应力/热流，形成 `HoT = kinetic - NS(kinetic)`。
4. 左右边界另取实际离散速度通量与相同状态解析半 Maxwell 通量的差值。
5. 冻结上述 HoT 和边界差值，调用 `monoNS` 从上一轮宏观状态继续推进。
6. 不用介观量逐轮重置宏观解，也不执行 Maxwell 差值 VDF lift。宏观量及总热流仅反馈下一轮碰撞。

结果中的主要剖面来自介观矩 `qf`；`macro.qf` 另存，不混用。第一阶 CIS 与无 lift GSIS 的离散固定点不保证相同，应直接检查剖面差异，不能仅凭下降步数宣称相同解的加速。纯 NS 在此稀薄激波中也不应严格等同于动理学解。

## 数据与验证

新版结果在 `results/unstructured`，原生日志在 `debug/unstructured_*.log`。原一维数组版已归档于 `debug/structured_array_20261007`，其步数与图不作为新版验证证据。旧的 `debug/check_pure.m` 等仅属于旧版，不用于新版验收。

- `cis/gsis_profile.dat`：x,y,rho,ux,uy,T,sigma_xx,sigma_xy,sigma_yy,qx,qy。
- `cis/gsis_history.dat`：迭代步、介观更新误差、宏观外层更新误差、两者最大值、累计求解秒数。
- `ns_history.dat`：迭代步、宏观更新误差、该步更新前方程残差最大绝对值。

更新误差为各宏观变量的 `sqrt(sum((q_new-q_old).^2)/max(sum(q_old.^2),Ncell))` 的最大值。分母以无量纲单位范数为下限，避免零横向速度产生除零；它不是方程残差。密度和纵向速度等非零量仍按自身范数归一化。

检查顺序为：网格几何闭合和正面积 → 纯 CIS/纯 NS 有限值与格心正密度正温度 → 稳态方程及守恒核查 → GSIS 前10步一致性、耦合收敛和剖面对比。Shakhov 分布不保证逐点保正，不通过截断 VDF 制造通过。

## 本版实测结果

2026-10-08更新：默认改为NS100，GSIS在第39步（含10步CIS启动）首次达到1e-6。依据见 `results/ns_steps_20261008/README.md`。下表和 `results/unstructured` 中尚未覆盖的历史结果属于NS20配置，不是新默认值的运行记录；旧ZIP同样保留NS20。当前源代码未重新运行覆盖历史数据。

MATLAB R2024b，本机单次运行；以下为求解循环 wall time，不是标准性能基准，也不是 CPU time。

| 方法 | 达到停止阈值的步数 | 最终更新误差 | wall time |
|---|---:|---:|---:|
| CIS | 1691 | 9.9830e-9 | 723.03 s |
| NS | 2046 | 9.9608e-9 | 16.08 s |
| GSIS | 148，包含10步CIS | 8.6245e-9 | 95.10 s |

CIS 最终动理学方程最大绝对残差3.6885e-8，轴向质量/动量/能量通量最大相对变化2.6647e-6；NS 最终方程残差3.7250e-6。两种独立求解器先通过验证，随后才运行 GSIS。前10步 CIS 与 GSIS 的更新历史逐值完全一致。分布有限，格心密度与温度为正。

固定 HoT 检查中，边界通量恢复误差为0。从收敛 CIS 状态推进100步 NS 后，宏观场变化约1.57%，方程残差约1.00e-2；**未证明 CIS 离散解是此 NS+HoT 算子的严格固定点**。未添加内部全通量修正强制两者一致。

GSIS 与 CIS 未平移的 rho/ux/T 归一化 L2 差异分别为0.945%、0.765%、0.757%。应力和热流仍有更明显差异，见 `comparison.png`，不能仅以更新误差小或迭代少宣称物理精度等同。这些是当前一阶离散的结果，不是空间独立性结论。

详细记录见 `results/unstructured/validation.txt`。轻量 ZIP 提供源码、原始网格、验证脚本、DAT、PNG及日志，不包含大型分布 MAT；按运行顺序执行后会重新生成 MAT。旧的一维数组版不包含在新版 ZIP 中。

## PDF BGK / Q1-DG solver

`monoBGKDG.m` 对应 `2D2V_BGK_Pseudo2D_DG_Source_Iteration.pdf` 的式 (1)-(17)。它是独立的 **2D2V 线性化 BGK、Q1 迎风 DG、源迭代** 求解器，没有 Shakhov 热流修正，也不调用 `monoNS`。

在 MATLAB 中切换到本目录执行：

```matlab
result = monoBGKDG;
```

默认直接读取现有 `mesh/nsw400.cas` 的节点、面和单元连接，保留其原始坐标 **x=[0,100]、y=[0,0.25]**，共 400 个矩形单元、1201 个面。默认参数为 `delta=0.01`（即 `delta*Lx=1`）、`Tw=[0.95,1.05]`、平均完整分布密度 1、零体源项、零初始扰动，以及 16×16 Gauss-Hermite 速度点。PDF 没有规定这些算例数值；这些是此实现提供的默认值。与 `nsw_init.m` 不同，此求解器默认不将网格缩放为 120×1。

每个单元、每个速度有 4 个独立自由度，节点顺序为 **左下、右下、左上、右上**。几何和邻居由 `.cas` 文件恢复；输入须为单行、轴对齐矩形网格。求解器保留 DG 单元间的双侧迹，不将自由度合并为连续节点。

模型与 PDF 使用同一归一化：

```text
M0 = exp(-vx²-vy²)/pi,  f = M0*(1+h)
v·grad(h) + delta*h = delta*P(h) + G
P(h) = rho + 2*vx*ux + 2*vy*uy + (vx²+vy²-1)*T
```

代码使用 PDF 建议的 Gram 矩阵离散投影 `B*(B'*W*B)^(-1)*B'*W`，在机器精度内守恒四个碰撞不变量。Gauss-Hermite 的 `velocity.weight` 已含 `M0`，求矩时不应再乘一次 `M0`。这里的 `h` 是 PDF 中的分布扰动，与原 `monoCIS.m` 的第三速度方向约化能量分布 `h` 含义不同。

左右采用完整 Maxwell 漫反射：在各壁面两个 Gauss 点上，用上一轮出射的 `M0*(1+h)` 求反射密度，再组装入射 `rho_wall*Mw/M0-1`。反射密度随解更新，壁温给定；没有用固定远场密度代替零质量通量条件。上下采用同一位置、`vy` 反号的上一轮分布。壁面分子入射值和碰撞增益均滞后一轮，内部迎风迹使用本轮已更新的上游解。

对每个速度，组装由 4×4 单元块组成的稀疏块三角矩阵并只分解一次。正 `vx` 的上游块位于低编号单元，负 `vx` 的上游块位于高编号单元；每轮回代与 PDF 的逐单元速度扫掠等价。内部面的交叉迹矩阵显式连接本单元入流边与邻居出流边，未用本单元自面矩阵替代。`debug/check_bgk_dg.m` 通过独立求积和显式逐单元扫掠核对这一实现。

默认用完整分布总质量约束消除封闭域的质量自由度。零源项时，按 `f <- alpha*f`（即 `h <- alpha*h+alpha-1`）归一化到给定平均密度；这个变换保留零源项方程和漫反射齐次关系，包括两端壁温不等的情况。有非零兼容源项时，通过常数密度校正施加质量约束，最终仍须满足未经校正的原方程残差。非零总质量生产的体源项会被拒绝。

### 修改参数

```matlab
options = struct;
options.delta = 0.1;
options.wallTemperature = [0.98,1.02];
options.velocityOrder = [20,20];
options.plot = false;
result = monoBGKDG(options);

% 如需更改物理尺寸：仍使用原文件的节点和拓扑，只仿射缩放坐标。
% 不设置 domain 时始终保留 .cas 中的原始坐标。
result = monoBGKDG(struct('domain',[0,1,0,1],'delta',1));
```

| 参数 | 默认值 | 含义 |
|---|---|---|
| `meshFile` | 本目录 `mesh/nsw400.cas` | Fluent 文本网格 |
| `domain` | `[]` | 原始坐标；可选 `[xmin,xmax,ymin,ymax]` |
| `delta` | `0.01` | 常数碰撞参数，须为正 |
| `wallTemperature` | `[0.95,1.05]` | 左、右壁面温度 |
| `velocityOrder` | `[16,16]` | x、y 方向 Gauss-Hermite 点数，均至少 3；支持奇数阶及零速度 |
| `meanDensity` | `1` | 完整分布的体平均密度 |
| `source`, `initial` | `0` | 标量、节点×速度数组，或 `@(x,y,vx,vy)` 函数 |
| `maxIter` | `30000` | 最大源迭代次数；达到上限会警告并返回 `converged=false` |
| `tolMoment` | `1e-9` | 宏观矩更新容差 |
| `tolKinetic` | `1e-9` | 稳态动理学残差容差 |
| `tolWall` | `1e-9` | 壁面数据变化、壁面质量通量、平均密度误差容差 |
| `checkEvery`, `reportEvery` | `20`, `200` | 残差检查与进度输出间隔；结束前必检查 |
| `plot`, `save` | `true`, `true` | 实时图与结果保存开关 |
| `outputDir` | `results/bgk_dg` | 输出目录，不覆盖原 FV 结果 |

函数式 `source/initial` 收到的 `x,y` 为 `(4*Ncell)×1` 列向量，`vx,vy` 为 `1×Nvelocity` 行向量，应利用逐元素运算和隐式扩展返回 `(4*Ncell)×Nvelocity` 数组。源项按 Q1 节点插值后积分。非零源项必须满足域内总质量生产为零；稳态存在还取决于其他方程相容条件。线性化内域配完整壁 Maxwell 分布遵循 PDF 的混合模型，算例默认使用较小温差。

### 收敛判据与输出

停止条件同时要求宏观矩更新、动理学方程残差、壁面入射数据更新、实际滞后壁面通量及质量约束误差达标。动理学残差会用**当前**分布重新计算碰撞增益和物理壁面入射值，不能用本轮已经精确求解的冻结右端残差代替。

动理学残差为 DG 弱残差经质量矩阵逆映射后的域内均方根，使用 `M0` 速度权重和单元质量矩阵计算。壁面残差比较本轮实际使用的入射值与当前出射分布所要求的入射值，而非强制 DG 内侧多项式迹等于数值入射迹。壁面质量通量同样使用本轮的滞后入射与当前出射，因此在中间迭代通常不为零；在收敛时趋于零。质量归一化对实际入射值的仿射变换也计入了诊断。

结果默认写入 `results/bgk_dg`：

- `bgk_dg.mat`：完整 `result`，含 `H`、网格、速度、节点矩、格心矩、历史和诊断。
- `bgk_dg_profile.dat`：`x,y,rho_pert,ux,uy,T_pert,heat_x,heat_y`。
- `bgk_dg_history.dat`：矩更新、动理学残差、壁面残差、壁面质量通量、质量误差、累计秒数、迭代步数。未检查残差的行填 `NaN`。
- `bgk_dg_status.txt`：是否收敛及最终指标。

`result.H` 为 `(4*Ncell)×Nvelocity`；`result.moments` 为 `4×Ncell×4`，最后一维依次是密度扰动、`ux`、`uy`、温度扰动。`result.cellMoments` 是这些量的 Q1 单元平均值。平衡态这四个扰动均为零；完整密度为 `1+rho_pert`，线性化温度为 `1+T_pert`。`heat_x/heat_y` 是 `integral[v*(|v|²-2)*h*M0]`，是诊断量，不反馈 BGK 碰撞。线性模型不对全部速度处的完整分布保正。

### 验证

```matlab
addpath('debug');
check_bgk_dg
```

验证包含：原始网格和等温平衡态保持、奇数速度网格及反射配对、离散碰撞守恒、与独立 2×2 Gauss 积分/逐单元扫掠的一致性、兼容源项的非平凡 Q1 制造解、不兼容质量源项拒绝，以及防止仅凭宏观矩更新误判收敛。若完整默认算例已保存，还会检查其收敛、有限值和壁面质量通量。日志在 `debug/bgk_dg_validation.log`，完整默认运行日志在 `debug/bgk_dg_convergence.log`。这些验证不代表速度阶数或网格独立性研究。

2026-10-15 在 MATLAB R2025b 上完成原始网格的默认算例：6989 次源迭代，求解循环约 147.72 秒，动理学残差 `9.8804e-10`，壁面残差 `4.7272e-11`，最大壁面质量通量 `9.5784e-13`，平均密度误差 `2.22e-16`。轴向质量、两方向动量、总能量通量的最大变化分别为 `8.10e-11`、`1.86e-11`、`3.47e-17`、`5.44e-10`。独立逐单元扫掠对比误差 `1.55e-15`，Q1 制造解最大误差 `2.78e-16`。验证还执行收敛解重启，检查绘图及 MAT/DAT/状态文件读写，并生成 `results/bgk_dg/bgk_dg_summary.png`。

收敛容差采用 `M0` 加权范数，不能解释为每个速度处的未加权 `h=f/M0-1` 误差都同样小。例如本次结果的同 x 位置上下节点差的加权 RMS 为 `3.36e-11`，高速度尾部的未加权最大差仍为 `3.53e-5`；这也记录在验证日志中。

## 独立 1D2V Couette / Poiseuille

在 MATLAB 本目录中直接运行 `mono1D2V` 或按 F5，默认分别计算 Couette 和 Poiseuille，并显示两张监视图。需要保留结果时：

```matlab
result = mono1D2V;
c = result.couette;
p = result.poiseuille;

% 单独计算；delta 越大越接近连续流，但普通源迭代会变慢。
c = mono1D2V(struct('flow','Couette','delta',10));
p = mono1D2V(struct('flow','Poiseuille','delta',1,'acceleration',0.01));

% 网格和速度求积加密；velocityOrder 顺序为 [Nvx Nvy]。
p = mono1D2V(struct('flow','Poiseuille','cells',128, ...
    'velocityOrder',[256 8],'plot',false,'save',false));

% 复现旧的整区间 Hermite 速度求积（可使用奇数阶及零法向速度）。
r = mono1D2V(struct('normalQuadrature','hermite','velocityOrder',[32 8]));

% 两壁的 y 向速度顺序为 [左壁 右壁]；可同时给定驱动，测试叠加流。
r = mono1D2V(struct('flow','Couette','wallVelocity',[0 0.02], ...
    'acceleration',0.005));
```

只有 **x 方向物理网格**，两壁位于 `x=±gap/2`，壁速及驱动沿 y。分子速度为 `(vx,vy)`：`vx` 负责法向输运，主要速度剖面为 `uy(x)`，无需 `.cas` 文件。采用线性化 BGK、P1 迎风 DG（每单元两个独立自由度）及普通源迭代。法向默认在正、负速度半区分别使用 Gauss–Legendre 求积（截断到 `|vx|<=8`，Gaussian 已计入权重），流向使用 Gauss–Hermite；可选旧的整区间 Hermite 法向求积。这样的混合速度求积用于分开处理壁面入射/出射分布，相关研究见 [Ambruș 等的 Poiseuille 速度求积对比](https://www.sciencedirect.com/science/article/abs/pii/S1877750316300370)。这是独立的小扰动算例，不是原激波 Shakhov/GSIS 的通道版本，也不含黏性加热。

主求解器在 `mono1D2V.m` 内；默认额外调用独立文件 `channel1D2VReference.m` 生成**同参数、同漫反射壁面的 BGK 参考解**。关闭 `reference` 可单独运行主求解器。

归一化与 `monoBGKDG` 一致：

```text
M0 = exp(-vx²-vy²)/pi,  f = M0*(1+H)
vx*dH/dx = delta*(P(H)-H) + 2*a*vy
P(H) = rho_pert + 2*ux*vx + 2*uy*vy + T_pert*(vx²+vy²-1)
```

离散碰撞采用 Gram 矩阵投影以守恒四个碰撞不变量。`velocity.weight` 已含 `M0`；速度矩不应重复乘 Maxwell 分布。两壁均为 `T=1` 的完全漫反射壁，入射扰动为 `rho_wall_pert+2*vy*Uwall`；左壁入射为 `vx>0`，右壁入射为 `vx<0`，壁面密度由上一轮出射分布按零法向质量通量求得。域平均密度固定为 1，温度和密度输出为相对平衡态的扰动。

Poiseuille 用均匀 y 向驱动 `a`，源项 y 动量矩为 `a`，等效于参考密度下 `-dp/dy=a` 的线性压力驱动；一维法向模型不显式建立入口、出口或流向压降。Couette / Poiseuille 的线性化建模背景可参考 [Chalmers 的 BGK 论文](https://publications.lib.chalmers.se/records/fulltext/208317/local_208317.pdf)。

| 参数 | 默认值 | 含义 |
|---|---|---|
| `flow` | `'both'` | `'Couette'`、`'Poiseuille'` 或分别运行两者 |
| `cells`, `gap` | `64`, `1` | 均匀法向网格，范围 `[-gap/2,gap/2]` |
| `delta` | `1` | 碰撞频率；本文件不额外定义 Kn 与 delta 的换算 |
| `velocityOrder` | `[32 32]` | `[Nvx Nvy]`；默认法向正负半区各 16 点，流向 32 点，共 1024 个速度点；`halfLegendre` 要求 Nvx 为至少 32 的偶数 |
| `normalQuadrature` | `'halfLegendre'` | 可选 `'hermite'`，此时允许奇数阶和零法向速度 |
| `wallVelocity` | `[]` | Couette 默认 `[-0.01 0.01]`，Poiseuille 默认 `[0 0]` |
| `acceleration` | `[]` | Couette 默认 `0`，Poiseuille 默认 `0.01` |
| `initial` | `0` | 标量、`2*cells × Nvelocity` 数组或 `@(x,vx,vy)` |
| `maxIter`, `tol` | `30000`, `1e-5` | 迭代上限和宏观量相对更新阈值 |
| `stopCriterion` | `'reference'` | 参考代码式宏观量停机；`'equation'` 恢复原绝对更新/残差判据 |
| `checkEvery`, `reportEvery` | `10`, `100` | 稳态检查和进度输出间隔 |
| `plot`, `save` | `true`, `true` | 实时绘图和结果保存 |
| `reference` | `true` | 独立特征线积分方程 BGK 参考解及误差比较 |
| `referenceCells`, `referenceOrder` | `1024`, `160` | 参考解物理网格数（向上调整为 cells 的整数倍）、正速度半区求积阶数 |
| `outputDir` | `results/channel_1d2v` | 独立结果目录 |

默认采用 NonOscillatory 参考代码的单元宏观量相对更新判据，密度、温度和速度模三个指标的最大值小于 `tol=1e-5` 时停止。CIS 监测 kinetic 矩，GSIS 监测混合后源项矩；具体定义见下文“收敛诊断与验证”。方程残差和壁面诊断继续记录。设置 `stopCriterion='equation'` 并显式指定 `tol` 可恢复原先的绝对更新/方程/壁面联合停机条件。达到上限会警告并返回 `converged=false`。默认网格用于试算，尚不能代表速度和空间网格独立性；尤其稀薄流对法向零速附近的求积敏感。

每个算例保存到 `results/channel_1d2v/couette` 或 `poiseuille`，再次运行同名算例会覆盖该目录下的结果文件：

- `channel.mat`：完整结果变量 `r`。`r.H` 为 `2*cells × Nvelocity`，相邻单元端点具有独立 DG 迹；`r.cellMoments` 为 `[rho_pert,ux,uy,T_pert]`。
- `r.kineticTrace`：`(cells+1) × Nvelocity` 的实际迎风面分布，包含两侧本次扫掠使用的壁面入射值，供 GSIS 宏观高阶项提取使用。
- `profile.dat`：`x,rho_pert,ux,uy,T_pert,Pxy,heat_x,heat_y,uy_NS_no_slip`。
- `history.dat`：前七列为 `iteration,momentUpdate,kineticResidual,wallResidual,wallMassFlux,massError,elapsedSeconds`；GSIS 再有 `sourceUpdate,coupledResidual`，最后均追加 `referenceMacroUpdate,relativeDensityUpdate,relativeTemperatureUpdate,relativeSpeedUpdate`。未检查的方程/壁面诊断填 `NaN`，相对更新每步记录。
- `face_flux.dat`：`x,mass,streamwise_momentum,normal_momentum,energy`，使用实际迎风面迹；稳态应满足 `d(Pxy)/dx=a`，其余三种通量为常数。
- `status.txt`：收敛状态、坐标约定、流量、气体侧壁面速度及滑移量。`r.flowRate=integral(uy dx)`；`r.slipVelocity=气体侧端点 y 速度-壁速`。
- `reference_comparison.dat`：开启 `reference` 时保存 `x,uy_DG,uy_BGK_reference,difference`；两者均使用单元平均速度比较，不混用格心采样值。
- `summary.png`：开启 `plot` 且图窗仍存在时保存。左上为蓝色 DG、橙色同模型 BGK 参考解、灰色点线 NS 无滑移连续极限曲线；右上为 DG 减 BGK 参考解的误差；左下为收敛；右下为剪切应力守恒。NS 曲线使用 `mu=1/(2*delta)`，在 `delta=1` 时不作为 BGK 精度基准。

独立参考解令 `H=2*vy*phi(x,vx)`，将方程降为 `vx*dphi/dx+delta*phi=delta*uy+a`。在参考网格内把宏观速度表示为分片常数，对每个速度的输运方程沿特征线做解析积分，再对速度求积，直接求解宏观速度的稠密积分方程；不调用 DG 矩阵或源迭代。该参考解仍有空间和速度求积误差，不能称为解析精确解。验证分别加密参考解与 DG 解。

独立验证运行：

```matlab
addpath('debug');
check_channel_1d2v
check_channel_reference
check_couette_physics
```

检查包括：用独立 Gauss 体积分装配包含全部速度耦合及当前漫反射壁密度的整体 DG 线性系统，与源迭代终态比较；匀速平移平衡态和零法向速度；Couette 反对称与 Poiseuille 对称；动量守恒、碰撞不变量、壁面滑移、无碰撞极限及驱动反号；防止只凭小更新误判收敛；绘图、收敛解重启和输出文件读写。日志在 `debug/channel_1d2v_validation.log`。

2026-10-27 在 MATLAB R2025b 上通过两组检查。原有小网格整体 DG 对照仍使用 Hermite 速度点，最大分布差异为 Couette `1.03e-12`、Poiseuille `4.71e-12`；这验证了离散系统求解，**不能证明速度求积已收敛**。

同模型参考解核查发现：旧 `32×8` 整区间 Hermite 配置的 Couette、Poiseuille 速度相对 L2 误差分别为 **3.11%、2.81%**。此前验证的 `64` 单元、`128×8` 半区混合求积在 `delta=1` 时误差分别为 **3.73e-5、2.67e-5**；进一步加密至 `128` 单元、`256×8` 时分别为 **1.54e-5、1.09e-5**。参考解由 512 加密至 1024 单元的速度相对变化小于 `6.74e-7`，继续加密至 2048 单元小于 `1.69e-7`；正半区速度阶数由 160 增至 320 的相对变化小于 `1.20e-6`。这是此前 `delta=1` 两算例的验证范围，不覆盖当前 `32×32` 默认速度网格或任意 delta。

此前 `128×8` 速度网格的 Couette、Poiseuille 分别在 16、39 次源迭代收敛；Poiseuille 流量为 `0.0153867739441`，独立 BGK 参考为 `0.0153867816836`，无滑移 NS 值为 `0.00166666666667`。精度证据在 `debug/channel_reference_validation.log`，该次完整运行记录在 `debug/channel_1d2v_default.log`。

### Couette 的物理检查与近壁结构

纯 Couette、开启 `plot` 时，求解结束后额外显示并保存 `couette_physics.png`。也可直接检查已保存的默认结果，无需求解：

```matlab
plotCouette1D2V
% 或传入返回值：
checks = plotCouette1D2V(result.couette);
```

该图上排依次为 `uy(x)`、`ux(x)`、剪切应力 `tau_xy=integral(vx*vy*H*M0 dv)`；下排为密度/温度扰动、相对中央直线拟合的偏离、左壁放大。密度和温度扰动为零表示实际密度和温度约为 1。中央直线拟合固定使用 `|x|<=0.2*gap`，并同时显示 DG 与独立 BGK 参考相对于这条直线的偏离。

分别报告气体侧端点速度与壁速之差（`actualWallSlip`），以及中央直线外推至壁面后的速度与壁速之差（`extrapolatedWallSlip`）。这两者通常不同；近壁弯曲可通过两条速度的差观察。中央拟合是诊断工具，拟合区间不等于已证明的纯主体区，其余量也不是严格的半空间 Knudsen 层函数。壁面影响在较强稀薄条件下可能相互重叠，相关近壁结构与滑移的区别可参考 [BGK Knudsen 层研究](https://pure.strath.ac.uk/ws/portalfiles/portal/64962346/strathprints005019.pdf)。

`check_couette_physics` 对默认 `delta=1` 和另存的 `delta=10` 进行检查，后者保存在 `results/channel_1d2v/delta10/couette`，不覆盖默认参数。实测 `delta=1` 时，最大 `|ux|`、`|rho_pert|`、`|T_pert|` 分别为 `2.64e-19`、`5.65e-19`、`2.22e-19`，剪切应力相对跨度 `7.68e-9`，左壁气体相对壁面滑移 `4.98009e-3`，中央拟合外推滑移 `5.52986e-3`。`delta=10` 时，中央区拟合相对 L2 偏差由 `2.77e-3` 减至 `4.80e-4`，左壁实际滑移为 `1.19344e-3`，近壁相对主体直线仍有可分辨偏离；DG 相对同模型参考的速度 L2 误差为 `2.04e-5`。

物理检查量、中央拟合和偏离值另存为 `couette_physics.mat/.txt/.dat`，验证日志为 `debug/couette_physics_validation.log`。

### Poiseuille 增大 delta 的对照

```matlab
r = runPoiseuilleDelta10;
```

保持 `acceleration=0.01`、静止漫反射壁面及宽度 1，分别计算 `delta=1` 和 `10`，两者使用相同的 128 个物理单元与 128×8 速度点。独立 BGK 参考采用 2048 单元、正半区 320 点速度求积；另外以 256 单元、256×8 速度点复核 `delta=10`。新结果在 `results/channel_1d2v/delta10/poiseuille`，原默认结果和 `delta=10` Couette 保留。

2026-10-27 MATLAB R2025b 实测：

| delta | 最大单元平均速度 | 气体侧壁面速度 | 气体侧壁面速度/峰值速度 | 流量 | 源迭代步数 |
|---|---:|---:|---:|---:|---:|
| 1 | 0.0174372 | 0.0101844 | 58.4% | 0.0153867830 | 39 |
| 10 | 0.0365591 | 0.00763897 | 20.9% | 0.0276863691 | 579 |

增大 delta 后，归一化速度形状更接近无滑移抛物线，但仍有明显滑移。固定加速度时，本模型 `mu=1/(2*delta)`，所以速度幅值也会变化。`delta=10` 无滑移 NS 流量为 `0.0166666667`，当前 BGK 流量比它高约 66.1%；不能把该参数下的 BGK 解视为已经等同于无滑移 NS 解。

`delta=10` 主网格相对独立 BGK 参考的速度 L2 差异为 `2.48e-5`，加密网格为 `9.50e-6`；DG 加密相对变化 `2.22e-5`，参考网格 1024→2048 单元的相对变化 `1.12e-5`，参考解自身仍有离散误差。最大 `|ux|/|rho_pert|/|T_pert|` 为 `2.73e-18`，最大守恒通量平衡缺陷 `9.05e-11`。

`delta_comparison.png` 给出原始速度、按各自最大单元平均速度归一化的剖面、剪切应力和收敛历史；速度曲线包含气体侧壁面迹。数值在 `delta_comparison.dat/.txt`，完整过程在 `run.log`；本次加密验证不覆盖其他参数。

## 1D2V BGK / HDG 的参考代码式 GSIS

当前宏观离散和耦合按桌面 `NonOscillatory/Synthetic_Acceleration.f90` 的方式适配到 **1D2V BGK**。保留二维 Gaussian 矩和热流系数，未照搬三维 Shakhov 的 3/2、4/3、5/4 等系数。按用户要求，稳定参数保持 **`[1 delta delta delta]`**，不采用参考 Fortran 当前启用的常数参数。

```matlab
% 原 CIS 入口和默认行为：
cis = mono1D2V(struct('delta',10,'plot',false,'save',false));
% 新的完整 GSIS 外层入口，默认也分别运行两个流：
gsis = monoGSIS1D2V(struct('delta',10,'plot',false,'save',false));
% 生成 delta=1/10、相同网格上的 CIS/GSIS 对照：
[summary,runs] = runGSISChannelComparison;
% 可单独加密空间网格：
runGSISChannelComparison(struct('deltas',10,'cells',128, ...
    'outputDir',fullfile(pwd,'results','gsis_channel_comparison','N128')));
```

默认空间 64 单元、速度 `32×32`、delta=1、`stopCriterion='reference'`、tol=`1e-5`，Couette 壁速 `[-0.01 0.01]`，Poiseuille 加速度 `0.01`，两者均为等温完全漫反射壁。`monoGSIS1D2V` 复用原 kinetic DG sweep，默认结果目录为 `results/channel_1d2v_gsis`；CIS 的结果目录保持 `results/channel_1d2v`。对照脚本默认关闭独立积分参考，以直接比较同网格的两个离散稳态，并保存 CSV、MAT 和日志。

### 宏观方程、面通量与变量

HDG 单元变量为 `[p ux uy T Lx Ly E]`，其中 `p=rho+T`，rho/T 均为扰动。使用两倍尺度的动量方程，并用质量方程从能量方程中减掉质量通量：

```text
d(ux)/dx = 0
d[p - Lx/delta]/dx = 2*ax
d[-Ly/delta]/dx = 2*ay
d[-E/delta]/dx = 0
Lx - d(ux)/dx = Pi_xx
Ly - d(uy)/dx = Pi_xy
E  - d(T)/dx  = Theta_x
```

物理量 `sigma_xx=-Lx/delta`、`sigma_xy=-Ly/delta=2*Pxy`、`qx=-E/delta`。每单元 7 个 P1 场、14 个自由度；全局求解 `p_hat,ux_hat,uy_hat,T_hat`，然后局部恢复。向外数值通量为：

```text
mass:             ux*n + s_p*(p-p_hat)
normal momentum:  (p_hat-Lx/delta)*n + s_ux*(ux-ux_hat)
shear momentum:   -Ly/delta*n + s_uy*(uy-uy_hat)
heat:             -E/delta*n + s_T*(T-T_hat)
```

质量通量中的 ux 使用单元内值；动量中的压力使用公共面变量。内部面通过两侧向外通量之和为零确定公共面变量，不指定为 kinetic 迎风面的宏观矩。

```matlab
mesh = struct('edges',linspace(-0.5,0.5,65)');
delta = 10;
bc = struct('trace',[0 0 -0.01 0; 0 0 0.01 0]); % [p ux uy T]
param = struct('stabilization',[1 delta delta delta]); % 默认值保持不变
[macro,aux,hdg] = hdg_solve_macro_1d2v( ...
    mesh,delta,[0 0],0,0,0,bc,param);
```

`macro.p` 为压力，`macro.rho=macro.p-macro.T`；对接 kinetic 的 `macro.nodal` 和 `macro.cellAverage` 仍为 `[rho ux uy T]`，`macro.pressureNodal` 和 `hdg.Uhat` 为 `[p ux uy T]`。`hdg.localFields` 第一块也为 p。`hdg.outwardFlux` 的动量通量已乘 2，因此与 kinetic 的 `Pxy` 比较时需除以 2。体源项输入仍是物理加速度 `[ax ay]`，求解器内部将动量源乘以 2。

`mesh.edges` 支持非均匀网格；Pi_xx/Pi_xy/Theta_x 支持常数、2*N 节点列向量、N×2 数组或函数。边界可为结构或回调：`bc.trace` 是左右两行 `[p ux uy T]`；NaN 分量须由同位置 `bc.flux` 指定向外数值通量。纯通量条件须另外处理零空间和相容性。缓存含 `pressure-reference-v1` 版本标记，旧密度形式缓存不能混用。

### 高阶项、壁面和外层回传

适配器从 H 提取 Rxx/Rxy/RT，使用单元内部导数 `M*Rx=G'*R`，对应参考 Fortran 的体积分导数；不加入 kinetic 迎风面项。`project_hot_derivatives_1d` 的 weak 模式仍可独立使用，但当前适配器明确选择 strong/interior 模式。

宏观壁面取相邻单元 `H([1 end],:)` 的矩，不使用完整漫反射数值迹：

```text
p_hat  = rho_kin,interior + T_kin,interior
ux_hat = 0
uy_hat = uy_kin,interior
T_hat  = T_kin,interior
```

kinetic 本身继续使用完全漫反射边界；宏观边界保留气体切向滑移。第一轮 kinetic sweep 后，适配器可直接调用：

```matlab
[macro,aux,hdg,hot] = gsis_macro_from_kinetic_1d2v(cis.poiseuille);
[macro,aux,hdg,hot] = gsis_macro_from_kinetic_1d2v( ...
    cis.poiseuille,struct('cache',hdg.cache));
```

完整 GSIS 每轮执行：由独立保存的 Qsource 组装碰撞增益，解 kinetic H，取 Qkin，解宏观方程，再更新下轮源项：

```text
Kn_cell = sqrt(pi)/(2*delta*dx)
omega   = min(1,1/Kn_cell)
Qsource = Qkin + omega*(Qmacro-Qkin)
```

宏观步骤不直接修正 H。保留原有闭合通道的密度均值为零规范，并对 Qsource 的密度均值作同样约束。`r.nodalMoments/cellMoments/flowRate` 始终来自实际 kinetic H；`r.synthetic` 另外保存宏观解、辅助量、面变量、高阶项、omega 和 `sourceMoments`。自定义稳定参数可通过 `options.macroParam.stabilization` 传入。

### 收敛诊断与验证

默认停机与参考 Fortran 的 `Calculate_Residual_rho_T_U` 对齐。每个单元计算 `R_i=int(1+rho_pert)dx`、`T_i=int(1+T_pert)dx`、`V_i=hypot(int(ux)dx,int(uy)dx)`；对这三个单元向量分别取 `norm(A_new-A_old)/norm(A_new)`，最大值小于 `1e-5` 时停止。CIS 使用 Qkin，GSIS 使用混合后的 Qsource，与参考代码的宏观量回传阶段一致。首次 OLD 取静止参考气体：密度/温度为单元长度，速度为零。零分母且无变化记为零，有变化记为 Inf。实现见 `reference_macro_update_1d2v.m`。

GSIS 的 `diagnostics.coupled` 是将当前混合源 Qsource 代入 kinetic DG 方程得到的残差；`diagnostics.kinetic` 则使用 Qkin 自身构造碰撞源，表示原 CIS 方程残差。默认模式下这些残差及壁面诊断只记录，不参与停机；因此宏观量更新小于 tol 不代表方程残差也小于 tol。`stopCriterion='equation'` 恢复原来的 kinetic 矩更新、源项更新、耦合残差、壁面残差、壁面质量通量和密度均值误差联合判据。历史文件保留原有列，在最后追加四列参考相对更新指标；`r.convergence` 记录所用判据、阈值和最终指标。

`runGSISChannelComparison` 的结果表记录完整外层迭代次数、两种残差、速度剖面和流量相对 CIS 的差异、实际稳定参数，以及停机判据、阈值和最终相对更新。相同宏观量停机标准并不保证最终离散解完全相同，也不表示相同运行时间。`demo_gsis_macro_1d2v` 仍是 CIS 收敛后仅执行一次宏观求解的子模块示例，其步数不是完整 GSIS 步数。

2026-11-05，统一参考停机标准后的实测结果（64 单元、32×32 速度点、tol=1e-5，稳定参数仍为 `[1 delta delta delta]`）：

| 流动 | delta | CIS 步数 | GSIS 步数 | GSIS 单元平均速度相对 CIS 的 L2 差异 |
|---|---:|---:|---:|---:|
| Couette | 1 | 11 | 11 | 0.001515% |
| Poiseuille | 1 | 23 | 23 | 0.009321% |
| Couette | 10 | 82 | 47 | 0.07736% |
| Poiseuille | 10 | 251 | 66 | 0.3003% |

完整结果保存于 `results/gsis_channel_comparison/reference_stopping/summary.csv`。相同停机条件下，CIS 后只做一次 HDG 的原始速度差异依次为 0.1890%、0.3441%、0.09709%、0.2395%，保存在 `results/gsis_macro_1d2v/converged_cis_one_step/reference_stopping/summary.csv`。停机对齐没有消除宏观解与 CIS 的差异。

```matlab
addpath('debug');
report_hdg_after_converged_cis_1d2v % CIS 收敛后只解一次 HDG，直接比较，不乘 omega
audit_macro_wall_reference_1d2v    % 对照参考壁面数据及边界施加
audit_hdg_bulk_interface_1d2v      % 独立高阶项、面公式和冻结宏观问题解析积分
check_hdg_one_step_refinement_1d2v % 原始单步差异的 32/64/128/256 单元检查
check_hdg_macro_1d2v          % 制造解、缓存、混合边界、解析流、二阶精度
report = audit_hdg_weak_form_1d2v; % 独立积分组装、通量、界面公式、2V 系数
check_gsis_channel_1d2v      % 不修改 H、下轮使用混合源、平移平衡态、停机检查
check_channel_1d2v           % 原 CIS 对整体 kinetic DG 解、物理和输出回归
check_reference_stopping_1d2v % 参考宏观相对更新、源项监测、零分母与默认阈值
```

独立核查的参考矩阵由三点 Gauss 积分逐基向量生成，不调用生产 F/G/M 或 Kronecker 装配。涵盖非均匀网格、delta=0.07/1/7/1000、默认及自定义稳定参数、法向/切向驱动、界面平均/跳跃公式和二维 Gaussian 恒等式。`fixed_point_audit.csv` 是将 kinetic 解和指定的测试面变量代入宏观方程的诊断，不能把这些测试面变量等同于 HDG 自由求出的内部面变量；参考式单元内部高阶导数下，辅助残差也不必为零。对完整耦合的判断应结合实际收敛、与 CIS 的稳态偏差和网格加密结果。

以下保留 2026-11-03 的严格方程停机历史结果：64 空间单元、32×32 速度点、tol=1e-10。现在复现须显式设置 `stopCriterion='equation'`、`tol=1e-10`；严格离散核查脚本也显式使用方程判据。

| 流动 | delta | CIS 步数 | GSIS 步数 | GSIS 单元平均速度相对 CIS 的 L2 差异 |
|---|---:|---:|---:|---:|
| Couette | 1 | 16 | 16 | 0.00151% |
| Poiseuille | 1 | 39 | 38 | 0.00938% |
| Couette | 10 | 160 | 87 | 0.0733% |
| Poiseuille | 10 | 579 | 139 | 0.323% |

四个 GSIS 耦合残差均小于 1e-10，原 CIS 残差则约为 1.28e-7、1.33e-6、6.68e-6、5.45e-5，不能将耦合收敛说成原 CIS 方程也收敛到 1e-10。delta=10 时，32→64→128 单元的速度相对 CIS 差异：Couette 为 0.1173%→0.0733%→0.0439%，Poiseuille 为 0.4298%→0.3227%→0.2256%；差异下降，但这组数据未体现二阶下降。空间加密同时按参考公式改变 omega。完整数值在 `results/gsis_channel_comparison/summary.csv`、`refinement.csv` 及相应子目录。

单步检查与完整耦合须分别报告：`report_hdg_after_converged_cis_1d2v` 先将 CIS 收敛，再仅解一次 HDG，直接比较原始宏观量，不乘 omega、不再推进 kinetic。相同 64 单元/32×32 速度点下，单元平均 uy 的相对 L2 差异为：delta=1 Couette 0.1890%、Poiseuille 0.3442%；delta=10 Couette 0.0926%、Poiseuille 0.2584%。rho/ux/T 的最大绝对差均小于 5e-17。该测试的 CSV/MAT/日志单独保存于 `results/gsis_macro_1d2v/converged_cis_one_step`，不使用上表的耦合稳态差异替代此测试。

进一步核查固定这份 Pi 和气体侧壁面速度，独立积分 `L'=-2*delta*a, uy'=L-Pi`。P1 kinetic 单元内的 Pi 为常数，因此可得到逐单元精确二次速度及其精确单元平均。64 单元时，Couette HDG 与该直接积分解的相对差约为 1e-14，但两者相对 CIS 的差异相同；Poiseuille delta=10 时，HDG 相对直接积分解的差异约 0.00291%，远小于相对 CIS 的 0.2584%。这将主要差异定位到冻结高阶项后的宏观子问题与 kinetic DG 之间，而非 HDG 线性系统或面公式的实现误差。

直接按 Fortran 的逐速度体积分顺序重算 Pi，负载差小于 1.2e-17；检查了两项 `Pi_xy/2` 在一维合并后必须得到完整 Pi_xy。实际 HDG 辅助方程和内部面平均/应力跳跃公式均满足到舍入误差。kinetic 剪应力矩的 DG 弱式残差约 1e-17，而其单元内强闭合有非零缺口；该缺口可由实际迎风面的 `J_sigma_hat-J_sigma_interior` 修正重构，恒等式误差约 1e-16，其中 `J_sigma=<2*vx^2*vy*H>=uy+Rxy`。参考式单元内导数没有自动恢复这些 DG 面贡献；不能由此推断仅替换 Pi 的导数就能让整套宏观离散精确保留 CIS。

原始单步加密检查（固定 32×32 速度点，不使用 omega）中，64→128→256 单元的相对速度差异：delta=1 Couette 为 0.1890%→0.1159%→0.0643%，Poiseuille 为 0.3442%→0.2040%→0.1113%；delta=10 Couette 为 0.0926%→0.0874%→0.0730%，Poiseuille 为 0.2584%→0.2264%→0.1791%。差异下降，但本组未达到二阶。完整记录为 `bulk_interface_audit.csv` 和 `one_step_refinement.csv`，稳定参数保持 `[1 delta delta delta]`。
