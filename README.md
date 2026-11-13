# MATLAB 有限体积正激波 GSIS 算例

新增二维顶盖驱动方腔入口 [`monoCavity2D2V.m`](monoCavity2D2V.m)（CIS）及 [`monoGSISCavity2D2V.m`](monoGSISCavity2D2V.m)（GSIS），配置、二维宏观方程和对比结果见下面的“二维 BGK 方腔”。

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

## 二维 BGK 方腔

2026-11-11：`monoCavity2D2V` 在 `[0,1]×[0,1]` 上求解线性化二维速度 BGK：

```text
vx*dH/dx + vy*dH/dy = delta*(P(H)-H)
M0=exp(-vx^2-vy^2)/pi
P(H) 由同一离散 Gram 投影保持质量、两方向动量和能量矩
```

顶壁沿 x 正向运动，`lidSpeed=0.01`；其余三壁静止。四壁完全漫反射、等温，壁温扰动为零，无体力。入射分布 `H_in=rhoWall+2*(vx*Uw_x+vy*Uw_y)`，壁面密度由同面出射分布及零法向质量通量确定。顶壁和侧壁在角点处各自按面法向提供入射条件，不对角点壁速作平均。不把气体速度强制设成壁速，因此可有滑移。

空间为矩形 Q1-DG，每单元 4 个节点，两个空间方向都离散；内部面按每个离散速度迎风。默认均匀 24×24 单元、32×32 referenceCubic 速度点，两速度都截断于 ±4、零附近加密。`meshStretch>0` 可向四壁加密，`quadrature='hermite'` 可切换为 Gauss-Hermite。`monoCavity2D2V` 默认使用 CIS，碰撞源和漫反射壁面密度滞后一轮，闭合方腔密度扰动体平均约束为零。二维 GSIS 使用下述独立入口。

停止时要求宏观相对更新和原 kinetic 方程相对残差都小于 `tol=1e-7`，最后一次 sweep 的壁面质量通量小于 `tol*abs(lidSpeed)`。密度/温度更新以完整参考量归一化，速度以当前速度向量范数归一化；方程残差使用速度权重和空间逆质量矩阵形成 L2 范数，除以共同首步残差。每步都用当前 H 自身的矩和当前漫反射映射重算原方程残差，未用冻结源局部求解残差代替。

```matlab
r = monoCavity2D2V(struct('delta',10,'cells',[32 32]));
addpath('debug');
check_cavity_2d2v;  % 独立 Gauss 积分整体 kinetic 线性系统、质量通量、反向驱动、对称性
summary = run_cavity_2d2v; % delta=1/10，24×24 与 32×32，固定 32×32 速度点
```

输出位于 `results/cavity_2d2v/delta_10_N32x32` 等独立目录。`cavity.png` 依次展示速度大小/箭头、流线、竖直中心线 ux、水平中心线 uy；`profile.csv` 为单元中心 `[x,y,rhoPert,ux,uy,TPert,Pxy]`，`history.csv` 为更新及方程残差历史，`cavity.mat` 保存完整分布与配置。总表为 `summary.csv`，空间加密对照为 `resolution.csv`。主涡位置由对称轴 x=0.5 上 ux 的负到正零点估计。空间加密差异只比较两张网格的中心线，不作为连续解或文献精度证明。

独立 2×2 小网格、4×4 Gauss-Hermite 速度点核查中，迭代 H 与一次性组装的完整 kinetic 方程直接解最大差约 4.4e-13；四壁零质量通量、碰撞守恒、反向驱动线性关系、静止平衡及方腔对称性均通过。

32×32 空间网格实测：delta=1、10 分别 22、160 步，原 kinetic 绝对残差 4.6426e-10、1.3064e-9，相对首步残差 1.6398e-8、9.1910e-8。两者均出现顺时针主涡，涡心近似为 (0.5,0.7010)、(0.5,0.7419)。24→32 空间加密时，竖直中心线 ux 的相对 L2 变化为 0.7625%、0.1235%，水平中心线 uy 为 0.9847%、0.1770%；迭代次数均保持 22/160。尚未进行速度网格独立性和外部文献比对。

指定的 delta=100、纯 CIS、1e-5 角涡测试使用 32×32 的双向 tanh 加密网格，`meshStretch=2`，速度仍为 32×32 referenceCubic。复现：

```matlab
r=monoCavity2D2V(struct('cells',32,'meshStretch',2,'delta',100, ...
    'iterationMethod','CIS','tol',1e-5,'maxIter',30000,'reportEvery',200, ...
    'outputDir',fullfile(pwd,'results','cavity_delta100_cis_tol1e5')));
addpath('debug'); analyze_cavity_delta100_cis;
```

`bottom_corners.png` 放大左右底角；上排为原始速度流线，下排为实际 kinetic 面质量通量重构的流函数。`vortex_candidates.csv` 记录原始单元平均速度双线性插值的椭圆型零点，并非仅凭插值流线判定涡旋。流函数通过边通量最小二乘重构，投影误差另存，不能把弱于此误差的等值线作为可靠角涡。该测试不使用 GMRES，也不使用 GSIS。`check_cavity_2d2v` 另独立核对了四个速度符号象限的显式扫描与直接输运解，最大差约 7.1e-15。

本次 delta=100 纯 CIS 在第 2803 步通过停止条件，宏观相对更新 5.0033e-7，原 kinetic 相对残差 9.9962e-6、绝对残差 3.3166e-7。可以分辨左右底角各一个逆时针小涡，原始速度插值的涡心约为 (0.02787,0.02811)、(0.97213,0.02811)，主涡顺时针，中心约为 (0.5,0.76205)。两角流函数为正、主涡为负，角涡流函数约 1.55e-8，边通量流函数投影的最大误差约 2.69e-10。最小单元宽度 0.005177。该单网格、单速度求积测试观察到一对角涡，尚未证明其精确位置/强度收敛，也未分辨更高阶涡列；不作为二维 GSIS 验证结果。

### 二维方腔 GSIS

`hdg_precompute_q1_2d2v` 和 `gsis_macro_from_kinetic_2d2v` 实现真正的二维 Q1-HDG 宏观子问题。每单元 8 个字段 `[p,ux,uy,T,Lxx,Lxy,Ex,Ey]`、每字段 4 节点，每面有 `[p,ux,uy,T]` 的 2 个公共节点；`p=rho+T`。物理应力和热流为 `sigma_xx=-Lxx/delta`、`sigma_xy=-Lxy/delta`、`qx=-Ex/delta`、`qy=-Ey/delta`。

```text
div(u) = 0
d_x(p-Lxx/delta) - d_y(Lxy/delta) = 0
-d_x(Lxy/delta) + d_y(p+Lxx/delta) = 0
div(-E/delta) = 0
Lxx = d_x ux - d_y uy + Pi_xx
Lxy = d_y ux + d_x uy + Pi_xy
Ex  = d_x T + Theta_x; Ey = d_y T + Theta_y
```

2026-11-12，当前二维默认算法按最新 `GSIS_Gray_20261110.zip` 的宏观离散及迭代结构对齐。已逐字节核对解压的 `Synthetic_Acceleration.f90`、`Velocity_Distribution.f90` 和主程序与原 ZIP 相同。参考程序是 Callaway 声子模型，当前仍为 2V BGK；物理矩定义、碰撞系数、密度零空间以及矩形 Q1 网格保留 BGK 设置，不能称为相同物理算例或逐系数移植。

对应关系为：参考宏观标量 T 对应此处压力 p，参考矢量 q 对应此处速度 u；质量数值通量使用单元内法向速度，动量通量压力使用公共 trace。四个高阶项使用体积分及 kinetic 迎风面的弱导数。特别地，高阶项的壁面入射部分使用本次输运扫描实际使用的滞后壁面密度（随 kinetic 密度规范同步平移），而宏观仿射壁面通量由扫描后的 H 重新漫反射得到。两者的取值时刻不同，与参考代码的 `Evaluate_Upwind_Face_Load` 和 `Calculate_Kinetic_Boundary_Flux_Moments` 分别对应。单独提供 H 而没有 `transportWallDensity` 的宏观调用会重构当前壁面入射量，`hdg.hotUsesTransportWallInflow` 标明实际路径。

默认稳定参数改为 `[1 1 1 1]`，与参考 `ST=1` 相同且不随 delta 改变。四壁采用 `F_HDG=F_kinetic+Jwall*(U-Ukinetic)`，由 BGK 离散漫反射规则和碰撞平衡态切线构造 Jacobian。通过体平均密度约束消除零空间，不删除守恒方程。每轮将整个 `Qmacro=[rho ux uy T]` 回传下一轮碰撞源，不改写 H。默认 `kineticConsistency=false`、`relaxation=1`，不添加内部残差匹配载荷或自动选择回传松弛。

二维 GSIS 入口默认 `tol=1e-5`，按返回宏观量的单元积分更新停止，移植参考 `Update_Macro_Residual` 的逐分量范数和小分母分支；BGK 使用全密度 `1+rho`、全温度 `1+T` 及单位热速度尺度下的 ux/uy。每个分量先计算单元积分向量 a，令 `r0=max(1,max_k norm(a_k)^2)`；若 `norm(a_k)^2>tol^2*r0`，使用 `norm(a_k-a_old_k)/norm(a_k)`，否则除以 `sqrt(r0)`；最后取四分量最大值。原 kinetic 残差、耦合残差和壁面质量通量继续记录，不参与这个默认停机判断。`r.convergenceCriterion='referenceCellIntegratedMacroUpdate'`，`r.diagnostics.originalConverged` 独立标记原方程是否达到阈值。`history.csv` 第 10-13 列保存四个停机分量，不能将宏观更新停止称为原方程已收敛。

原有八矩残差匹配和回传松弛仍可显式设置，属于另一个实验版本；`kineticConsistency=true` 或 `requireOriginalResidual=true` 会额外要求原方程和耦合方程残差通过。这些选项不属于当前参考默认算法。若需旧稳定参数，显式传入 `macroParam.stabilization=[1 delta delta delta]`。

```matlab
g=monoGSISCavity2D2V(struct('cells',32,'meshStretch',2,'delta',100, ...
    'tol',1e-5,'reportEvery',20));
addpath('debug'); check_gray_defaults_cavity_2d2v; check_gsis_cavity_2d2v;
summary=run_gray_aligned_cavity_2d2v; % 当前默认值，独立保存 delta=100/1000
```

最新默认配置复测：32×32 Q1、meshStretch=2、32×32 referenceCubic 速度、顶盖速度 0.01、零初值、tol=1e-5。

| delta | 宏观更新停止步数 | 宏观更新量 | 原 kinetic 相对残差 | 耦合相对残差 | 实际 kinetic 速度相对保存的 CIS 差异 |
|---|---:|---:|---:|---:|---:|
| 100 | 14 | 5.7450e-6 | 1.3012e-2 | 1.1323e-3 | 0.15622% |
| 1000 | 7 | 9.0239e-6 | 1.3310e-2 | 3.9222e-4 | 未运行同配置 CIS |

这两组只通过参考宏观更新停止，原 kinetic 方程没有达到 1e-5。宏观壁面通量关系误差分别为 1.63e-19、1.08e-19；这是宏观线性子问题误差，不是 kinetic 方程残差。完整结果、逐步四分量更新、H 和最后扫描入射密度位于 `results/cavity_gsis_gray_defaults_2d2v/`，总表为 `summary.csv`。独立检查覆盖实际扫描入射量的高阶面矩、宏观壁面取值时刻、四法向通量、Jacobians、凝聚与整体求解，以及完整回传但不修改 H。

以下为历史结果，使用旧的 delta 稳定参数和旧停机配置，不代表当前默认算法。2026-11-11，同一 delta=100、32×32 角点加密空间网格、32×32 速度点、零初值、tol=1e-5：

| 方法 | 外迭代步数 | 原 kinetic 相对残差 | 速度相对 CIS 差异 | 底角小涡 |
|---|---:|---:|---:|---:|
| CIS | 2803 | 9.9962e-6 | — | 2 |
| GSIS，参考结构、omega=1 | 18 | 6.0336e-2 | 0.48689% | 2 个速度零点，但原方程未通过 |
| GSIS，八矩一致性修正、omega=0.1 | 100 | 2.5682e-6 | 0.019544% | 2 |

后一版本原方程绝对残差 8.5209e-8，耦合相对残差 9.5552e-6。左右角涡中心约为 (0.02780,0.02819)、(0.97220,0.02819)，主涡约为 (0.5,0.76197)，保持与 CIS 相同的旋转结构。该差异包含 CIS 本身有限停止阈值带来的误差，不等同于相对连续解误差。

正式对比位于 `results/cavity_gsis_comparison_delta100/summary.csv`、`one_macro.csv`、`comparison.mat`，两个版本分别保存在 `reference/` 和 `consistent/` 子目录，含 `cavity.png` 和 `bottom_corners.png`。独立验证覆盖：二维方程的 Gauss 积分组装、四个法向、壁面 Jacobian 有限差分、整体与凝聚解、八矩残差匹配、CIS 固定点、回传松弛公式及不修改 H。一般网格/速度独立性、外部文献精度和更高阶角涡仍需另测。

2026-11-12，按用户指定只运行 delta=1000 的参考结构版 GSIS：`macroParam.kineticConsistency=false`、`macroParam.relaxation=1`，不加入内部弱残差修正和回传松弛；沿用 BGK 适配版稳定参数 `[1 delta delta delta]`，不等同于声子参考程序的三个 ST 数值。顶壁速度 0.01，双速度 32×32 referenceCubic，tol=1e-5。

```matlab
cfg=struct('cells',64,'meshStretch',3,'delta',1000,'tol',1e-5, ...
    'macroParam',struct('kineticConsistency',false,'relaxation',1), ...
    'outputDir',fullfile(pwd,'results','cavity_gsis_delta1000_reference'));
r=monoGSISCavity2D2V(cfg);
% 本次还运行了 cells=32、meshStretch=2 的较粗网格；两张网格都完成后：
addpath('debug'); analyze_gsis_delta1000;
```

32×32/stretch2 用 9 步达到耦合停止条件，仍只有左右各一个角部椭圆型速度零点；64×64/stretch3 用 18 步，最小单元宽度约 5.1094e-4，识别到每侧两个候选角涡。较大底角涡中心约 (0.03620,0.03656)、(0.96380,0.03656)，逆时针；新出现的较小一对约 (0.001290,0.001249)、(0.998710,0.001249)，顺时针。`corner_multiscale.png` 直接放大原始速度场，不将重构流函数作为额外速度修正。候选小涡的位置与数量尚未通过进一步网格/速度/停止阈值核查。

较细网格原 kinetic 相对残差仍为 2.2635e-2，而耦合相对残差为 6.5546e-6；原方程阈值未通过。流函数全域边通量投影误差约 6.14e-8，比第二级涡的重构流函数量级大，因此不能仅凭流函数图或一次双线性速度插值就确认更高阶涡的物理精度。结果单独保存于 `results/cavity_gsis_delta1000_reference/summary.csv`，两种网格子目录含原始分布、残差历史、候选中心及多尺度放大图。本次未运行 CIS，也未运行一致性修正版。

### 五次映射空间网格的 delta=1000 角涡检查

2026-11-12 启动、2026-11-13 复核：使用用户截图的五次坐标映射 `x(s)=10s^3-15s^4+6s^5`，两方向对称加密，保留 Q1 四边形、32×32 referenceCubic 速度网格、顶盖速度 0.01、Gray 默认稳定参数 `[1 1 1 1]`、完整回传及无内部残差修正。入口 `addpath('debug'); run_quintic_cavity_delta1000(64)`；独立原始速度复核及绘图为 `python3 debug/analyze_quintic_cavity_delta1000.py`。输出单独存入 `results/cavity_delta1000_quintic_20261112/`。

新增可选参数 `meshMapping='quintic'`，默认仍为原 `tanh` 映射；`minIter` 默认为 0；`tolCoupled` 默认为 Inf，不改变默认宏观停机。此次显式设 `minIter=60`、`tolCoupled=1e-9`，不改变宏观更新公式及小分母分支。64×64 五次映射最小单元宽度 3.725849e-5，最大 0.02927781；原 tanh/stretch3 最小宽度约 5.10936e-4。第 24 步首次通过宏观更新阈值，但耦合残差约 0.04657；第 160 步才达到耦合相对残差 9.634259e-10，宏观更新 1.857143e-13，原方程相对残差仍为 0.005147204。

单元中心原始 kinetic 速度的双线性插值，分别经 Newton 法和消元多项式求根复核，只识别到左右底角各两层候选角涡，没有第三层。左侧涡心约 (0.03609525,0.03654385)、(0.001181812,0.001181895)，右侧镜像；较大涡逆时针、较小涡顺时针。此结果仍有原 kinetic 方程缺口，不能据此证明物理上不存在第三层。流函数边通量投影误差约 1.9421e-8，不宜用更弱的重构流函数极值单独确认微涡。报告图为 `corner_multiscale_independent.png`，检查覆盖原有整体 DG 对照、五次映射平衡态、CIS/GMRES 一致性和显式耦合停机选项，日志为 `validation.log`。

### 最新 Gray 参考程序的平台复测

2026-11-11，直接从用户的 `GSIS_Gray_20261110.zip` 提取源码及原网格，在独立目录编译运行。保持原 Callaway 声子模型、200 三角形、P3、80×40 角度点、`tau_R=tau_N=0.01`、四壁热化温度 `[0,1,0,0]`、`ST=[1,1,1]`。它是热传导方腔，不是上述气体顶盖方腔。macOS 使用压缩包自带的 UMFPACK/PARDISO 接口；限制 OpenMP 至 4 线程、OpenBLAS 为 1 线程。输运迭代及宏观模块公式未改变，源码哈希与诊断补丁见 `results/gray_reference_platform/manifest.json` 和 `diagnostic_only.patch`。

只读诊断分别将宏观回传源、当前 VDF 自身矩代入参考程序已有的完整输运残差例程；临时替换后逐值恢复源变量。另从压缩包编译未加诊断的原始主程序和输运模块，原版按默认 1e-5 在第 8 步停止，与诊断版前 8 步的宏观更新/总量最大差约 1.1e-19，确认诊断没有改变迭代轨迹。为看清平台，诊断版只延后停止到 100 步。

GSIS 第 8 步宏观相对更新 5.0966e-6；第 100 步降至 2.7498e-15。原 kinetic 方程弱残差的绝对范数停在 1.091418e-5，除以初始原方程残差 2.421259e-2 后为 4.507647e-4；耦合方程同一归一化残差为 3.4572e-14。40→100 步原残差平台的相对变化仅约 6.5e-13。宏观源与 VDF 自身矩的组合相对 L2 差异保持 1.8081e-4（0.01808%）。原参考实现也存在耦合固定点与原 kinetic 离散方程不完全一致的现象。

归一化必须区分：参考自带例程用 `sqrt(sum(w*R^2)/sum(w*((A*H)^2+b^2)))`，第 100 步原方程该值为 8.220845e-6；不是上面初态归一化的 4.507647e-4。绝对范数是参考弱残差的速度加权系数范数，不是我们气体代码的逆质量矩阵 L2 范数，因此平台数值不能直接跨模型比较。

辅助纯 CIS 控制只跑到 500 步，尚未达到默认宏观停机阈值；其原/耦合残差一致，源与 VDF 矩的相对差约 8.4e-17。该运行用于核查诊断自洽，不作为收敛 CIS 参考解，也未计算相对收敛 CIS 的剖面误差。

复现及分析：`python3 debug/run_gray_reference_platform.py --cis-max 500 --gsis-steps 100`、`python3 debug/check_gray_native_audit.py`、`python3 debug/analyze_gray_reference_platform.py`。结果为 `results/gray_reference_platform/reference_platform.png`、`summary.csv`、`summary.json`，原始日志及逐步残差分别在 `gsis/`、`cis/`，原版交叉核查在 `native_gsis/`。

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
| `normalQuadrature` | `'halfLegendre'` | 可选 `'hermite'`（允许奇数阶和零法向速度）或参考程序的 `'referenceCubic'`（要求偶数阶） |
| `normalCutoff` | `8` | `halfLegendre` / `referenceCubic` 法向截断；对齐参考三次映射时显式设为 `4` |
| `tangentialQuadrature` | `'hermite'` | 可选 `'referenceCubic'`，将流向速度也改为零附近加密的有限区间求积 |
| `tangentialCutoff` | `4` | 仅用于流向 `referenceCubic`；Hermite 不使用该截断 |
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

### 参考程序的三次映射速度网格（2026-11-09）

`NonOscillatory/Velocity_Mesh.f90` 的 `VELM_FLAG=3` 使用均匀辅助坐标加三次映射；原 `control.in` 选择的正是这个分支。两个速度方向的截断均为 `a=4`、指数 `beta=3`。对每个方向：

```text
xi_i = (2*i - 1 - N)/(N - 1),  i=1,...,N
v_i  = a*xi_i^3
dv_i = 6*a/(N - 1)*xi_i^2
```

`xi_i` 等距，物理速度 `v_i` 在零附近加密；参考代码直接累加所有端点权重，没有再做梯形端点减半。移植到本求解器后，一维权重为 `dv_i*exp(-v_i^2)`，二维权重为两方向一维权重的乘积除以 `pi`，已含 Maxwell 分布。

```matlab
opt = struct('velocityOrder',[32 32], ...
    'normalQuadrature','referenceCubic','normalCutoff',4, ...
    'tangentialQuadrature','referenceCubic','tangentialCutoff',4, ...
    'delta',0.1,'reference',false,'plot',false,'save',false);
c = mono1D2V(opt);       % 分别运行两种流动的 CIS
g = monoGSIS1D2V(opt);   % 相同网格的 GSIS
```

每个方向仍为 32 点，共 1024 点；最小非零速度为 `4/31^3≈0.000134269`，`|v|<0.1` 内有 10 点。旧法向 32 点半区 Legendre 的最小非零速度约为 `0.0423963`，该区间内只有 2 点。默认配置保留旧网格，以上选项显式启用新网格，旧结果不重标为新网格结果。

运行 `addpath('debug'); check_reference_cubic_velocity` 可对照原网格、仅法向三次映射、双方向三次映射。使用相同 Nx64、32×32 速度和相对更新阈值 `1e-5`，双方向三次映射的 CIS 速度相对已保存独立 BGK 数值参考的偏差为：

| δ | Couette：原网格 → 三次映射 | Poiseuille：原网格 → 三次映射 | CIS / GSIS 步数：Couette；Poiseuille |
|---|---|---|---|
| 0.1 | 4.3399% → 0.02646% | 0.63995% → 0.01503% | 5 / 5；8 / 8 |
| 1 | 0.21183% → 0.004149% | 0.15793% → 0.003195% | 11 / 11；23 / 23 |
| 10 | 0.01007% → 0.009548% | 0.03136% → 0.03125% | 82 / 47；251 / 66 |

参考自身仍有数值误差；该对照同时改变了法向节点映射与截断范围，不能视为单独改变截断的研究。三次映射的 δ=0.1 Poiseuille CIS 流量 `Qp≈2.03272032`，相对文献 BGK 表值 `2.0327` 的有符号偏差约 `+0.00100%`。本次原网格回归保持六组结果不变，Maxwell 求积矩、碰撞不变量、零法向流和壁面质量通量检查通过。

迭代次数未改变，速度求积的改善也没有消除耦合误差：δ=10 时 GSIS-CIS 的速度相对差为 Couette `0.08241%`、Poiseuille `0.31906%`。稳定参数、宏观混合、边界和停止准则均沿用原设置。结果保存在被 Git 忽略的 `results/reference_cubic_velocity32/`。

完整流量检查：运行 `addpath('debug'); check_cubic_poiseuille_flowrates`，比较 δ=`[0.1 0.5 1 2 5 10]` 的旧网格和双方向三次映射网格。主表固定 Nx64、32×32 速度、相对更新 `1e-5`，同时记录 CIS/GSIS 流量、文献偏差和步数；另保存三次映射 CIS 的 `1e-9` 停机敏感性检查。三次映射主表为 `results/cubic_poiseuille_flowrates32/cubic_results.csv`，完整对照为同目录 `summary.csv`，严格停机表为 `strict_cis_1e9.csv`；这些生成文件由 Git 忽略。六组新网格 CIS 的最大文献偏差约 `0.08255%`，GSIS 约 `0.3062%`；CIS 收紧至 `1e-9` 后最大偏差约 `0.07258%`，两种停机结果不能混作同一张主表。

六页报告的完整复测：MATLAB 运行 `addpath('debug'); run_cubic_channel_presentation_checks`，覆盖 δ=`[0.1 1 10]`、Couette/Poiseuille、Nx=`[32 64 128]` 的 18 组空间配置，以及 Nx64 下速度 `[32 32]`、`[64 32]`、`[128 32]`、`[128 64]` 的 24 组速度配置（两组集合共享 Nx64 的六个主算例），均运行 CIS/GSIS。保持相对更新 `1e-5`、稳定参数及混合方式不变。随后运行 `python3 debug/plot_channel_presentation.py --cubic` 和 `python3 debug/build_cubic_channel_results_pdf.py`；新图表位于 `output/channel_cubic_figures_20261109/`，当前六页 PDF 为 `output/pdf/1D2V_CIS_GSIS_全指标与停机敏感性_20261109.pdf`，第 6 页包含大小 δ 的停机敏感性和耦合一致性诊断。该工作流保留旧报告与旧数据，不显示独立 BGK 曲线。18 组空间对照的迭代次数与旧网格完全相同；当前 32×32 CIS 相对三次映射 128×64 的速度变化，在 δ=0.1 时为 Couette `0.02682%`、Poiseuille `0.01520%`，δ=1 时为 `0.0006264%`、`0.0004453%`，δ=10 时小于 `0.000007%`。这些是有限网格间差异，不代表相对精确解误差。

停机敏感性：`addpath('debug'); check_cubic_stopping_residuals` 对六组主算例分别测试相对宏观更新阈值 `1e-5/1e-7/1e-9`，每次从 H=0 开始，统一使用各算例 CIS 首轮残差作为分母。δ=10 Poiseuille 的 CIS 原方程相对残差为 `2.969e-4 → 2.997e-6 → 3.026e-8`，GSIS 耦合方程为 `2.999e-4 → 2.787e-6 → 2.997e-8`，但 GSIS 对原 CIS 方程的相对残差为 `4.425e-3 → 4.346e-3 → 4.346e-3`。严格停机时，该原方程残差与混合源相对 kinetic 自身宏观量的源项缺口范数相符。因此收紧停止阈值可以降低迭代残差，却不能消除当前混合源耦合对原 CIS 方程造成的非零残差；不能统称“减小阈值没用”。完整数据在被 Git 忽略的 `results/cubic_stopping_residuals32/summary.csv`。

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

### 与 Su 等 2019 年 BGK Couette 表 1 的剪切量对照（非 SIS 复现）

2026-11-12，比较来源为 W. Su, P. Wang, H. Liu, L. Wu, *Journal of Computational Physics* 378 (2019) 573-590，DOI `10.1016/j.jcp.2018.11.015`，表 1（印刷页 580，PDF 第 8 页）。使用表中每组第一行的 Jiang & Luo 高精度积分方程值，另存 `paper_table1.csv`。完整模型并非同一模型：论文写三维分子速度，当前求解器为二维分子速度；本次仅针对线性等温 BGK Couette 的剪切子问题。令相对 Maxwellian 的扰动 `H=2*vt*phi(x,vn)`，则密度、温度扰动为零，第三维高斯速度可以积分消去，两者均得到 `vn*d_x phi=delta*(U-phi)`、`U=int(phi*exp(-vn^2)/sqrt(pi))dvn`。在壁面和归一化匹配的条件下，只有这里对比的速度及剪切应力量可按该约化作基准对照，不能推广为完整热力学模型、其他流动或二维方腔的等价。本比较不使用 Shakhov 或 Boltzmann 分子势结果。

算法也不同：论文使用专用于 Couette 的 SIS，先用平均剪切应力解标量宏观方程，再按式（25）直接修正分布函数的速度矩和剪切应力矩；当前程序用多宏观量 HDG 子问题，将宏观解回传碰撞源而不改写 H。因此这不是对论文 SIS 的复现，也不是 SIS 与 GSIS 的公平效率对照。本节的迭代步数只属于本程序的 CIS/GSIS 配对运行；三项数值接近文献不能单独证明当前 GSIS 耦合或固定点正确。

间距 1、两壁速度 `[-1,1]`、无外力、完全漫反射、零初值。基准采用 128 个均匀 P1 单元、32×32 三次映射速度点、双向截断 ±4；分别测试单元宏观更新阈值 `1e-5/1e-8`。GSIS 显式设置全 1 稳定参数、`kineticConsistency=false`、全量宏观回传、不修改 H；此算例配置不改变 1D 求解器默认值。文献表 1 使用非均匀空间点、不同法向速度规则和逐点最大相对更新阈值 `1e-10`，故本文迭代次数只在本次 CIS/GSIS 之间比较，不要求复现论文 SIS 的次数。

表 1 三列分别换算为 `uy_wall/(2*Uw)`、`d_x uy_center/(2*Uw)`、`-mean(Pxy)/(2*Uw)`；文献应力定义为 `P12=2*<vx*vy*H>/Uw`。壁面速度由完整边界分布（规定入射加单元出射）取矩，不能用仅弱满足入射条件的单元内端点速度代替；后者另存 `interiorEndpointVelocityOver2`。中心梯度取中心两侧 P1 单元内导数的平均。主要比较量均来自实际 kinetic H，宏观回传源另外记录。相对偏差为 `(本次数值-文献值)/文献值`，误差列取其绝对值，单位 %。

128 单元、32×32 速度、tol=1e-8 的结果如下，三项误差均为 GSIS 相对文献：

| delta | CIS/GSIS 步数 | 壁面速度误差 % | 中心梯度误差 % | 剪切应力误差 % |
|---|---:|---:|---:|---:|
| 0.01 | 5/5 | 0.14299 | 0.88152 | 0.00001951 |
| 0.1 | 7/7 | 0.003088 | 0.061058 | 0.000001844 |
| 1 | 17/13 | 0.004231 | 0.039719 | 0.002130 |
| 10 | 143/14 | 0.025393 | 0.054714 | 0.023975 |
| 100 | 6319/16 | 0.013431 | 0.034297 | 0.014136 |

停止阈值的影响：delta=100 时，CIS 中心梯度误差由 tol=1e-5 下 1.3103% 降至 tol=1e-8 下 0.0004349%；同样改变阈值，GSIS 的梯度误差为 0.034682%/0.034297%，原 kinetic 相对残差为 7.8555e-4/7.8616e-4。这表明 CIS 的早停偏差和未加内部一致性修正 GSIS 的残留偏差不同，不能用统一停机阈值代替原方程检查。

512 单元、64×32 速度、tol=1e-8 的 GSIS，五个 delta 的三项误差最大值为 0.02043%；delta=100 的壁面速度/中心梯度/剪切应力误差为 0.006555%/0.012528%/0.006665%，但原 kinetic 相对残差仍为 2.3847e-4，尚不能称为充分满足原离散 kinetic 方程。单因素核查：delta=0.01 固定 128 单元，仅把法向速度从 32 改为 64，GSIS 中心梯度误差由 0.88152% 降为 0.016481%；delta=100 固定 32×32 速度，把单元数从 128 改为 512，梯度误差由 0.034297% 降为 0.012528%，在 512 单元下继续改为 64×32 几乎不再改变这一项。它们分别显示零速附近积分分辨率与空间/宏观离散的影响。

复现命令（共 29 次实际求解）：

```matlab
addpath('debug');
compare_su2019_bgk_couette('baseline');
compare_su2019_bgk_couette('refinement');
compare_su2019_bgk_couette('isolation');
plot_su2019_bgk_couette;
```

结果位于 `results/su2019_bgk_couette/`：`baseline_summary.csv`、`refinement_summary.csv`、`isolation_summary.csv` 含文献值、实际值、带符号及绝对相对偏差、迭代数与原/耦合残差；各组 MAT 保存完整 H。`provenance.json` 记录定义与原 PDF 校验值。`paper_comparison.png` 比较三项文献相对误差，`stopping_and_residual.png` 分开显示迭代次数和原方程残差；图例 V 表示法向×切向速度点数。没有修改生产求解器公式。

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

2026-11-10 更新：适配器从 H 提取 Rxx/Rxy/RT，使用弱导数 `M*Rx=-G*R+[-Rhat_left;Rhat_right]`。先逐速度构造完整 kinetic 迎风迹，再取高阶矩；包含内部面跳跃和两侧漫反射入射分布。

宏观壁面采用用户提供的 `GSIS_Gray_20261110` 中的仿射通量结构，系数改为当前 BGK 的离散速度积分：

```text
F_HDG(U,Uhat) = F_kinetic(H) + Jwall*(U-Ukinetic)
U = [p ux uy T]
F_kinetic = integral (v*n)*[1, 2*vx, 2*vy, v^2-2]*Htrace*M0 dv
```

`F_kinetic` 包含全部出射与漫反射入射速度。`Jwall` 使用与碰撞项相同的离散 Gram 平衡投影，对出射分布作宏观线性化，同时调整入射壁面密度以保持零法向质量通量；壁速和壁温固定。所有守恒方程保留，增加体平均密度扰动为零的约束，消除闭合通道的零空间。宏观面变量由通量方程求出，不直接指定为 kinetic 单元矩。

仅迁移上述壁面和高阶导数仍会留下内部 HDG 稳定通量与 kinetic 迎风通量的离散差异。当前另外使用 `kinetic_hdg_consistency_1d2v` 冻结七个矩的弱载荷修正：在当前 H 对应的宏观状态处，四个守恒残差等于原 kinetic-DG 残差的矩，三个应力/热流辅助残差等于对应矩的负值。宏观算子及稳定参数保持不变。因而原 kinetic 残差为零时，宏观修正也为零。这是 BGK 适配时增加的一致性处理，不是对 Gray 源码的逐项照搬。

第一轮 kinetic sweep 后，适配器可直接调用：

```matlab
[macro,aux,hdg,hot] = gsis_macro_from_kinetic_1d2v(cis.poiseuille);
[macro,aux,hdg,hot] = gsis_macro_from_kinetic_1d2v( ...
    cis.poiseuille,struct('cache',hdg.cache));
```

完整 GSIS 每轮执行：由独立保存的 Qsource 组装碰撞增益，解 kinetic H，取 Qkin，解宏观方程，再更新下轮源项：

```text
Qsource = Qmacro
```

宏观步骤不直接修正 H。保留原有闭合通道的密度均值为零规范，并对 Qsource 的密度均值作同样约束。`r.nodalMoments/cellMoments/flowRate` 始终来自实际 kinetic H；`r.synthetic` 另外保存宏观解、辅助量、面变量、高阶项、兼容字段 `omega=1` 和 `sourceMoments`。自定义稳定参数可通过 `options.macroParam.stabilization` 传入。

### 收敛诊断与验证

默认停机与参考 Fortran 的 `Calculate_Residual_rho_T_U` 对齐。每个单元计算 `R_i=int(1+rho_pert)dx`、`T_i=int(1+T_pert)dx`、`V_i=hypot(int(ux)dx,int(uy)dx)`；对这三个单元向量分别取 `norm(A_new-A_old)/norm(A_new)`，最大值小于 `1e-5` 时停止。CIS 使用 Qkin，GSIS 使用回传的 Qsource。首次 OLD 取静止参考气体：密度/温度为单元长度，速度为零。零分母且无变化记为零，有变化记为 Inf。实现见 `reference_macro_update_1d2v.m`。

GSIS 的 `diagnostics.coupled` 是将当前回传源 Qsource 代入 kinetic DG 方程得到的残差；`diagnostics.kinetic` 则使用 Qkin 自身构造碰撞源，表示原 CIS 方程残差。默认模式下这些残差及壁面诊断只记录，不参与停机；因此宏观量更新小于 tol 不代表方程残差也小于 tol。`stopCriterion='equation'` 恢复原来的 kinetic 矩更新、源项更新、耦合残差、壁面残差、壁面质量通量和密度均值误差联合判据。历史文件保留原有列，在最后追加四列参考相对更新指标；`r.convergence` 记录所用判据、阈值和最终指标。

`runGSISChannelComparison` 的结果表记录完整外层迭代次数、两种残差、速度剖面和流量相对 CIS 的差异、实际稳定参数，以及停机判据、阈值和最终相对更新。相同宏观量停机标准并不保证最终离散解完全相同，也不表示相同运行时间。`demo_gsis_macro_1d2v` 仍是 CIS 收敛后仅执行一次宏观求解的子模块示例，其步数不是完整 GSIS 步数。

2026-11-10 当前耦合实测：64 个 P1 空间单元，双向 32 点 referenceCubic（截断 ±4、零附近加密），零初值，Couette 壁速 ±0.01、Poiseuille 加速度 0.01，稳定参数 `[1 delta delta delta]`。相同参考宏观更新阈值 `1e-5`：

| 流动 | delta | CIS 步数 | GSIS 步数 | kinetic 速度相对 CIS 的 L2 差异 |
|---|---:|---:|---:|---:|
| Couette | 0.1 | 5 | 5 | 0.00000219% |
| Poiseuille | 0.1 | 8 | 5 | 0.0000599% |
| Couette | 1 | 11 | 9 | 0.0002089% |
| Poiseuille | 1 | 23 | 8 | 0.001491% |
| Couette | 10 | 82 | 8 | 0.008740% |
| Poiseuille | 10 | 251 | 14 | 0.03022% |

阈值收紧至 `1e-11` 时，delta=10 的 CIS/GSIS 步数分别为 Couette 204/20、Poiseuille 678/30；GSIS 原 kinetic 方程绝对残差分别为 `1.3595e-13`、`1.5716e-13`，除以各自共同首步 CIS 残差后为 `9.3374e-12`、`1.1728e-11`。对应 kinetic 速度差异为 `7.9167e-9%`、`2.9895e-8%`。在所测阈值范围内，未再观察到原有残差平台。这里比较同一离散方程的固定点和迭代步数，不代表已获得文献或连续解精度、运行时间加速比或网格无关性。

复现：`addpath('debug'); check_gray_gsis_coupling_1d2v`。原始数据与完整历史位于 `results/gray_consistent_coupling32/summary.csv`、`comparison.mat`；收敛 CIS 后单步宏观检查位于 `converged_cis_one_macro.csv`。`check_gray_gsis_consistency_1d2v` 独立检查漫反射 Jacobian、完整 kinetic 残差的七个矩、以及带壁面和密度约束的整体线性系统。

关闭内部弱残差修正的对照：`check_gray_gsis_coupling_1d2v(false)`，结果单独保存于 `results/gray_uncorrected_coupling32`，使用相同的壁面通量、弱高阶导数、完整宏观回传及稳定参数。单次调用可设置 `options.macroParam.kineticConsistency=false`；默认仍为 true。脚本核查关闭时宏观右端只含原有体源及高阶项，额外弱载荷严格为零。

关闭时，1e-5 下 Couette 的 CIS/GSIS 步数为 5/5、11/9、82/8，速度差异为 0.000161%、0.01479%、0.1007%；Poiseuille 为 8/5、23/8、251/13，速度差异为 0.02178%、0.1177%、0.3165%（delta 依次为 0.1/1/10）。收紧到 1e-11 后，delta=10 的 Couette/Poiseuille 速度差异仍为 0.10848%/0.28846%，原 kinetic 方程绝对残差为 2.4570e-5/1.8028e-4，而耦合残差为 6.8869e-13/2.0398e-12。因而壁面迁移本身能减少迭代步数，但没有消除原离散方程的残差平台。

稳定参数对照：`check_gray_gsis_coupling_1d2v(false,[1 1 1 1])` 仅在本次测试中覆盖惩罚系数，保持求解器默认值不变，结果保存于 `results/gray_uncorrected_coupling32_S1_1_1_1`。当前切向动量通量为 `-Ly/delta*n+s_uy*(uy-uy_hat)`，代入 `sigma_xy=-Ly/delta` 后，惩罚项仍为 `s_uy*(uy-uy_hat)`，因此默认 `s_uy=delta` 没有被辅助变量缩放抵消。全 1、关闭内部修正、阈值 1e-11 时，delta=10 Couette/Poiseuille 的 GSIS 步数为 20/28，速度相对 CIS 差异为 0.08779%/0.22934%，原 kinetic 绝对残差为 2.3630e-5/1.7342e-4。相比 `[1 delta delta delta]` 偏差减小，但残差平台仍存在。

下列为 2026-11-05 的历史实现结果（旧内部导数、单元矩壁面和混合回传），不代表当前耦合：

| 流动 | delta | CIS 步数 | GSIS 步数 | GSIS 单元平均速度相对 CIS 的 L2 差异 |
|---|---:|---:|---:|---:|
| Couette | 1 | 11 | 11 | 0.001515% |
| Poiseuille | 1 | 23 | 23 | 0.009321% |
| Couette | 10 | 82 | 47 | 0.07736% |
| Poiseuille | 10 | 251 | 66 | 0.3003% |

完整结果保存于 `results/gsis_channel_comparison/reference_stopping/summary.csv`。相同停机条件下，CIS 后只做一次 HDG 的原始速度差异依次为 0.1890%、0.3441%、0.09709%、0.2395%，保存在 `results/gsis_macro_1d2v/converged_cis_one_step/reference_stopping/summary.csv`。停机对齐没有消除宏观解与 CIS 的差异。

```matlab
addpath('debug');
report_hdg_after_converged_cis_1d2v % CIS 收敛后只解一次当前 HDG，直接比较
check_hdg_macro_1d2v          % 制造解、缓存、混合边界、解析流、二阶精度
report = audit_hdg_weak_form_1d2v; % 独立积分组装、通量、界面公式、2V 系数
check_gsis_channel_1d2v      % 不修改 H、下轮使用完整宏观源、平移平衡态、停机检查
check_gray_gsis_consistency_1d2v % 壁面 Jacobian、七矩残差、宏观整体组装
check_gray_gsis_coupling_1d2v % 0.1/1/10，三档停机阈值，完整 CIS/GSIS 对比
check_channel_1d2v           % 原 CIS 对整体 kinetic DG 解、物理和输出回归
check_reference_stopping_1d2v % 参考宏观相对更新、源项监测、零分母与默认阈值
```

独立核查的参考矩阵由三点 Gauss 积分逐基向量生成，不调用生产 F/G/M 或 Kronecker 装配。涵盖非均匀网格、delta=0.07/1/7/1000、默认及自定义稳定参数、法向/切向驱动、界面平均/跳跃公式和二维 Gaussian 恒等式。当前 `fixed_point_audit.csv` 将 kinetic 解及通量相容的基准面变量代入带一致性载荷的宏观弱式，独立核查七矩残差匹配；不要求内部 HDG 通量逐面等于 kinetic 通量。`audit_macro_wall_reference_1d2v`、`audit_hdg_bulk_interface_1d2v`、`check_hdg_one_step_refinement_1d2v` 用于历史单元矩壁面/内部高阶导数数据分析，未迁移至当前仿射通量边界。

以下保留 2026-11-03 的严格方程停机历史结果：64 空间单元、32×32 速度点、tol=1e-10。复现旧数值须使用当时的实现，仅调整当前代码的停机参数不会恢复旧结果。

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
