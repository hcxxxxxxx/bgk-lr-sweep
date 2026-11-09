# 测试算例设置

历史主结果核对日期：2026-10-29。通道测试包含 Couette 与 Poiseuille 各两组，碰撞参数分别为 delta=1 和 delta=10。下表列出四组已保存主结果的实际设置；直接运行的默认值、加密验证和其他求解器在后文列出。

2026-11-05 起，CIS/GSIS 默认统一为参考代码的单元宏观量相对更新停机，`stopCriterion='reference'`、`tol=1e-5`。下表和 JSON 保留历史快照；复现历史结果须显式使用 `stopCriterion='equation'`。新结果见文末。

主结果配置取自各自 `channel.mat` 的 `r.config`。完整配置和结果标识另存于 [CASE_SETTINGS.json](CASE_SETTINGS.json)，该文件是本次核对的快照，不是求解器自动读取的配置文件。

## 四组通道主算例

| 设置 | Couette delta 1 | Poiseuille delta 1 | Couette delta 10 | Poiseuille delta 10 |
|---|---|---|---|---|
| 主程序 | `mono1D2V` | `mono1D2V` | `mono1D2V` | `mono1D2V` |
| `flow` | `'Couette'` | `'Poiseuille'` | `'Couette'` | `'Poiseuille'` |
| `delta` | 1 | 1 | 10 | 10 |
| `wallVelocity`，左壁、右壁的 y 向速度 | `[-0.01 0.01]` | `[0 0]` | `[-0.01 0.01]` | `[0 0]` |
| `acceleration`，沿 y | 0 | 0.01 | 0 | 0.01 |
| `gap` | 1 | 1 | 1 | 1 |
| `cells`，物理单元数 | 64 | 64 | 128 | 128 |
| `velocityOrder`，`[Nvx Nvy]` | `[128 8]` | `[128 8]` | `[128 8]` | `[128 8]` |
| `normalQuadrature` | `'halfLegendre'` | `'halfLegendre'` | `'halfLegendre'` | `'halfLegendre'` |
| `reference` | `true` | `true` | `true` | `true` |
| `referenceCells` | 1024 | 1024 | 1024 | 2048 |
| `referenceOrder`，正半轴速度点数 | 160 | 160 | 160 | 320 |
| `initial` | 0 | 0 | 0 | 0 |
| `maxIter` | 30000 | 30000 | 30000 | 30000 |
| `tol`，绝对阈值 | `1e-10` | `1e-10` | `1e-10` | `1e-10` |
| `checkEvery` | 10 | 10 | 10 | 10 |
| `reportEvery` | 100 | 100 | 100 | 200 |
| `plot`，主求解监视图 | `true` | `true` | `false` | `true` |
| `save` | `true` | `true` | `true` | `true` |
| 分布数组 `H` 大小 | `128×1024` | `128×1024` | `256×1024` | `256×1024` |
| 保存结果的收敛步数 | 16 | 39 | 160 | 579 |

Couette delta=10 的主求解设置 `plot=false`，物理诊断图由 `plotCouette1D2V` 单独生成。它和是否保存数值结果是两个独立开关。

四组输出目录分别是：

| 算例 | 输出目录 |
|---|---|
| Couette delta=1 | `results/channel_1d2v/couette` |
| Poiseuille delta=1 | `results/channel_1d2v/poiseuille` |
| Couette delta=10 | `results/channel_1d2v/delta10/couette` |
| Poiseuille delta=10 | `results/channel_1d2v/delta10/poiseuille` |

`outputDir` 指向上表目录的父目录，程序再按 `flow` 追加 `couette` 或 `poiseuille`。同一目录下再次运行会覆盖该算例的结果文件。

## 通道共同物理设置

采用稳态、二维速度的线性化 BGK 模型：

```text
M0 = exp(-vx²-vy²)/pi
f = M0*(1+H)
vx*dH/dx = delta*(P(H)-H) + 2*a*vy
P(H) = rho_pert + 2*vx*ux + 2*vy*uy + T_pert*(vx²+vy²-1)
```

| 项目 | 设置 |
|---|---|
| 物理区域 | `x ∈ [-0.5,0.5]`，只离散 x |
| 坐标方向 | x 为壁面法向，y 为流向；主要剖面为 `uy(x)` |
| 参考状态 | 实际密度、温度均为 1，平均速度为零 |
| 初始分布 | `H=0`，即完整分布 `f=M0` |
| 左右壁面温度 | 实际温度均为 1；温度扰动均为 0 |
| 壁面法向速度 | 0，壁面间距固定 |
| 壁面反射 | 完全漫反射，漫反射比例为 1 |
| 入射半空间 | 左壁 `vx>0`，右壁 `vx<0` |
| 入射扰动 | `H_in = n_w + 2*vy*Uwall` |
| 壁面密度扰动 | 从出射分布按离散零质量通量求得；当前纯剪切/力驱动等温线性化解理论上 `n_w=0` |
| 总质量约束 | 域平均密度为 1，通过常数密度扰动平移消除自由度 |
| 加速度源项 | `G=2*a*vy`，只在 Poiseuille 中非零 |

`rho_pert`、`T_pert` 为扰动量：实际密度和线性化温度为 `1+rho_pert`、`1+T_pert`。本模型采用二维速度温度基函数；不含第三速度方向的独立约化能量分布、黏性加热或 Shakhov 修正。上表保存结果为 CIS；新增的 `monoGSIS1D2V` 入口使用同一 kinetic 模型与 HDG 宏观方程耦合，见下文。

delta=1、10 是碰撞强弱对照的试算值，壁速 0.01 和加速度 0.01 用于小扰动试算。新增 GSIS 按参考代码取 `Kn_cell=sqrt(pi)/(2*delta*dx)` 计算混合系数；这些数值尚未对应到某组有量纲气体和壁面参数。

## 通道共同数值设置

| 项目 | 设置 |
|---|---|
| 物理网格 | 均匀网格，无需 Fluent 文件 |
| 空间离散 | P1 迎风 DG，每单元、每速度 2 个独立自由度 |
| 主算例法向速度 | 128 点，其中正负半轴各 64 点，范围 `[-8,0]`、`[0,8]` |
| 法向求积 | 各半轴 Gauss–Legendre，权重包含法向 Gaussian |
| 流向速度 | 8 点 Gauss–Hermite |
| 总速度点数 | `128*8=1024` |
| 求矩 | 权重已包含 `M0`，不重复乘 Maxwell 分布 |
| 碰撞离散 | Gram 矩阵投影，离散守恒质量、两方向动量、能量 |
| 主迭代 | 稳态源迭代，碰撞增益和壁面密度滞后一轮 |
| 输运求解 | 每个 `vx` 共用一个预分解的块三角 DG 矩阵，多个 `vy` 右端同时求解 |
| 历史主结果停止条件 | 矩更新、当前动理学方程残差、壁面关系误差、壁面质量通量、平均密度误差同时低于 `tol`；现在用 `stopCriterion='equation'` 复现 |

`cells` 是单元数；DG 节点自由度数为 `2*cells`。`velocityOrder(1)` 包含正负全部法向速度点，而 `referenceOrder` 只统计参考解的正半轴速度点，两种计数方式不同。

独立 BGK 参考解由 `channel1D2VReference` 计算：利用 `H=2*vy*phi(x,vx)`，对分片常数的宏观速度沿特征线积分，直接解积分方程。正速度区间为 `[0,8]`，负速度由对称性处理；参考解只用于比较，不参与 DG 迭代。比较值均按 DG 单元平均速度对齐。

NS 无滑移曲线使用 `mu=1/(2*delta)`。静止壁面的 Poiseuille 形式为 `uy_NS=a*delta*(0.25-x²)`。该曲线是连续极限对照，同模型 BGK 数值参考才用于当前稀薄流的误差核对。

## 默认运行与主结果复现

在项目目录直接执行：

```matlab
result = mono1D2V;
```

这会计算两组 delta=1 的当前默认算例：64 单元、32×32 速度点、参考宏观量停机、tol=1e-5，与历史表格的速度网格和停机设置不同。单独写 `mono1D2V(struct('flow','Couette','delta',10))` 时，`cells` 仍默认 64；复现历史结果请使用下列完整配置。

下面的四组调用按表中设置复现。它们写入现有结果目录；只想试算时将 `save` 设为 `false`。

```matlab
root = pwd;  % MATLAB 当前目录应为项目目录
cfg = struct('gap',1,'cells',64,'delta',1, ...
    'velocityOrder',[128 8],'normalQuadrature','halfLegendre', ...
    'reference',true,'referenceCells',1024,'referenceOrder',160, ...
    'initial',0,'maxIter',30000,'tol',1e-10,'stopCriterion','equation', ...
    'checkEvery',10,'reportEvery',100,'plot',true,'save',true, ...
    'outputDir',fullfile(root,'results','channel_1d2v'));

cfg.flow = 'Couette'; cfg.wallVelocity = [-0.01 0.01]; cfg.acceleration = 0;
c1 = mono1D2V(cfg);

cfg.flow = 'Poiseuille'; cfg.wallVelocity = [0 0]; cfg.acceleration = 0.01;
p1 = mono1D2V(cfg);

cfg.cells = 128; cfg.delta = 10;
cfg.outputDir = fullfile(root,'results','channel_1d2v','delta10');
cfg.flow = 'Couette'; cfg.wallVelocity = [-0.01 0.01]; cfg.acceleration = 0;
cfg.plot = false;
c10 = mono1D2V(cfg);
plotCouette1D2V(c10);

cfg.flow = 'Poiseuille'; cfg.wallVelocity = [0 0]; cfg.acceleration = 0.01;
cfg.referenceCells = 2048; cfg.referenceOrder = 320;
cfg.reportEvery = 200; cfg.plot = true;
p10 = mono1D2V(cfg);
```

需要完整的 Poiseuille delta 对照和加密检查时执行 `r=runPoiseuilleDelta10`。它会额外计算一组 128 单元的 delta=1 基线，以及 delta=10 加密解，配置如下表。

## 验证与加密配置

| 入口 | 设置变化 | 用途和保存位置 |
|---|---|---|
| `check_channel_1d2v` | 4 单元，`[6 4]` 整区间 Hermite，`tol=1e-12` | 与独立整体 DG 线性系统对照，`save=false` |
| 同一入口的平移平衡态 | 4 单元，`[7 5]` Hermite，两壁 `Uy=0.007`，以 `H=2*0.007*vy` 初始化 | 验证零法向速度点及匀速平衡保持 |
| 同一入口的对称性测试 | 32 单元，`[32 6]` Hermite，`tol=1e-12` | Couette 反对称、Poiseuille 对称与驱动反号 |
| 同一入口的无碰撞极限 | Couette，`delta=1e-9`，4 单元，`[32 6]` Hermite | 半空间分布及应力极限 |
| 同一入口的停止判据测试 | Poiseuille，`delta=100`，4 单元，`[6 4]` Hermite，`maxIter=1`，`tol=1e-4` | 故意停止并要求 `converged=false`；不是已收敛的 delta=100 物理算例 |
| `check_channel_reference` 的旧求积对照 | 两种流，delta=1，64 单元，`[32 8]` 整区间 Hermite | 复现旧速度求积误差，`save=false` |
| 同一入口的 DG 加密 | 两种流，delta=1，128 单元，`[256 8]` 半区混合求积 | 参考解仍为 1024 单元、正半轴 160 点，`save=false` |
| 同一入口的参考解加密 | 参考网格 512→1024→2048；正半轴速度阶数 160→320 | 独立参考解的自身收敛检查 |
| `check_couette_physics` | 读取默认 Couette；另算 128 单元的 delta=10 Couette | 保存 delta=10 主结果及两组物理诊断图 |
| `runPoiseuilleDelta10` 的对照基线 | delta=1，128 单元，`[128 8]`；参考为 2048 单元、正半轴 320 点，`reportEvery=200` | 重新计算的基线，`plot=false, save=false`；与默认已存 64 单元的 Poiseuille 分开 |
| 同一入口的加密解 | delta=10，256 单元，`[256 8]`，`reference=false` | 与主解使用的 2048/320 参考解比较，`plot=false, save=false` |
| 同一入口的参考网格核对 | delta=10，参考网格 1024→2048，正半轴均 320 点 | 检查参考解的空间离散误差 |

1D2V 验证入口：

```matlab
addpath('debug');
check_channel_1d2v
check_channel_reference
check_couette_physics
runPoiseuilleDelta10
```

重启读写测试位于 `tmp`。这些测试会把 `initial` 设为已有分布并仅推进一步，属于结果读写和保持性检查；上表四个物理主算例均从 `initial=0` 开始。

## 已保存结果核对

| 算例 | 是否收敛 | 步数 | 动理学残差 | 相对同模型数值参考的速度 L2 差异 |
|---|---|---:|---:|---:|
| Couette delta=1 | 是 | 16 | `8.06e-11` | `3.73e-5` |
| Poiseuille delta=1 | 是 | 39 | `6.18e-11` | `2.67e-5` |
| Couette delta=10 | 是 | 160 | `9.61e-11` | `2.04e-5` |
| Poiseuille delta=10 | 是 | 579 | `9.76e-11` | `2.48e-5` |

L2 数值是相对差异，例如 `2.48e-5` 对应约 `0.00248%`；它包含被比较的独立参考解自身的离散误差。当前记录支持上述设置下的验证，改变 delta、驱动或求积后应重新检查。

## 其他求解器的设置

项目还保留两个不同的模型入口，边界和归一化不与上述通道算例混用。

| 项目 | 2D2V BGK 温差条带 | Ma 3 正激波 |
|---|---|---|
| 主入口 | `monoBGKDG` | `monoCIS`、`monoNS` |
| 模型 | 二维速度线性化 BGK | 三维单原子气体的二维速度约化 Shakhov / NS |
| 物理区域 | `x=[0,100], y=[0,0.25]` | `x=[-60,60], y=[0,1]` |
| 网格 | `mesh/nsw400.cas`，400 个矩形单元 | 同一个文件，缩放坐标，400 单元 |
| 空间离散 | Q1-DG，每单元每速度 4 个自由度 | 一阶有限体积 |
| delta | 0.01 | `sqrt(pi/2)`；碰撞频率还依赖密度、温度 |
| 速度求积 | 16×16 Gauss–Hermite | vx 为 `[-18,18]` 的 96 点，vy 为 `[-15,15]` 的 80 点，等权矩形求积 |
| 驱动及初始态 | 左右壁温 `[0.95 1.05]`，体源和初始扰动均为 0 | Ma=3，上游 rho=T=1，下游按 Rankine–Hugoniot；初始为中心两侧不连续平衡态 |
| 边界 | 左右漫反射；上下镜面反射 | 左右远场入射分布；上下镜面反射 |
| 主要迭代参数 | 最多 30000 步，矩/方程/壁面阈值均 `1e-9`，`checkEvery=20` | CIS 最多 6000 步、更新阈值 `1e-8`；纯 NS 最多 4000 步、更新阈值 `1e-8`、CFL=5 |
| GSIS | 不使用 | `demo_method='GSIS'; monoCIS`，前 10 步 CIS，第 10 步初始化宏观状态，此后每轮 100 步 NS；外层阈值 `1e-6` |

`monoBGKDG` 当前源码默认 `wallModel='nonlinear'`、`reportEvery=200`、`plot=true`。现存 `results/bgk_dg/bgk_dg.mat` 的配置则记录 `reportEvery=500`、`plot=false`，且没有 `wallModel` 字段；该字段缺失不能作为保存配置的显式取值。其保存结果为 6989 步收敛。

正激波表中为当前源码设置；现有 `results` 目录只有通道和 BGK 条带结果，不能将 README 中旧激波记录视为本次重新生成的结果。各算法的详细说明仍见 [README](README.md) 和 [求解流程指南](SOLVER_GUIDE.md)。

## 新增参考代码式 1D2V GSIS（2026-11-03）

入口为 `monoGSIS1D2V`，当前默认 64 空间单元、32×32 速度点、delta=1、`stopCriterion='reference'`、tol=1e-5。Couette/Poiseuille 的 kinetic 驱动、漫反射壁面和速度求积沿用 `mono1D2V`。宏观采用压力变量 `p=rho+T`、两倍尺度动量方程、单元内部高阶项导数，以及单元内壁面矩（法向速度置零）。稳定参数仍为 `[1 delta delta delta]`。

每轮只将 `Qkin+omega*(Qmacro-Qkin)` 放入下一轮 kinetic 碰撞源，`omega=min(1,1/Kn_cell)`；宏观步不修改分布 H。停机对齐参考代码：计算每单元完整密度积分、完整温度积分、速度积分的模，分别取新旧两轮差的 L2 范数除以当前量的 L2 范数，三者最大值小于 1e-5 时停止。CIS 监测 kinetic 矩，GSIS 监测混合源矩。耦合残差、原 CIS 残差和壁面误差继续记录，不参与默认停机。结果剖面和流量始终来自 kinetic H。

```matlab
r = monoGSIS1D2V(struct('delta',10,'plot',false,'save',false));
[summary,runs] = runGSISChannelComparison;
```

对照脚本默认使用 64 单元，分别计算 delta=1/10 两个流的 CIS 和 GSIS，参考停机结果位于 `results/gsis_channel_comparison/reference_stopping`。2026-11-05 实测 CIS/GSIS 步数依次为：Couette delta=1：11/11，Poiseuille delta=1：23/23，Couette delta=10：82/47，Poiseuille delta=10：251/66。

`results/gsis_channel_comparison` 根目录保留旧的严格方程停机对照和 `refinement.csv`，不与新的默认结果混用。细节和验证入口见 [README](README.md)。
