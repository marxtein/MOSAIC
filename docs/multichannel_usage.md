# 多通道天线脚手架 使用说明

> 分支: `feature/multichannel-antenna-scaffolding`  ·  代码位置: `src/`
> 阅读前置: 知道 `runme.m` / `run_paper.m` 的常规跑法。

---

## 0. 我能用它做什么

1. **保持现状不变** —— 旧单频天线工作流（论文复现）走原 `ant_interface` + `poincare_perturb` + `cal_char_freq_perturb`，**什么都不用改**。
2. **多频天线扩展** —— 用新的 `ant.channels(k)` 结构和 `eval_perturb_field` 多通道求值器，把任意数量的相干扰动 `(omega_k, n_k, m_modes_k, profiles_k)` 同时叠加到导心轨道。
3. **守恒量诊断** —— 每步逐通道算 `E'_k = E - (omega_k/n_k)·q·P_zeta`，看 multi-channel 下守恒破坏程度。

---

## 1. 核心概念

### 通道（channel）

一个通道 = 一组相干扰动，由以下字段唯一描述：

| 字段 | 含义 |
|---|---|
| `omega`         | 角频率（MOSAIC 归一化单位 = 2π·f[Hz]·`unit.gtc_utime`） |
| `n`             | 环向模数 |
| `m_modes`       | 该通道的极向模数数组 `[m_min ... m_max]` |
| `num_modes`     | `numel(m_modes)` |
| `psi`           | 径向网格（列向量） |
| `phi_real{i}` / `phi_imag{i}`     | 第 i 个 m 谐波的标量势 δφ 的实/虚部，逐 m 元胞 |
| `apara_real{i}` / `apara_imag{i}` | δA‖ 的实/虚部 |
| `dptdp_real{i}` / `dptdp_imag{i}` | ∂δφ/∂ψ 的实/虚部 |
| `dapdp_real{i}` / `dapdp_imag{i}` | ∂δA‖/∂ψ 的实/虚部 |

### 多通道 `ant.channels`

`ant.channels(1..N)` 是一个 struct 数组。`numel(ant.channels)==1` 时所有新代码与旧内联 m 循环 **逐位等价**。`>1` 时叠加。

顶层 `ant.omega`、`ant.n`、`ant.m_modes`、`ant.phi_real` 等字段 **保留** = `ant.channels(1)` 的镜像；只读旧字段的调用者不受影响。

---

## 2. 文件清单 + 用途

| 文件 | 角色 | 何时调用 |
|---|---|---|
| `ant_make_channel.m`        | 单通道打包                     | 你自己构造通道时 |
| `ant_attach_channel.m`      | 把 legacy `ant.*` 包成 `channels(1)` | 跑完 `ant_interface` 想升级时 |
| `ant_merge_channels.m`      | 多通道拼接                     | 把多个单通道 ant 拼成 multi-channel ant |
| `eval_perturb_field.m`      | 多通道场求值器                  | `push_orbit_multichannel` 内部调用，也可手动调 |
| `push_orbit_multichannel.m` | RK2 多通道扰动推进器             | 跑多通道轨道时 |
| `test_synthetic_channel.m`  | 解析 Gaussian 单 m 通道构造器     | 测试用 |
| `test_multi_n_resonance.m`  | **Test A** 驱动                  | 验证无天线多 n 共振线叠加 |
| `test_multichannel_conservation.m` | **Test B** 驱动           | 验证多通道 E'_k 守恒 |

---

## 3. Quick start

### 3a. 一行不改跑论文复现

```matlab
% 原 runme.m Section 4 不动
ps_option   = 'Poincare';
ant_option  = 'on';
PoinE       = [26 34.5 30];
theta_rand  = 1;
run run_paper.m;
```

新脚手架对此路径**零影响**。

### 3b. 拿一份单通道 ant，升级到 multi-channel 接口

```matlab
% Step 1: 跑原 ant_interface 拿 legacy ant
ant = ant_interface('RSAE_ant.mat', 1, -13, -11);
ant.omega = 0.0027;         % MOSAIC 归一化（82.03 kHz）
ant.n     = 4;              % 如 RSAE_ant.mat 不带 n，手动补

% Step 2: 升级为 channels（不破坏 legacy 字段）
ant = ant_attach_channel(ant);

% 现在 ant.channels(1) 可用；ant.omega/.n/.m_modes 仍在
disp(numel(ant.channels));   % => 1
```

### 3c. 用 eval_perturb_field 求场（手动）

```matlab
% 在轨道里某一步，给定 (psi, theta, zeta, t)
amp_mod = 1/20;
t_now   = tstep * (istep + 0.5*irk);   % 与 poincare_perturb 约定一致
F = eval_perturb_field(ant.channels, psi, theta, zeta, t_now, amp_mod);

% 总场（替代原内联 m 循环结果）
dptdp    = F.dptdp;
dptdt    = F.dptdt;
dptdz    = F.dptdz;
apara    = F.apara;
dapdp    = F.dapdp;
dapdt    = F.dapdt;
dapdz    = F.dapdz;
paparapt = F.paparapt;

% 逐通道分量（诊断守恒量用）
disp(F.per_channel(1).dptdp);   % channel-1 的 dptdp 贡献
```

---

## 4. 构造多通道天线 — 三种方式

### 方式 A: 从多个 `.mat` 文件各加载一通道

```matlab
% 通道 1: n=3 RSAE
ant_n3 = ant_interface('RSAE_n3.mat', 1, -10, -8);
ant_n3.omega = 2*pi*59460*unit.gtc_utime;
ant_n3.n     = 3;
ant_n3 = ant_attach_channel(ant_n3);

% 通道 2: n=4 RSAE
ant_n4 = ant_interface('RSAE_n4.mat', 1, -13, -11);
ant_n4.omega = 2*pi*66070*unit.gtc_utime;
ant_n4.n     = 4;
ant_n4 = ant_attach_channel(ant_n4);

% 通道 3: n=6 TAE
ant_n6 = ant_interface('TAE_n6.mat', 1, -18, -15);
ant_n6.omega = 2*pi*97500*unit.gtc_utime;
ant_n6.n     = 6;
ant_n6 = ant_attach_channel(ant_n6);

% 合并
ant = ant_merge_channels(ant_n3, ant_n4, ant_n6);
disp(ant.num_channels);   % => 3
```

### 方式 B: 解析合成（测试/调试用）

参见 `test_synthetic_channel.m` 和 `test_multichannel_conservation.m`：

```matlab
psi_grid  = linspace(0, eq.psiw, 200).';
ch1 = test_synthetic_channel(3,  59.46, 5,  psi_grid, ...
                             0.55*eq.psiw, 0.10*eq.psiw, 1e-5);
ch2 = test_synthetic_channel(4,  66.07, 7,  psi_grid, ...
                             0.55*eq.psiw, 0.10*eq.psiw, 1e-5);
ch3 = test_synthetic_channel(6,  97.50, 11, psi_grid, ...
                             0.55*eq.psiw, 0.10*eq.psiw, 1e-5);

ant_wrap = @(c) struct('channels', c);
ant = ant_merge_channels(ant_wrap(ch1), ant_wrap(ch2), ant_wrap(ch3));
```

### 方式 C: 手动 `ant_make_channel`

完全自己控制 profiles。详见 `test_synthetic_channel.m:33-42` 写法。

---

## 5. 多通道轨道推进

```matlab
% 准备 zpart 初值（14 x np 矩阵；按 poincare_perturb 约定）
np = 6;
zpart = zeros(14, np);
zpart(1, :) = [...psi];
zpart(2, :) = [...theta];
zpart(3, :) = [...zeta];
% slot 4 = rho_||, slot 6 = sqrt(mu)
% 见 test_multichannel_conservation.m:75-95 完整初始化

% 推进
B = push_orbit_multichannel(zpart, ...
        4000, ...      % mstep
        10, ...        % tstep (omega_cp^-1)
        particle.charge, particle.mass, ...
        1, 1, ...      % eq_curr, eq_gradB
        ant, ...       % multi-channel ant
        1.0);          % amp_mod

% 输出
B.t        % (mstep x 1)
B.E        % (mstep x np)         total energy, MOSAIC norm
B.Pzeta    % (mstep x np)         P_zeta
B.Eprime   % (mstep x np x nc)    per-channel E'_k(t)
```

### 推进结果 → keV

```matlab
global unit
E_keV       = B.E       * unit.energy_norm / 1000;
Eprime_keV  = B.Eprime  * unit.energy_norm / 1000;
```

---

## 6. 跑测试

### Test A — 无天线多 n 共振线叠加

```bash
matlab -sd src -batch "test_multi_n_resonance"
```
输出: `output/test_multi_n_resonance.{fig,png}`

### Test B — 三通道 E'_k 守恒诊断

```bash
matlab -sd src -batch "test_multichannel_conservation"
```
输出: `output/test_multichannel_conservation.{fig,png}` + 命令行 `E'_k` 最大漂移表。

### 调振幅

`test_multichannel_conservation.m` 顶部:
```matlab
amp_phi = 2e-5;
```
经验值：
- `1e-7`：等价无扰动，看 RK2 积分本身的数值噪声
- `2e-5`：扰动可见但全程 bounded（推荐起点）
- `5e-4`：合成 Gaussian 包络下轨道发散，**不要用**（梯度过大）
- 真实本征模 `phim_x` 可设回 `1e-3` ~ `1e-1` 量级

---

## 7. 守恒量诊断

每个通道 `k` 单独有
```
E'_k(t) = E(t) - (omega_k / n_k) * q * P_zeta(t)
```

- **单通道开启时** ：`E'_1` 严格守恒（精度 = RK2 积分误差）
- **多通道同时开启** ：无单一守恒量；每个 `E'_k` 都会漂；漂幅反映该通道对粒子的耦合强度
- 注意 ：若粒子轨道远离某通道的径向窗口（如 ψ_peak ± few·σ 之外），该通道贡献近零，相应 `E'_k` 漂幅会被 `E`、`P_zeta` 同源涨落主导（Test B 现象）

诊断建议步骤：

1. 单开 ch1，看 `B.Eprime(:,:,1)` 漂幅应 ≲ 1e-4·E（RK2 噪声）
2. 单开 ch2，确认 `B.Eprime(:,:,2)` 类似
3. 三通道同开，比较各 `Eprime(:,:,k)` 漂幅
4. 若想强制单通道生效：把其余 channels 的 `phi_real`/`apara_real` 全置 0

---

## 8. 常见陷阱

| 症状 | 原因 | 修法 |
|---|---|---|
| `ant.channels` 字段不存在 | `ant_interface` 直接用，没跑 `ant_attach_channel` | 加一行 `ant = ant_attach_channel(ant)` |
| `eval_perturb_field` 输出全 0 | `ant.channels(k).num_modes==0` 或 `m_modes` 越出 `phim_x` 列 | 检查 `ant_interface(file, amp, mmin, mmax)` 中 mmin/mmax 是否在 `f.m1`..`f.m2` 内 |
| 轨道指数发散到 1e+200+ | `amp_phi` 过大 → 合成包络 `∂φ/∂ψ` 太陡 | 振幅降一两个量级；或换真实本征模 profile |
| `E'_k` 三通道漂幅几乎一样 | 合成包络 + 大振幅下 `E(t)`、`P_zeta(t)` 自身波致涨落主导 | 用真实窄径向窗口本征模；或单独跑每通道作 baseline 对比 |
| `unit.gtc_utime` 不存在 | `physics_unit.m` 未跑 | 先 `run physics_unit;`（依赖 `plasma.b0` 已设） |

---

## 9. 何时升级出货码

`poincare_perturb.m` 与 `cal_char_freq_perturb.m` 仍走原内联 m 循环。当：

- 主要工作流已迁到 multi-channel
- Test A/B 在真实多 n 数据上回归验证通过
- 单通道 paper 复现已用新路径 byte-equivalent 对比

可考虑做一次 swap：把 `poincare_perturb.m` 的 `for i = 1:ant.num_modes` 块替换为 `eval_perturb_field(ant.channels, ...)` 调用。在此之前**不要碰**。

---

## 10. 与未来谱方法的关系

见 `wiki/concepts/mosaic-multichannel-antenna.md` 与本 repo 的 `docs/status_report.pdf` §6。摘要：

- **Path 1 (CoM Cheb 采样)** — 与多频正交，本脚手架不需改动即兼容。
- **Path 3 (Floquet + SEM)** — 单频严格；多频走 quasi-periodic Floquet，求解器接口 `floquet_matrix(channels)` 写好后扩展。
- **Path 2 (Time-harmonic Galerkin)** — 多频严重冲突，留接口不落地。

---

## 11. API 速查

### `ch = ant_make_channel(omega, n_tor, m_modes, psi_grid, phi_real, phi_imag, apara_real, apara_imag, dptdp_real, dptdp_imag, dapdp_real, dapdp_imag)`

打包单通道。`phi_real`、`phi_imag` 等是 cell array，每元胞对应一个 m。

### `ant = ant_attach_channel(ant, n_tor)`

把 legacy `ant.*` 包成 `ant.channels(1)`。`n_tor` 可选（如 `ant.n` 已设则可省）。

### `ant_multi = ant_merge_channels(ant1, ant2, ...)`

每输入必须已含 `.channels`。返回总 `ant`，`.channels(1..N)` + 顶层镜像 channel 1。

### `F = eval_perturb_field(channels, pdum, tdum, zdum, t_now, amp_mod)`

返回：
- `F.dptdp, F.dptdt, F.dptdz` — 总 δφ 及导数
- `F.apara, F.dapdp, F.dapdt, F.dapdz, F.paparapt` — δA‖ 及导数
- `F.per_channel(k).<same fields>` — 逐通道分量

### `B = push_orbit_multichannel(zpart0, mstep, tstep, qpart, apart, eq_curr, eq_gradB, ant_multi, amp_mod)`

精简 RK2 多通道扰动推进器。返回结构：
- `B.t, B.psi, B.theta, B.zeta, B.rho_par` — 轨道
- `B.E, B.Pzeta` — 守恒量分量
- `B.Eprime(mstep, np, nc)` — 逐通道 E'_k(t)
- `B.confine{m}` — `'inside'` / `'lost'`

---

## 12. 验证清单 (建议落地顺序)

- [x] Test A 在解析平衡上跑通 → 三 n 共振线叠加图
- [x] Test B 在合成 Gaussian 通道上跑通 → 多通道 E'_k 诊断框架可用
- [ ] **下一步** ：在真实 RSAE/TAE eigenfunction (`.mat` 文件) 上跑 Test B，看通道间漂幅是否分化
- [ ] 单通道 byte-equivalent 回归 ：把 `poincare_perturb.m` 输出与 `push_orbit_multichannel(ant)` (N=1) 对比，应 < 1e-12 相对误差
- [ ] 多 n 真实数据上完整复现某历史 Poincare 图
