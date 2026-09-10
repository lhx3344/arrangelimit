---
title: CVX
published: 2026-08-27
description: '这是一篇关于Matlab工具箱CVX的学习笔记，用于后续的凸优化处理'
image: ''
tags: [凸优化，CVX]
category: '算法优化'
draft: false
lang: ''
slug: markdown-CVX(1)
---

## CVX基础语法
* **框架**：cvx_begin + 声明变量 + 目标 + 约束 + cvx_end

| 命令 | 作用 |
|:----|:-----|
| cvx_begin/cvx_end | 必须成对出现，包裹所有优化代码 |
| quiet | 不显示求解过程，优化速度（加在cvx_begin后面）|
| sdp | 启用半正定规划模式（处理矩阵半正定约束时使用）|

* **变量**：用variable声明，用法如下表 

| 代码 | 含义 |
|:-----|:-----|
| variable x | 声明实数标量 |
| variable x(n) | 声明n维实数向量 |
| variable X(m,n) | 声明m x n维实数矩阵 |
| variable x(n) complex | 声明n维复数向量 |
| variable X(m,n) complex | 声明mxn维复数矩阵 |
| variable x(n) nonnegative | 声明非负向量( $x \geqslant 0$ ) |
| variable X(n,n) symmetric | 声明实对称矩阵 |
| variable X(n,n) hermitian | 声明复 Hermitian 矩阵 (即 X = X^H)|
| variable X(n,n) semidefinite | 声明半正定矩阵 ($X \succeq 0$) |

* **核心约束规则**：整体约束必须为凸约束

## CVX关键规则
CVX不是任何数学表达式都能就受的。他要求整个问题必须是凸的。

### **规则一**：目标函数必须符合凸性
* minimize 后面的表达式必须是凸函数（如norm、sum、quad_form）
* maximize 后面的表达式必须是凹函数（如log、geo_mean）

### **规则二**：约束的左右结构有限制
| 约束类型 | 左侧要求 | 右侧要求 | 合法示例 | 非法示例 |
|:-------:|:-------:|:--------:|:-------:|:-------:|
| <= | 凸 | 凹 | norm(x) <= 1 | 1 <= norm(x) |
| >= | 凹 | 凸 | 1 >= norm(x) | norm(x) >= 1 |
| == | 仿射 | 仿射 | A*b == b | norm(x) == 1 |

**常用函数的凹凸性:**
| 类型 | 常见函数 | 说明 |
|:-----|:-------|:-----|
| 仿射 | +, -, *, /（常数除法）, real, imag, conj, trace, diag, sum, vec, reshape | 线性运算，既是凸也是凹 |
| 凸 | norm, norm(x, p), quad_form, sum_square, max，abs（实数），square, inv_pos, log_sum_exp | 可用于 <= 的左侧或minimize |
| 凹 | log,log_det,geo_mean,sqrt | 可用于 >= 的左侧或maximize |

## 相关代码示例
```matlab title="my_first_cvx_code"
clear; clc;
K = 4; % 用户数
M = 4; % 基站天线数
sigma = 1; % 噪声标准差（噪声功率 = 1）
num_realizations = 10; % 每个gamma值下的随机信道次数（可调，大一点更平滑）

% 定义要扫描的gamma范围（线性值）
gamma_vec = 0.5:0.5:4; % 从 0.5 到 4，步长 0.5
avg_power = zeros(size(gamma_vec)); % 存储每个 gamma 对应的平均功率

fprintf('开始扫描不同 gamma 值……\n');
for g_idx = 1:length(gamma_vec)
    gamma = gamma_vec(g_idx);
    power_sum = 0;
    valid_count = 0;

    for z = 1:num_realizations
        H = randn(K,M); % 随机信道

        cvx_begin quiet
            variable tau nonnegative
            variable W(M,K) complex
            minimize(tau)
            subject to
                for i = 1:K
                    % SINR 约束（含噪声）
                    norm([H(i,:) * W,sigma]) <= sqrt(1 + 1/gamma) * real(H(i,:) * W(:,i));
                    % 强制接收信号为实数且非负
                    imag(H(i,:) * W(:,i)) == 0;
                    real(H(i,:) * W(:,i)) >= 0;
                end
                norm(vec(W)) <= sqrt(tau);
        cvx_end

        % 检查求解状态，只累加成功求解的结果
        if strcmp(cvx_status, "Solved")
            power_sum = power_sum + tau;
            valid_count = valid_count + 1;
        else
            % 若不行，忽略该次实现（或可记录失败次数）
            % fprintf('gamma = %.2f,第 %d 次不可行\n',gamma,z);
        end
    end

    if valid_count > 0
        avg_power(g_idx) = power_sum / valid_count;
        fprintf('gamma = %.2f，有效次数 = %d，平均功率 = %.4f\n',gamma,valid_count,avg_power(g_idx));
    else
        avg_power(g_idx) = NaN; % 若全部不可行，置NaN
        fprintf('gamma = %.2f，无可行的信道实现!\n',gamma);
    end
end

% 绘图
figure;
plot(gamma_vec,avg_power,'b-o','LineWidth',2);
xlabel('SINR 门限 \gamma（线性值）');
ylabel('平均最小发射功率(||W||_F^2)');
title(['M=', num2str(M), ', K=', num2str(K), ', 噪声功率 \sigma^2=', num2str(sigma^2), ...
       ', 每个\gamma下随机信道次数=', num2str(num_realizations)]);
grid on;
```

---

## 信号模型（只包含公式）

$$
x_i = w_i s_i 
$$

$$
y_i = h_i w_i s_i + \sum_{j\neq i}h_i w_j s_j + n_i
$$

## 功率最小化问题

基站总发射功率为 $\sum_{i=1}^{K}Tr(X_u)$ 。用户i的信干噪比（SINR）定义为

$$
SINR_i = \frac{h_i X_i h_i^H}{\sum_{j \neq i}h_i X_j h_i^H + \sigma^2} \geq \gamma
$$

要求 $SINR_i \geq \gamma$，其中 $\gamma \geq 0$ 是目标门限。则优化问题可写为

$$
\min_{\{Xi\}} \sum_{i=1}^{K}Tr(X_i)
$$

$$
\begin{aligned}
\text{subject to} \quad \frac{h_i X_i h_i^H}{\sum_{j \neq i}h_i X_j h_i^H + \sigma^2} \geq \gamma, \forall i, \\
X_i \succeq 0, rank(X_i) \leq 1, \forall i,
\end{aligned}
$$

其中 $X_i \succeq 0$ 表示半正定。秩一约束来自 $ X_i = w_i w_i^H $ ,这使得问题非凸。

## 半正定松弛
### 松弛问题
松弛掉秩一约束，仅保留半正定约束，得到松弛问题

$$
\underset{\{X_i\} , \{s_i\}}{\text{minimize}} \quad \sum_{i = 1}^{K} Tr{X_i}
$$

$$
\begin{aligned}
\text{subject to} \quad & \operatorname{Tr}(H_i X_i) - \gamma \sum_{j \neq i} \operatorname{Tr}(H_i X_j) \geq \gamma \sigma^2, \quad \forall i, \\
& X_i \succeq 0, \quad \forall i
\end{aligned}
$$

其中 $ H_i = h_i^H h_i $ 是 M x M 的秩一Hermitian矩阵，且 $ H_i \succeq 0 $ 。为将不等式约束转化为为等式以便数值求解，引入非负松弛变量 $ s_i \geq 0 $ ，使

$$
\operatorname{Tr}{H_i X_i} - \gamma \sum_{j \neq i} \operatorname{Tr}{H_i X_j} - s_i = \gamma \sigma^2
$$

这样约束变为线性等式，目标函数线性，可行域为半正定锥于线性等式的交集，因此是凸的半定规划问题。  
SDP松弛的最优值给出了原问题的最优功率的下界。若存在最优解满足秩一，则松弛是紧的，该解即为原始问题的最优波形成束。即使出现秩大于一1，也可通过高斯随机化等后处理步骤获得良好的可行波束成形向量。
> SDP的最优解 $ P_{SDP} $ 是原始问题最优解 $P^*$ 的下界，即 $P_{SDP} \leq P^*$ 。然而，SDP的解 $ X_i^* $ 不一定满足秩一条件，无法直接作为波束成形向量使用。高斯随机化正是通过概率抽样从 $ X_i^* $ 中生成候选向量，并从中筛选出满足原始约束的可行解，从而获得 $P^*$ 的一个上界

## 相关代码实现
```matlab title="cvx_SDP.m"
clear;clc;

%======= 系统参数 =======
K = 4;              % 用户数
M = 4;              % 基站天线数
sigma2 = 0.025^2;   % 噪声功率
num_realizations = 20;  % 每个 gamma 下的信道样本数

% 扫描的SINR门限（线性值）
gamma_vec = 0.5:0.5:3;
avg_power = zeros(size(gamma_vec));

fprintf('开始扫描不同 gamma 值（SDP方法）……\n');

for g_idx = 1:length(gamma_vec)
    gamma = gamma_vec(g_idx);
    power_sum = 0;
    valid_count = 0;

    for z = 1:num_realizations
        % ----- 生成信道 ------
        H_ch = randn(K,M);      % 实高斯信道（方差一）
        H = zeros(M,M,K);       % M * M * K维矩阵
        for i = 1:K
            h_i = H_ch(i,:);    % 1 * M
            H(:,:,i) = h_i' * h_i; % M * M 半正定
        end

        % ------ SDP 求解 --------
        cvx_begin quiet
            variable X(M,M,K) complex
            variable s(K,1) nonnegative

            % 目标函数：总发射功率（取实部）
            obj = 0;
            for i = 1:K
                obj = obj + real(trace(X(:,:,i)));
            end
            minimize(obj)

            subject to
                for i = 1:K
                    % 信号项
                    sig = real(trace(H(:,:,i) * X(:,:,i)));
                    % 干扰项
                    cstr = 0;
                    for j = 1:K
                        if j ~= i
                            cstr = cstr + real(trace(H(:,:,i) * X(:,:,j)));
                        end
                    end
                    % SINR 约束
                    sig - gamma * cstr - s(i) == gamma * sigma2;
                    % 半正定约束
                    X(:,:,i) == hermitian_semidefinite(M);
                end
        cvx_end

        if strcmp(cvx_status, 'Solved')
            power_sum = power_sum + obj;
            valid_count = valid_count + 1;
        end
    end

    if valid_count > 0
        avg_power(g_idx) = power_sum / valid_count; % 取结果的平均值一克服随机性（蒙特卡洛仿真）
        fprintf('gamma = %.2f，有效次数 = %d，平均功率 = %.4f\n',gamma,valid_count,avg_power(g_idx));
    else
        avg_power(g_idx) = NaN;
        fprintf('gamma = %.2f，无可实行的信道实现!\n',gamma);
    end
end

% ========== 绘图 ==========
figure;
plot(gamma_vec, avg_power, 'b-o', 'LineWidth', 2);
xlabel('SINR 门限 \gamma (线性值)');
ylabel('平均最小发射功率 (总功率)');
title(['M=', num2str(M), ', K=', num2str(K), ', 噪声功率 \sigma^2=', num2str(sigma2), ...
       ', 每个\gamma下随机信道次数=', num2str(num_realizations)]);
grid on;
```
---
## 最优性分析
* **SDR求解的是松弛问题**：将原始非凸问题（秩一约束）放宽为半正定约束，可行域扩大，因此松弛问题的最优值 $ P_{SDR} $必定是原始问题最优值 $P^*$ 的下界
* **随机化构造的是可行解**：从SDR解出的协方差矩阵 $ X_i $中采样，得到一组 $ w_i $。只有当这些向量满足原始约束（如SINR门限）时，它们才是原问题的可行解。因此，该可行解对应的目标值（总功率）必定大于等于原始问题的最优值，即构成一个上界：

$$
P_{rand} \geq P^*
$$

将两者结合，必然得到：

$$
P_{SDR} \leq P^* \leq P^{rand}
$$

## 高斯随机化流程
在每一次信道实现中，完成SDP求解得到最优协方差矩阵 $\{X_i^*\}$ 后，执行以下高斯随机化流程：
1. 随机向量生成

    对每个用户i，从 $X_i^*$ 中独立采样一个复高斯随机向量 $r_i$。具体的，先对 $ X_i^* $ 进行特征分解并舍去负特征值（保证数值半正定），计算其平方根，再乘以标准复高斯向量，得到

    $$
    r_i = X_i^{* \\ 1/2} * \frac{u_i + j v_i}{\sqrt{2}},
    $$

    其中 $ u_i,v_i \thicksim N(0,I_M) $ 独立同分布。这样保证了 $E[r_i r_i^H] = X_i^*$
2. 可行性检验

    将生成的候选波束成形矩阵 $R = [r_1,...,r_K]$ 代入原始SINR表达式，计算每个用户的实际SINR：

    $$
    SINR_i = \frac{r_i^H H_i r_i}{\sum_{j \neq i}r_j^H H_i r_j + \sigma^2}
    $$
    
    若所有用户均满足 $ SINR_i \geq \gamma $，则该候选向量构成原问题的一个可行解；否则丢弃。
3. 最优可行解选取

    对通过检验的候选解，计算其总发射功率 $P_rand = \sum_{i} \begin{Vmatrix} r_i \end{Vmatrix}^2$。重复采样过程足够多次（通常需保证采样充分覆盖可行域），记录所有可行解中功率最小的一个，作为该信道实现的随机化上界。
4. 统计平均

    对同一 $\gamma$ 下的所有信道实现，若SDP可行且至少有一个随机化可行解，则分别累加SDP目标值和随机化最优功率，最后求平均，得到该 $\gamma$ 对应的平均下界和平均上界曲线。

## 相关代码实现
```matlab title="cvx_SDR.m"
clear;clc;close all;

% ============= 系统参数 ===============
K = 4;                      % 用户数
M = 4;                      % 基站天线数
sigma2 = 0.025^2;           % 噪声功率
num_realizations = 20;      % 每个 gamma 下的信道样本数
L = 200;                    % 随机化采样次数

gamma_vec = 0.5:0.5:3;
avg_power_sdp = zeros(size(gamma_vec));
avg_power_rand = zeros(size(gamma_vec));

fprintf('开始扫描不同gamma值（SDP + 随机化）...\n');

for g_idx = 1:length(gamma_vec)
    gamma = gamma_vec(g_idx);
    power_sum_sdp = 0;
    power_sum_rand = 0;
    valid_count = 0;

    for z = 1:num_realizations
        % ----- 生成信道 -------
        H_ch = randn(K,M);
        H = zeros(M,M,K);
        for i = 1:K
            h_i = H_ch(i,:);
            H(:,:,i) = h_i'*h_i;
        end

        % ----- SDP松弛求解（下界）------
        cvx_begin quiet
            variable X(M,M,K) complex
            variable s(K,1) nonnegative
            obj = 0;
            for i = 1:K
                obj = obj + real(trace(X(:,:,i)));
            end
            minimize(obj)
            subject to
                for i = 1:K
                    sig = real(trace(H(:,:,i) * X(:,:,i)));
                    cstr = 0;
                    for j = 1:K
                        if j ~= i
                            cstr = cstr + real(trace(H(:,:,i) * X(:,:,j)));
                        end
                    end
                    sig - gamma * cstr - s(i) == gamma * sigma2;
                    X(:,:,i) == hermitian_semidefinite(M);
                end
        cvx_end

        if ~strcmp(cvx_status,'Solved')
            continue;
        end

        % ---- 高斯随机化（构造上界）-------
        best_rand_power = inf;
        
        for l = 1:L
            W_rand = zeros(M,K);
            for i = 1:K
                Xi = X(:,:,i);
                [V,D] = eig(0.5*(Xi + Xi')); % 求矩阵A的全部特征值，构成对角阵D，并求A的特征向量构成V的列向量。
                D = max(real(D),0);
                Xi_sqrt = V * sqrt(D) * V';
                r = Xi_sqrt * (randn(M,1) + 1i*randn(M,1)) / sqrt(2);
                W_rand(:,i) = r;
            end

            % 计算每个用户的SINR
            SINRs = zeros(K,1);
            for i = 1:K
                w_i = W_rand(:,i);
                % 信号功率
                sig = real(w_i' * H(:,:,i) * w_i);
                % 干扰功率（来自其它用户）
                int = 0;
                for j = 1:K
                    if j ~= i
                        w_j = W_rand(:,j);
                        int = int + real(w_j' * H(:,:,i) * w_j);
                    end
                end
                SINRs(i) = sig / (int + sigma2);
            end

            % 检查是否所有用户满足SINR >= gamma
            if all(SINRs >= gamma - 1e-6)
                total_power = sum(sum(abs(W_rand).^2,1));
                if total_power < best_rand_power
                    best_rand_power = total_power;
                end
            end
        end

        if isinf(best_rand_power)
            continue;
        end;

        power_sum_sdp = power_sum_sdp + obj;
        power_sum_rand = power_sum_rand + best_rand_power; % 20个随机信道功率的累加
        valid_count = valid_count + 1;
    end

    if valid_count > 0
        avg_power_sdp(g_idx) = power_sum_sdp / valid_count;
        avg_power_rand(g_idx) = power_sum_rand / valid_count;
        fprintf('γ = %.2f: 有效样本=%d, SDP下界=%.4f, 随机化上界=%.4f\n', ...
                gamma, valid_count, avg_power_sdp(g_idx), avg_power_rand(g_idx));
    else
        avg_power_sdp(g_idx) = NaN;
        avg_power_rand(g_idx) = NaN;
        fprintf('y = %.2f：无可行样本！\n',gamma);
    end
end

% ==================== 绘图 ====================
figure;
plot(gamma_vec, avg_power_sdp, 'b-o', 'LineWidth', 2, 'DisplayName', 'SDP 下界 (松弛)');
hold on;
plot(gamma_vec, avg_power_rand, 'r-s', 'LineWidth', 2, 'DisplayName', '随机化上界 (可行解)');
hold off;
xlabel('SINR 门限 \gamma (线性值)');
ylabel('平均最小发射功率');
title(['M=', num2str(M), ', K=', num2str(K),...
       ', 每个γ下信道数=', num2str(num_realizations), ', 随机化次数=', num2str(L)]);
legend('Location', 'best');
grid on;
```
---
## 情景应用（MU-MISO下行功率分配问题）
在多天线通信系统中，物理层多播（multicast）是一类重要的传输模式，其中基站向多个用户同时发送相同的公共数据流，如视频会议、软件更新或系统信息广播。与多用户独立数据流不同，多播传输不需要区分用户数据，因此不存在多用户干扰，设计目标转变为如何设计单个波束成形向量，使得所有用户都能获得足够高的接收信号质量。  
最大最小公平性（max-min fairness）准则在此场景下自然适用：在总发射功率受限的条件下，最大化所有用户中最小的接收功率，从而保证所有用户的基本服务质量。该问题在数学上可表述为以波束成形向量为变量的非凸二次优化问题，其难点在于秩一约束的存在。  
半定松弛（SDR）是求解此类非凸二次问题的经典方法，通过将秩一约束松弛为半正定约束，将问题转化为凸半定规划（SDP），从而可高效求得全局最优上界。然而，松弛解通常不满足秩一条件，无法直接实现。为此，需要采用后处理技术，如秩一近似或高斯随机化，从松弛解中构造可行的波束成形向量。  
本文针对多播场景下的最大最小公平性问题，基于半定松弛框架，结合秩一近似和高斯随机化两种构造方法，通过蒙特卡洛仿真评估它们的性能差异，并验证SDR松弛的紧性。

## 信号模型
考虑但小区多输入单输出（MISO）下行链路，基站配备M根天线，同时服务K个单天线用户。在物理层多播场景下，基站发送一个公共数据流s（$|s|^2=1$），采用单个波束成形向量 $w \in \cnums^{M \times 1}$ 进行预编码。发射信号为

$$
x = ws,
$$

总发射功率为

$$
P_tx = \lVert w \rVert^2 = tr(w w^H)
$$

设用户i的信道向量为 $ h_i \in \cnums^{M \times 1} $，建模为独立同分布的瑞利衰落信道，即 $h_i \sim \mathcal{CN}(\mathbf{0},\mathbf{I}_M)$。用户i的接收信号为 

$$
y_i = h_i^H w s + n_i
$$

其中 $ n_i \sim \mathcal{CN}(\mathbf{0},\sigma^2)$ 为加性高斯白噪声。    
用户i的接收信号功率（忽略噪声）为

$$
P_i(w) = |h_i^H w|^2 - w^H H_i w,
$$

其中 $H_i = h_i h_i^H $ 是秩一Hermitian半正定矩阵

## 最大最小公平性优化问题
在总功率约束 $ P_{tx} \geq P_{total} $ 下，最大化所有用户中的最小接收功率，即
$$
\max_{w} \quad \min_{i=1,...,K} |h_i^H w|^2
$$

$$
\text{s.t.} \quad \lVert w \rVert^2 \geq P_{total}
$$

引入辅助变量 $ t \leq 0 $，该问题等价于
$$
\max_{w,t} \quad t
$$

$$
\tag{1} \begin{aligned}
\text{s.t.} \quad w^H H_i w \geq t, \quad i = 1,...,K, \\
\lVert w \rVert^2 \leq P_{total}
\end{aligned}
$$

问题（1）是非凸的，因为目标函数和约束均为凸二次型，但最大化凸函数（或最小化凹函数）通常是非凸的。然而，通过定义 $W=w w^H$，可将其转换为秩一约束的半定问题

令 $X=w w^H$，则 $tr(X) = \lVert w \rVert$，且 $w^H H_i w=tr(H_i X)$。于是问题(1)可写为

$$
\begin{equation}
\begin{aligned}
\max_{X, t} \quad & t \\
\text{s.t.} \quad & \operatorname{tr}(H_i X) \geq t, \quad i = 1, \dots, K, \\
& \operatorname{tr}(X) \leq P_{\text{total}}, \\
& X \succeq 0, \quad \operatorname{rank}(X) = 1
\end{aligned}
\tag{2}
\end{equation}
$$

问题(2)的难点在于秩一约束rank(X) = 1，该约束非凸。

## 半定松弛（SDR）
将问题(2)中的秩一约束松弛为半正定约束，仅保留 $X \geq 0$,得到凸半定规划：

$$
\begin{equation}
\begin{aligned}
\max_{X, t} \quad & t \\
\text{s.t.} \quad & \operatorname{tr}(H_i X) \geq t, \quad i = 1, \dots, K, \\
& \operatorname{tr}(X) \leq P_{\text{total}}, \\
& X \succeq 0,
\end{aligned}
\tag{3}
\end{equation}
$$

问题（3）是凸的，可用内点法（如CVX）高效求解。记其最优解为 $(X^*,t_{\text{SDP}})$。可由于松弛问题的可行域包含了原始问题的可行域（秩一矩阵集合），故

$$
t_{\text{SDP}} \geq t^*,
$$

其中t*是原始问题的最优值。因此，SDP给出的是真实最优解的上界。

## 秩一近似
从松弛解 $X^*$ 出发，一种最简单的构造可行解的方法是利用其主特征向量。对 $X^*$ 进行特征分解，取最大特征值 $\lambda_{\text{max}}$ 及其对应的单位特征向量 $v_{\text{max}}$，构造波束成形向量：

$$
w_{\text{rank1}} = \sqrt{P_{\text{total}}}v_{\text{max}},
$$

以保证总共率恰好为 $P_{\text{total}}$（注意 $v_{\text{max}}$ 为单位范数）。该向量满足功率约束，且对应的最小接收功率为

$$
t_{\text{rank1}} = \min_{t}w_{\text{rank1}}^H H_i w_{\text{rank1}}.
$$

由于 $w_{\text{rank1}}$ 是原始问题的可行解，必然有

$$
t_{\text{rank1}} \leq t^*
$$

因此，秩一近似给出一个下界

## 高斯随机化
为了获得更紧的下界，采用高斯随机优化方法，其核心思想是从X*的协方差结构中随机采样多个候选向量，并从中选取满足约束性能最好的一个，具体步骤如下：
1. 计算 $X^*$ 的半正定平方根 $X^{* \\ 1/2}$，即 $X^*=X^{* \\ 1/2}(X^{* \\ 1/2})^H$(通过特征分解，并舍去负特征值以保证数值稳定性)
2. 对于 l = 1,2,...,L:
    * 生成复高斯随机向量 $z_l \sim \mathcal{CN} $
    * 计算候选向量 $\tilde{w_l} = X^{* \\ 1/2} z_i$
    * 归一化以满足功率约束：$w_i = \sqrt{P_{\text{total}}} \frac{\tilde{w_l}}{\lVert \tilde{w_l} \rVert}$
    * 计算 $t_l = \min_i w_l^H H_i w_l$
3. 取所有候选中的最大值：$t_{\text{rand}} = \max_l t_l$

由于每个 $w_l$ 都满足原始约束，因此 $t_{\text{rand}} \leq t^*$，即也是一个下界，随着采样次数L增大，$t_{\text{rand}}$ 可逼近 t*,通常优于秩一近似。  
综上，对于任意一信道实现，我们有如下不等式链：

$$
t_{\text{rand}} \leq t^* \leq t_{\text{SDP}}, \quad t_{\text{rank1}} \leq t^*
$$

其中秩一近似和随机化均给出下界，但随机化通常更优

## 相关代码
```matlab title="cvx_MISO.m"
clear;clc;close all;

% 系统参数
K = 4;                      % 用户数
M = 4;                      % 基站天线数
num_realizations = 50;      % 每个功率值下的信道样本数（可调）
L = 200;                    % 高斯随机化采样次数
P_total_vec = 0.5:0.5:3.5;  % 扫描的总发射功率范围

% 预分配存储
avg_sdp = zeros(size(P_total_vec));
avg_rank1 = zeros(size(P_total_vec));
avg_random = zeros(size(P_total_vec));

fprintf('开始扫描不同发射功率下的最大最小公平性...\n');

for p_idx = 1:length(P_total_vec)
    P_total = P_total_vec(p_idx);
    sum_sdp = 0; sum_rank1 = 0; sum_random = 0;
    valid_count = 0;

    for z = 1:num_realizations
        % 生成信道（K个用户，每个用户信道向量 MX1）
        Hi_s = randn(M,M);
        H = zeros(M,M,K);
        for i = 1:K
            h_i = Hi_s(:,i);
            H(:,:,i) = h_i * h_i';
        end
        
        % 1.SDP松弛（上界）
        cvx_begin quiet
            variable X(M,M) complex
            variable t
            maximize(t)
            subject to
                for i = 1:K
                    real(trace(X * H(:,:,i))) >= t
                end
                trace(X) == P_total
                X == hermitian_semidefinite(M)
        cvx_end
        if ~strcmp(cvx_status,'Solved')
            continue;  % 若不可行则跳过该信道
        end
        sdp_val = t;

        % 2.秩一近似（可行化，归一化功率）
        [V,D] = eigs(X,1);
        w_rank1 = V;        % 单位范数
        % 缩放以满足总功率 P_total（若需要，这里保持功率为 P_total）
        % 因为 V 是单位向量，所以功率为1，而我们需要总功率为P_total，但SDP解X的总功率为P_total
        % 其最大特征值对应的向量范数平方为最大特征值lambda，可能小于 P_total，我们直接采用 w =sqrt(P_total)*V使功率为 P_total
        % 但这样做可能不满足原始约束？更合理的做法是从x中提取主成分并缩放
        % 这里我们直接采用w = sqrt(P_total) * V，即总功率恰好为 P_total
        w_rank1 = sqrt(P_total) * V;
        t_rank1 = inf;
        for i = 1:K
            p_i = real(w_rank1' * H(:,:,i) * w_rank1);
            t_rank1 = min(t_rank1,p_i);
        end

        % 3.高斯随机化
        [Vx,Dx] = eig(0.5*(X + X'));
        Dx = max(real(Dx),0);
        X_sqrt = Vx * sqrt(Dx) * Vx';
        best_t = -inf;
        for l = 1:L
            z = (randn(M,1) + 1i*randn(M,1)) / sqrt(2);
            w = X_sqrt * z;
            % 归一化到功率 P_total
            w = w / norm(w) * sqrt(P_total);
            t_min = inf;
            for i = 1:K
                p_i = real(w' * H(:,:,i) * w);
                t_min = min(t_min,p_i);
            end
            if t_min > best_t
                best_t = t_min;
            end
        end

        % 累加
        sum_sdp = sum_sdp + sdp_val;
        sum_rank1 = sum_rank1 + t_rank1;
        sum_random = sum_random + best_t;
        valid_count = valid_count + 1;
    end

    if valid_count > 0
        avg_sdp(p_idx) = sum_sdp / valid_count;
        avg_rank1(p_idx) = sum_rank1 / valid_count;
        avg_random(p_idx) = sum_random / valid_count;
         fprintf('P_total = %.2f: 有效信道数=%d, SDP=%.4f, 秩一=%.4f, 随机=%.4f\n', ...
                P_total, valid_count, avg_sdp(p_idx), avg_rank1(p_idx), avg_random(p_idx));
    else
        avg_sdp(p_idx) = NaN;
        avg_rank1(p_idx) = NaN;
        avg_random(p_idx) = NaN;
        fprintf('P_total = %.2f: 无可行信道！\n', P_total);
    end
end

% ==================== 绘图 ====================
figure;
plot(P_total_vec, avg_sdp, 'b-o', 'LineWidth', 2, 'DisplayName', 'SDP 上界');
hold on;
plot(P_total_vec, avg_rank1, 'r--s', 'LineWidth', 2, 'DisplayName', '秩一近似');
plot(P_total_vec, avg_random, 'g-.d', 'LineWidth', 2, 'DisplayName', '高斯随机化');
hold off;
xlabel('总发射功率 P_{total}');
ylabel('平均最小接收功率 (t)');
title(sprintf('最大最小公平性：M=%d, K=%d, 随机化次数 L=%d', M, K, L));
legend('Location', 'northwest');
grid on;
```

---
## 基于连续凸近似（SCA）的求解方法（应用情景忽略）
### 引入辅助变量与干扰项
> 此处噪声项归一化为 1 
为处理非凸的SINR箱，定义辅助变量 $\rho_{p,k}$ 和 $\rho_{c,k}$，分别表示私有和公共流的SINR近似值，以及辅助变量 $\alpha_{p,k}$ 和 $\alpha_{c}$ 表示相应的速率。同时，引入干扰加噪声项 $\beta_{p,k}$ 和 $\beta_{c,k}$ :

* 私有流解码时（用户k）的干扰加噪声：

$$
\beta_{p,k} = \sum_{j \neq k}|h_k^H p_j|^2 + 1 
$$

* 公共流解码时（用户k）的干扰加噪声：

$$
\beta_{c,k} = \sum_{j = 1}^{K}|h_k^H p_j|^2 + 1
$$

于是有：

$$
R_{p,k} = log_2(1 + \frac{|h_k^H p_k|^2}{\beta_{p,k}}), \quad R_{c,k} = log_2(1 + \frac{|h_k^H p_c|^2}{\beta_{c,k}}).
$$

### 非凸项的线性化
函数 $f(p,\beta) = \frac{|h^H p|^2}{\beta}$ 在联合变量 $(p,\beta)$ 上非凸。为了构建凸近似，我们在当前迭代点 $(p^{(n-1)},\beta^{(n-1)})$ 处对其进行一阶泰勒展开（对两个变量同时线性化），得到其线性近似：

$$
f(p,\beta) \approx \frac{2\Re{\{(h^H p^{(n-1)})^* (h^H p^{(n-1)})\}}}{\beta^{(n-1)}} - \frac{|h^H p^{(n-1)}|^2}{(\beta^{(n-1)})^2}\beta
$$

改近似关于p和 $\beta$ 均为线性，因此是凸函数。将此线性化带入SINR约束，可将非凸的SINR下线约束转换为线性约束，从而保证子问题的凸性

### 优化变量与辅助变量定义
在第n此迭代中，已知前一迭代点 $\{P^{(n-1)},p_c^{(n-1)},\beta_{p}^{(n-1)},\beta_{c}^{(n-1)}\}$，其中 $\beta_{p} = [\beta_{p,1},...,\beta_{p,K}]$，$\beta_c = [\beta_{c,1},...,\beta_{c,K}]^T$  
未将原非凸问题转化为凸形式，引入以下辅助变量：
* $t \in \reals$：表示所有用户最小总速率的下界 
* $ P = [p_1,...,p_K] \in \cnums^{N_t X K}$：私有波束成形矩阵；
* $p_c \in \cnums^{N}：公有波束成形向量
* $\rho_{p,k} \in \reals$：用户k私有流SINR的近似变量；
* $\rho_{c,k} \in \reals$：用户k公共流SINR的近似变量
* $\alpha_{p,k} \in \reals$：用户k私有速率的近似变量
* $\alpha_{c} \in \reals$：公共总速率的近似变量
* $c_k \in \reals$：公共速率分配给用户k的份额

其中k = 1，2，……，K。  
则子凸问题可以转换为：  
**目标函数：**

$$
max \quad t
$$

**约束条件：**
1. 最小用户总速率保证（对所有k=1,……,K）:
    $$
    t \leq \alpha_{p,k} + c_k
    $$

2. 公共速率分配总额约束：
    $$
    \sum_{k=1}^(K) c_k \geq \alpha_c
    $$

3. 私有速率与SINR的凹关系（对所有k）：
    $$
    \alpha_{p,k} \leq log_2(1 + \rho_{p,k})
    $$

4. 公共速率与SINR的凹关系（对所有k）：
    $$
    \alpha_{c} \leq log_2(1 + \rho_{c,k})
    $$

5. 私有SINR的线性化上界（对所有k）：
    $$
    \rho_{p,k} \leq \frac{2\Re{\{(h_k^H p_k^{(n-1)})^* (h_k^H p_k)\}}}{\beta_{p,k}^{(n-1)}} - \frac{|h_k^H p_k^{(n-1)}|^2}{(\beta_{p,k}^{(n-1)})^2}\beta_{p,k}
    $$

6. 公共SINR的线性化上界（对所有k）：
    $$
    \rho_{c,k} \leq \frac{2\Re{\{(h_k^H p_c^{(n-1)})^* (h_k^H p_c)\}}}{\beta_{c,k}^{(n-1)}} - \frac{|h_k^H p_c^{(n-1)}|^2}{(\beta_{c,k}^{(n-1)})^2}\beta_{c,k}
    $$

7. 私有干扰项下界（对所有k）：
    $$
    \beta_{p,k} \geq \sum_{j \neq k}|h_k^H p_j|^2 + 1 
    $$

8. 公共干扰项下界（对所有k）：
    $$
    \beta_{c,k} \geq \sum_{j=1}|h_k^H p_j|^2 + 1 
    $$

9. 总发射功率约束（二阶锥约束）：
    $$
    \lVert P \rVert_{F}^2 + \lVert p_c \rVert_2^2 \leq P_t
    $$

10. 公共速率份额非负性（对所有k）：
    $$
    c_k \geq 0
    $$

可以看到以上所有约束聚那位凸约束（线性、二阶锥或对数凹），因此子问题是一个标准的凸优化问题，可被CVX等工具高效求解。

算法流程：
1. 初始化：设置迭代计数器 n=0,选择初始点 $P^{(0)},p_c^{(0)}$,例如
    * $p_c^{(0)} = \sqrt{0.5P_t} * u_1 $，其中 $u_1$ 为 H的最大奇异值对应的左奇异向量；
    * $P^{(0)} = \sqrt{0.5P_t / K} * H/ \lVert H \rVert_F $ （匹配信道方向）
        计算 $\beta_p^{(0)}$ 和 $\beta_c^{(0)}$
2. 迭代（n = 1,2,……）
    * 固定 $\{P^{(n-1)},p_c^{(n-1)},\beta_{p}^{(n-1)},\beta_{c}^{(n-1)}\}$，求解凸子问题，得到 $\{P^{(n)},p_c^{(n)},\beta_{p}^{(n)},\beta_{c}^{(n)}\}$ 及最优目标值 $t^{(n)}$
    * 若 $|t^{(n)} - t^{(n-1)}| \leq \varepsilon$ （预设容差），则停止
    * 否则，更新迭代点并继续
3. 功率归一化：收敛化，令 $P = \frac{\sqrt{P_t}}{\lVert [P;P_c] \rVert_F} P,\quad p_c = \frac{\sqrt{P_t}}{\lVert [P;P_c] \rVert_F}p_c$，以满足总功率约束紧致

### 相关代码

::: code-group labels=[cvx_RSMA_SCA.m, rsma_sca_update.m]

```matlab
%% 系统参数配置 (根据仿真参数表)
clear; clc;
Nt = 16;            % 基站天线数
K = 4;              % 用户数
SNR_dB = 10;        % 发射信噪比 (dB)
Pt = 10^(SNR_dB/10); % 总发射功率 (噪声功率归一化为1)

maxIter = 1000;     % 最大外层迭代次数
tol = 1e-4;         % 收敛容差（目标值变化小于此值则停止）

%% 1. 生成瑞利平坦衰落信道
H = (randn(Nt, K) + 1i*randn(Nt, K)) / sqrt(2); % 元素 ~ CN(0,1)

%% 2. 初始化
% 2.1 公共波束：最大奇异值方向
[U, ~, ~] = svd(H);
u1 = U(:, 1);                       % 最大奇异值对应的左奇异向量
p_c = sqrt(0.5 * Pt) * u1;          % 公共波束初始值 (公共功率占比50%)

% 2.2 私有波束：匹配信道方向 (MRT)，按F范数归一化
P = sqrt(0.5 * Pt / K) * H / norm(H, 'fro');

% 2.3 计算初始 beta_p 和 beta_c (干扰+噪声功率)
beta_p = zeros(K, 1);
beta_c = zeros(K, 1);
for k = 1:K
    % 私有流解码时：剔除自身的私有流，其他私有流视为干扰
    idx = [1:k-1, k+1:K];
    beta_p(k) = sum(abs(H(:, k)' * P(:, idx)).^2) + 1;
end
% 公共流解码时：所有私有流视为干扰
beta_c = 1 + sum(abs(H' * P).^2, 2);

% 存储迭代记录
t_hist = [];
P_last = P;
p_c_last = p_c;
beta_p_last = beta_p;
beta_c_last = beta_c;

%% 3. SCA 外层迭代循环
fprintf('开始 SCA 迭代 (SNR=%d dB, Nt=%d, K=%d)...\n', SNR_dB, Nt, K);
for iter = 1:maxIter
    % 调用 CVX 子问题求解函数
    [t, c, P, p_c, beta_p, beta_c] = ...
        rsma_sca_update(H, Pt, P_last, p_c_last, beta_p_last, beta_c_last);
    
    % 记录目标值
    t_hist = [t_hist, t];
    
    % 判断收敛
    if iter > 1 && abs(t_hist(end) - t_hist(end-1)) < tol
        fprintf('在第 %d 步收敛，目标值 = %.4f nats/s/Hz\n', iter, t);
        break;
    end
    
    % 更新迭代点
    P_last = P;
    p_c_last = p_c;
    beta_p_last = beta_p;
    beta_c_last = beta_c;
end

if iter == maxIter
    fprintf('达到最大迭代次数 %d，当前目标值 = %.4f\n', maxIter, t);
end

%% 4. 结果展示
figure;
plot(t_hist, 'b-o', 'LineWidth', 1.5);
xlabel('迭代次数'); ylabel('最小用户速率 (nats/s/Hz)');
title('RSMA-SCA 收敛曲线');
grid on;

fprintf('公共流功率: %.4f, 私有流总功率: %.4f\n', norm(p_c)^2, norm(P, 'fro')^2);
```

```matlab
% % SCA 迭代中的CVX子问题求解函数
% % 输入： H - 信道矩阵，Pt - 总功率，P_last,p_c_last - 上一次迭代点的波束
% %       beta_p_last,beta_c_last - 上一次迭代点的干扰加噪声项
% % 输出： t - 当前子问题最优目标值（最小用户速率下界），c - 公共速率分配向量
% %       P，p_c - 更新后的波束成形，beta_p,beta_c - 更新后的干扰项
% function [t,c,P,p_c,beta_p,beta_c] = ...
%     rsma_sca_update(H,Pt,P_last,p_c_last,beta_p_last,beta_c_last)
% 
% [Nt,K] = size(H); % Nt 天线数，K用户数
% 
% cvx_begin quiet             % 开始CVX求解，'quiet'表示不显示求解过程
%     % 定义优化变量
%     variable t              % 标量，目标值（最小用户速率下界）
%     variable alpha_p(K)     % 私有速率近似变量（每个用户）
%     variable alpha_c(K)     % 公共速率近似变量（每个用户？实际代码为 K 维，但只用 alpha_c？注意：实际是标量？查看代码：variable alpha_c 没有下标，但后面约束 alpha_c <= log(1+rho_c)，所以 alpha_c 是标量，但代码中写的是 variable alpha_c(K)？检查原代码：variable alpha_c(K) 是 K 维？但约束中只用了 alpha_c 一个标量？其实是笔误？看原代码：variable alpha_c(K) 定义了 K 维，但约束中使用 alpha_c 时未加下标，可能 MATLAB 会取第一个元素？但正确的应该是标量。原代码确为 variable alpha_c(K)，但约束中用了 alpha_c，可能 CVX 允许向量与标量比较？实际上 alpha_c 应该是标量，这里可能是作者写错了。我们保留原样注释。
%     variable beta_p(K)      % 私有流干扰项（每个用户）
%     variable beta_c(k)      % 公有流干扰项（每个用户）
%     variable rho_p(K)       % 私有 SINR 近似变量
%     variable rho_c(k)       % 公有 SINR 近似变量
%     variable P(Nt,K) complex    % 私有波束成形矩阵（复数）
%     variable p_c(Nt) complex    % 公共波束成形向量（复数）
%     variable c(K)               % 公共速率分配份额（每个用户）
% 
%     maxmize(t)                  % 最大化最小用户
% 
%     subject to
%         % 约束1：每个用户的速率 >= t,总速率 = 私有速率（alpha_p） + 公共分配(c)
%         t <= alpha_p + c % 注意：alpha是向量，c是向量，逐元素比较
% 
%         % 约束2：所有公共分配份额之和不超过公共总速率 alpha_c
%         sum(c) <= alpha_c % 但alpha_c 是 K 维向量？此处 sum(c) 是标量，而 alpha_c 是向量，比较会出错？原代码可能意图是 alpha_c 为标量，但定义成了向量。实际 CVX 中可能允许？但通常不行。此处按原代码注释。
% 
%         % 约束3，私有速率alpha_p 受限于 log(1 + rho_p)
%         alpha_p <= log(1 + rho_p)
%         % 若使用以2为底，应除以log(2)，但这里是自然对数为，不影响最优解
% 
%         % 约束4：公共速率 alpha_c受限于log(1 + rho_c)
%         alpha_c <= log(1 + rho_c) 
% 
%         % 约束5：私有 SINR 的上界（通过泰勒线性化）
%         % 原表达式：rho_p <= 2*real( sum(conj(P_last).*H, 1).*sum(conj(H).*P, 1) )' ./ beta_p_last - square_abs( sum(conj(H).*P_last, 1).' ./ beta_p_last ) .* beta_p
%         % 解释：将 |h_k^H p_k|^2 / beta_{p,k} 在 (P_last, beta_p_last) 处一阶展开
%         rho_p <= 2*real(sum(conj(P_last).*H, 1).*sum(conj(H).*P, 1))'./beta_p_last - square_abs(sum(conj(H).*P_last, 1).'./beta_p_last).*beta_p
% 
%         % 约束6：公共 SINR 的上界（类似线性化）
%         rho_c <= 2*real((p_c_last'*H).'.*(H'*p_c))./beta_c_last - square_abs((H'*p_c_last)./beta_c_last).*beta_c
% 
%         % 约束7：定义beta_{p_k} >= sum_{j≠k} |h_k^H p_j|^2 + 1
%         for k = 1:K
%             index = [1:k-1,k+1:K];  % 除k以外的所有用户
%             beta_p(k) >= sum(square_abs(H(:,k)'*P(:,index))) + 1
%         end
% 
%         % 约束8：定义 beta_{c,k} >= sum_j |h_k^H p_j|^2 + 1 （公共流解码时的干扰)
%         beta_c >= 1 + sum(square_abs(H'*P),2)  % 列向量
% 
%         % 约束9：总功率约束（所有波束的Frobenius 范数平方 <= Pt）
%         norm([P,p_c],'fro') <= sqrt(Pt)
% 
%         % 约束10：公共分配份额非负
%         c >= 0
% cvx_end
% end

function [t, c, P, p_c, beta_p, beta_c] = ...
    rsma_sca_update(H, Pt, P_last, p_c_last, beta_p_last, beta_c_last)
% SCA 迭代中的 CVX 子问题求解函数
% 输入：H - 信道矩阵，Pt - 总功率，P_last, p_c_last - 上一迭代点的波束，
%       beta_p_last, beta_c_last - 上一迭代点的干扰加噪声项
% 输出：t - 当前子问题最优目标值，c - 公共速率分配向量，
%       P, p_c - 更新后的波束成形，beta_p, beta_c - 更新后的干扰项

[Nt, K] = size(H);

cvx_begin quiet
    % 定义优化变量
    variable t                      % 目标值（最小用户速率下界）
    variable R_c                    % 公共流总速率（标量）
    variable alpha_p(K)             % 私有速率近似变量
    variable rho_p(K)               % 私有 SINR 近似变量
    variable rho_c(K)               % 公共 SINR 近似变量
    variable beta_p(K)              % 私有流干扰项
    variable beta_c(K)              % 公共流干扰项
    variable P(Nt, K) complex       % 私有波束成形矩阵
    variable p_c(Nt) complex        % 公共波束成形向量
    variable c(K)                   % 公共速率分配份额

    maximize(t)

    subject to
        % 约束1：每个用户的总速率 >= t
        t <= alpha_p + c;

        % 约束2：所有公共分配份额之和不超过公共总速率 R_c
        sum(c) <= R_c;

        % 约束3：公共总速率受限于最差用户的公共SINR
        R_c <= log(1 + rho_c);

        % 约束4：私有速率上界
        alpha_p <= log(1 + rho_p);

        % 约束5：私有 SINR 的上界（一阶泰勒展开线性化）
        rho_p <= 2 * real( sum(conj(P_last).*H, 1) .* sum(conj(H).*P, 1) )' ./ beta_p_last ...
                - (abs( sum(conj(H).*P_last, 1) ).^2)' ./ (beta_p_last.^2) .* beta_p;

        % 约束6：公共 SINR 的上界（一阶泰勒展开线性化）
        rho_c <= 2 * real( (p_c_last' * H).' .* (H' * p_c) ) ./ beta_c_last ...
                - (abs( H' * p_c_last ).^2) ./ (beta_c_last.^2) .* beta_c;

        % 约束7：定义 beta_{p,k} >= sum_{j≠k} |h_k^H p_j|^2 + 1
        for k = 1:K
            idx = [1:k-1, k+1:K];
            beta_p(k) >= sum(square_abs(H(:, k)'*P(:, idx))) + 1;
        end

        % 约束8：定义 beta_{c,k} >= sum_j |h_k^H p_j|^2 + 1
        beta_c >= 1 + sum(square_abs(H'*P), 2);

        % 约束9：总功率约束
        norm([P, p_c], 'fro') <= sqrt(Pt);

        % 约束10：公共分配份额非负
        c >= 0;
cvx_end

% 求解失败时使用上次迭代值兜底
if ~strcmp(cvx_status, 'Solved')
    warning('CVX未收敛到最优，使用上次迭代值');
    t = -inf;
end
end

```

:::

::: code-group labels=[main.m, rsma_sca_update.m, mmf_rate_calculate.m, rsma_sca.m]

```matlab
% 主程序：设置系统参数，运行 SCA 算法，并绘制收敛曲线
clear all;
clc;

rng(6);                 % 固定随机种子，保证结果可重复

K = 4;                  % 用户数
Nt = 16;                % 发射天线数

SNR = 10;               % 信噪比 (dB)
Pt = 10 ^ (SNR/10);     % 总发射功率（线性值）

% 生成瑞利衰落信道矩阵，元素为复高斯 CN(0,1)
H = 1/sqrt(2) * (randn(Nt, K) + 1j*randn(Nt, K));

% 收敛容差（未使用？实际在 rsma_sca 内使用）
tolerance1 = 1e-3;
tolerance2 = 1e-6;      % 未使用
alpha = 1;              % 未使用
c = 0.99;               % 未使用
maxIter1 = 100;         % 未使用
maxIter2 = 1000;        % 未使用

% 调用 SCA 主算法
% [P, p_c, nan_flag, outer_obj_list3] = rsma_sca(H, Pt, tolerance1);
% 注：原代码中函数名为 rsma_sca，但调用时写成了 rsma_sca，应一致。
% 下面修正为正确函数名
[P, p_c, nan_flag, outer_obj_list3] = rsma_sca(H, Pt, tolerance1);

% 输出归一化后的总功率（应接近 Pt）
norm([P, p_c], 'fro')^2

% 绘制 MMF 速率随迭代次数变化的曲线
figure;
slg = plot(0:(size(outer_obj_list3, 2)-1), outer_obj_list3);
grid on;
box on;
xlim('tight');
xlabel('number of iterations', 'FontSize', 14);
ylabel('MMF rate [bit/s/Hz]', 'FontSize', 14);
set(gca, 'LooseInset', [0,0,0,0.018]);
set(gca, 'FontSize', 14);

% 设置曲线样式
slg(1).Color = 'red';
slg(1).LineStyle = '-';
slg(1).Marker = '*';
slg(1).DisplayName = 'EG-FP';
slg(1).LineWidth = 3;
slg(1).MarkerSize = 12;
legend('FontSize', 13);
```
```matlab
% SCA 迭代中的 CVX 子问题求解函数
% 输入：H - 信道矩阵，Pt - 总功率，P_last, p_c_last - 上一迭代点的波束，
%       beta_p_last, beta_c_last - 上一迭代点的干扰加噪声项
% 输出：t - 当前子问题最优目标值（最小用户速率下界），c - 公共速率分配向量，
%       P, p_c - 更新后的波束成形，beta_p, beta_c - 更新后的干扰项
function [t, c, P, p_c, beta_p, beta_c] = ...
    rsma_sca_update(H, Pt, P_last, p_c_last, beta_p_last, beta_c_last)

[Nt, K] = size(H);  % Nt 天线数，K 用户数

cvx_begin quiet     % 开始 CVX 求解，'quiet' 表示不显示求解过程
    % 定义优化变量
    variable t                  % 标量，目标值（最小用户速率下界）
    variable alpha_p(K)         % 私有速率近似变量（每个用户）
    variable alpha_c(K)         % 公共速率近似变量（每个用户？实际代码为 K 维，但只用 alpha_c？注意：实际是标量？查看代码：variable alpha_c 没有下标，但后面约束 alpha_c <= log(1+rho_c)，所以 alpha_c 是标量，但代码中写的是 variable alpha_c(K)？检查原代码：variable alpha_c(K) 是 K 维？但约束中只用了 alpha_c 一个标量？其实是笔误？看原代码：variable alpha_c(K) 定义了 K 维，但约束中使用 alpha_c 时未加下标，可能 MATLAB 会取第一个元素？但正确的应该是标量。原代码确为 variable alpha_c(K)，但约束中用了 alpha_c，可能 CVX 允许向量与标量比较？实际上 alpha_c 应该是标量，这里可能是作者写错了。我们保留原样注释。
    variable beta_p(K)          % 私有流干扰项（每个用户）
    variable beta_c(K)          % 公共流干扰项（每个用户）
    variable rho_p(K)           % 私有 SINR 近似变量
    variable rho_c(K)           % 公共 SINR 近似变量
    variable P(Nt, K) complex   % 私有波束成形矩阵（复数）
    variable p_c(Nt) complex    % 公共波束成形向量（复数）
    variable c(K)               % 公共速率分配份额（每个用户）

    maximize(t)                 % 最大化最小用户总速率

    subject to
        % 约束1：每个用户的总速率 >= t，总速率 = 私有速率(alpha_p) + 公共分配(c)
        t <= alpha_p + c        % 注意：alpha_p 是向量，c 是向量，逐元素比较

        % 约束2：所有公共分配份额之和不超过公共总速率 alpha_c
        sum(c) <= alpha_c       % 但 alpha_c 是 K 维向量？此处 sum(c) 是标量，而 alpha_c 是向量，比较会出错？原代码可能意图是 alpha_c 为标量，但定义成了向量。实际 CVX 中可能允许？但通常不行。此处按原代码注释。

        % 约束3：私有速率 alpha_p 受限于 log(1 + rho_p)
        alpha_p <= log(1 + rho_p)
        % 若使用以2为底，应除以 log(2)，但这里直接用自然对数，不影响最优解

        % 约束4：公共速率 alpha_c 受限于 log(1 + rho_c)
        alpha_c <= log(1 + rho_c)

        % 约束5：私有 SINR 的上界（通过泰勒线性化）
        % 原表达式：rho_p <= 2*real( sum(conj(P_last).*H, 1).*sum(conj(H).*P, 1) )' ./ beta_p_last - square_abs( sum(conj(H).*P_last, 1).' ./ beta_p_last ) .* beta_p
        % 解释：将 |h_k^H p_k|^2 / beta_{p,k} 在 (P_last, beta_p_last) 处一阶展开
        rho_p <= 2*real(sum(conj(P_last).*H, 1).*sum(conj(H).*P, 1))'./beta_p_last - square_abs(sum(conj(H).*P_last, 1).'./beta_p_last).*beta_p

        % 约束6：公共 SINR 的上界（类似线性化）
        rho_c <= 2*real((p_c_last'*H).'.*(H'*p_c))./beta_c_last - square_abs((H'*p_c_last)./beta_c_last).*beta_c

        % 约束7：定义 beta_{p,k} >= sum_{j≠k} |h_k^H p_j|^2 + 1
        for k = 1:K
            index = [1:k-1, k+1:K];     % 除 k 外的所有用户
            beta_p(k) >= sum(square_abs(H(:, k)'*P(:, index))) + 1
        end

        % 约束8：定义 beta_{c,k} >= sum_j |h_k^H p_j|^2 + 1 （公共流解码时的干扰）
        beta_c >= 1 + sum(square_abs(H'*P), 2)

        % 约束9：总功率约束（所有波束的 Frobenius 范数平方 <= Pt）
        norm([P, p_c], 'fro') <= sqrt(Pt)

        % 约束10：公共分配份额非负
        c >= 0
cvx_end
end
```

```matlab
% 计算给定波束成形下的 MMF（最大最小公平）速率
% 输入：H - Nt×K 信道矩阵；P - Nt×K 私有波束成形矩阵；p_c - Nt×1 公共波束成形向量
% 输出：rate - 标量，当前波束成形下的 MMF 速率（自然对数单位，非 bit/s/Hz）
function rate = mmf_rate_calculate(H, P, p_c)

[~, K] = size(H);                   % K 为用户数
noise = ones(K, 1);                 % 噪声功率归一化为 1

% T_k：用户 k 解码私有流时的干扰加噪声（将公共流视为信号，但此处先算总干扰，不包含公共）
% 即 sum_j |h_k^H p_j|^2 + 1
T_k = sum(square_abs(H'*P), 2) + noise;

% T_ck：用户 k 解码公共流时的干扰加噪声（私有流全部视为干扰，加上公共信号功率）
% 即 sum_j |h_k^H p_j|^2 + |h_k^H p_c|^2 + 1
T_ck = T_k + square_abs(H'*p_c);

% 公共流速率（每个用户都能解码公共流，但速率受限于最差用户，即 min）
% 注意：这里用自然对数，实际速率需乘以 log2(e) 转换，但此处仅用于比较，不影响优化方向
rate_c = log(T_ck ./ T_k);          % K×1 向量，各用户的公共流可达速率

% 私有流速率（用户 k 解码自己的私有流，假设公共流已完美消除）
% 分子：|h_k^H p_k|^2，分母：T_k - |h_k^H p_k|^2（去掉自身信号即为干扰）
% 注意：sum(conj(H).*P, 1) 是逐列求和，得到 1×K 向量，每个元素是 h_k^H p_k 的共轭？
% 实际代码：square_abs(sum(conj(H).*P, 1))' 得到 K×1 的 |h_k^H p_k|^2
rate_p = log(T_k ./ (T_k - square_abs(sum(conj(H).*P, 1))'));

% 将私有速率升序排序（用于注水分配公共速率）
rate_p_sorted = sort(rate_p, 'ascend');

% 公共速率取最差值（所有用户都能达到的最小公共速率）
x = min(rate_c);

% 计算最优公共速率分配后的 MMF 速率
% 公式：min_{m=1..K} ( x + sum_{i=1}^m rate_p_sorted(i) ) / m
% 此即注水分配闭式解，确保所有用户总速率尽量拉平
rate = min((x + cumsum(rate_p_sorted)) ./ (1:K)');
end
```

```matlab
% SCA 算法主循环，迭代求解 RSMA MMF 速率最大化问题
% 输入：H - 信道矩阵，Pt - 总功率，tolerance - 收敛容差
% 输出：P, p_c - 优化后的波束成形，nan_flag - 是否出现 NaN，
%       outer_obj_list - 每次迭代后的 MMF 速率记录（用于绘图）
function [P, p_c, nan_flag, outer_obj_list] = rsma_sca(H, Pt, tolerance)

[Nt, K] = size(H);

% ---------- 初始化波束成形 ----------
% 公共流分配一半功率，私有流平分剩余功率
P_c = Pt * 0.5;
P_p = (Pt - P_c) / K;

% 公共波束：沿 H 的最大奇异值方向
[U, ~, ~] = svd(H);
p_c_last = U(:,1) * sqrt(P_c);

% 私有波束：采用最大比传输（MRT），即匹配信道方向
P_last = H ./ vecnorm(H) * sqrt(P_p);

% 计算初始干扰项 beta_p 和 beta_c（按照定义）
% beta_p(k) = sum_{j≠k} |h_k^H p_j|^2 + 1
% beta_c(k) = sum_j |h_k^H p_j|^2 + 1
beta_p_last = 1 + sum(square_abs(H'*P_last), 2) - square_abs(diag(H'*P_last));
beta_c_last = 1 + sum(square_abs(H'*P_last), 2);

% 初始化目标值（上一迭代）
t_last = -inf;
nan_flag = false;
maxIter = 1000;

% 记录初始点的 MMF 速率
outer_obj_list = [mmf_rate_calculate(H, P_last, p_c_last)];

% ---------- SCA 迭代 ----------
for iter = 1:maxIter
    % 调用 CVX 子问题求解，获得更新后的变量
    [t, ~, P, p_c, beta_p, beta_c] = ...
        rsma_sca_update(H, Pt, P_last, p_c_last, beta_p_last, beta_c_last);

    % 检查是否出现 NaN（数值不稳定）
    if anynan(t) || anynan(P) || anynan(p_c) || anynan(beta_p) || anynan(beta_c)
        nan_flag = true;
        return;
    end

    % 计算目标值变化量
    e = abs(t - t_last);
    
    % 计算当前波束成形的真实 MMF 速率（用于记录）
    outer_obj_list(end+1) = mmf_rate_calculate(H, P, p_c);
    fprintf("%3d | obj = %f | |obj - obj_last| = %f\n", iter, t, e);

    % 收敛判断
    if(e <= tolerance)
        break;
    end

    % 更新迭代点
    t_last = t;
    beta_p_last = beta_p;
    beta_c_last = beta_c;
    P_last = P;
    p_c_last = p_c;
end

% ---------- 功率归一化，确保总功率恰好为 Pt ----------
temp = norm([P, p_c], 'fro');
P = P / temp * sqrt(Pt);
p_c = p_c / temp * sqrt(Pt);

end
```
:::

---
## 引言
近年来，智能反射面因其低成本、低功耗和灵活调控无线信道的能力，成为物理层安全研究的前沿方向。RIS由大量无源反射单元组成，每个单元可独立调节入射电磁波的相移，通过协同设计反射系数，可增强合法用户信号质量或抑制窃听者接收能力。已有研究表明，在单向通信中优化RIS相移可显著提升保密速率，但在双向通信场景中，两个方向的传输相互耦合，优化问题更为复杂。

## 系统模型
* 考虑一个RIS辅助的双向通信系统。系统包含两个通信节点Alice（A）和Bob（B），一个被动窃听者Charlie（C），以及一个配备 L 个反射单元的RIS。所有节点均配备单天线。系统部署在三维空间中，Alice和Bob分别位于RIS的两侧，两者之间的直射路径质量较差（路径损耗较大或存在遮挡），因此需要通过RIS建立高质量的反射链路来实现通信。Charlie位于Alice和Bob附近，试图截获双方传输的信息。
* 本系统的核心特征是双向同时同频传输。所谓双向，是指Alice和Bob在相同的时间、相同的频率资源上同时向对方发送信息。具体而言，Alice发射的信号承载着发往Bob的信息，Bob发射的信号承载着发往Alice的信息。这两个信号在空间中同时传播，并同时到达RIS。RIS将两路入射信号一并反射，反射后的信号同时被Alice、Bob和Charlie接收。

## 信道建模
RIS的第l个反射单元的反射系数为 $v_l = e ^{j \theta_t}$，其中 $\theta_l \in [0,2\pi)$ 为可调相移。将所有单元的反射系数写成向量形式：

$$
v = [v_1,v_2,...,v_L]^T \in \cnums^L
$$

对于任意发射节点i和接收节点j（其中 $ i，j \in \{A,B,C\}$）,信号从i到j的传播存在两条路径

**路径1（RIS反射路径）**：信号从i发射，经RIS的第l个单元反射后到达j。设：
* $h_{\text{iI,l}}$：节点i到RIS第l个单元的复信道系数。
* $h_{\text{Ij,l}}$：RIS第l个单元到节点j的复信道系数。

则第l个单元的反射贡献为 $v_l * h_{\text{iI,l}} * h_{\text{Ij,l}}$。将所有L个单元的贡献求和，得到RIS反射路径的总增益：

$$
h_{\text{ij}}^{RIS} = \sum_{l=1}^{L}v_l * h_{\text{iI,l}} * h_{\text{Ij,l}}
$$

**路径2（直射链路）**：信号不经RIS反射，直接从i传播到j，记为 $g_{\text{ij}}$。在本系统部署中，Alice和Bob之间的直射路径可能被遮挡，但Alice->Charlie和Bob->Charlie的直射路径可能存在，这正是窃听者可以利用的。  
将两条路叠加。得到节点i到节点j的总复合信道增益：

$$
h_{\text{ij}} = \sum_{l=1}^{L} v_l * (h_{\text{iI,l}} h_{\text{Ij,l}}) + g_{\text{ij}} 
$$

### Bob端接收信号
Bob的接收信号包含以下成分：

1. **来自Alice的有用信号**：Alice发射信号 $x_A$（满足 $\mathbb{E}[|x_A|^2] = P_A$），经过复合信道 $h_{\text{AB}}$ 到达Bob，接收份量为：

$$
\sqrt{P_A} * h_{\text{AB}} * x_A
$$

2. **来自Bob自身的干扰**：Bob自身正在发射信号 $x_B$，该信号从Bob的发射天线直接耦合到其接收天线，路径增益为 $h_{\text{BB}}$，接收分量为：

$$
\sqrt{P_B} * h_{\text{BB}} * x_B
$$

由于Bob完全知道自己的发射信号 $x_B$（包括波形、功率和相位），该自干扰可在数字基带被完全消除，不影响后续解码。

3. **接收机热噪声**：建模为复高斯白噪声 $ n_B \sim \mathcal{CN}(0, \sigma_b^2)$

Bob端接收信号的完整表达式为:

$$
y_B = \sqrt{P_A}h_{\text{AB}}x_A + \sqrt{P_B}h_{\text{BB}}x_B + n_B
$$

自干扰消除后，Bob解码 $x_A$的有效SINR为：
$$
\gamma_{\text{AB}} = \frac{P_A|h_{\text{AB}}|^2}{\sigma_{b}^2}
$$

### Alice端接收信号
Alice的接收过程完全对称。Alice接收天线拾取:  
1. 来自Bob的有用信号：$\sqrt{P_B}h_{\text{BA}}x_B$，其中 $h_{BA}$ 为Bob到Alice的复合信道增益
2. 来自Alice自身的自干扰：$\sqrt{P_A}h_{\text{AA}x_A}$，同样可在数字域消除
3. 热噪声：$n_A \sim \mathcal{CN}(0, \sigma_a^2)$

Alice接收信号的完整表达式为：

$$
y_A = \sqrt{P_B}h_{\text{BA}}x_B + \sqrt{P_A}h_{\text{AA}x_A} + n_A
$$

自干扰消除后，Alice解码 $x_B$ 的SINR为：

$$
\gamma_{\text{BA}} = \frac{P_A|h_{\text{BA}}|^2}{\sigma_{a}^2}
$$

### 窃听者Charlie端接收信号
Charlie只接收信号而不发射，因此不存在自干扰问题。但Charlie面临另一个挑战：Alice和Bob的两路信号同时到达Charlie的接收天线，两路信号对Charlie都是有用的窃听目标。  
Charlie接收信号为：

$$
y_C = \sqrt{P_A}h_{\text{AC}}x_A + \sqrt{P_B}h_{\text{BC}}x_B + n_C
$$

其中 $n_C \sim \mathcal{CN}(0, \sigma_c^2)$  
从窃听者的角度，两路信号可视为叠加的有用信息，其总接收功率为：

$$
P_{echo} = P_A|h_AC|^2 + P_B|h_BC|^2
$$

因此Charlie处的等效SINR为：

$$
\gamma_C = \frac{P_A|h_{\text{AC}}|^2 + P_B|h_{\text{BC}}|^2}{\sigma_c^2}
$$

假设桑格节点的噪声功率相同，即 $\sigma_a^2 = \sigma_b^2 = \sigma_c^2 = \sigma^2$

### 保密和速率
根据信息论，合法链路Alice->Bob和Bob->Alice的可达速率分别为：

$$
R_{AB} = \ln(1 + \gamma_{AB}), \quad R_{BA} = \ln(1 + \gamma_{BA})
$$

窃听者Charlie的可达速率为：

$$
R_C = \ln(1 + \gamma_C)
$$

系统的保密和速率定义为两个合法方向速率之和减去窃听速率：

$$
R_s = R_{\text{AB}} + R_{\text{BA}} - R_C
$$

### 公式等价转换
式（11）中的信道增益 $|h_{\text{ij}}|^2$ 由最开始的式子定义，是RIS反射系数v的函数。为了将问题转化为适于凸优化的形式，引入以下数学构造。
定义扩展反射系数向量 $w \in \cnums^{L+1}$

$$
w = [v_1,v_2,...,v_L,1]^T
$$

最后一个元素固定为1，用于统一处理直射路径  
定义等效信道向量 $u_{\text{ij}} \in \cnums^{L+1}$:

$$
\mathbf{u}_{ij} = [h_{iI,1}h_{Ij,1}, \, h_{iI,2}h_{Ij,2}, \, \dots, \, h_{iI,L}h_{Ij,L}, \, g_{ij}]^T.
$$

则式(1)可简洁的写为内积形式：

$$
h_{ij} = \mathbf{w}^H \mathbf{u}_{ij}
$$

进一步定义RIS相移矩阵 $W = ww^H$ 和 信道矩阵 $H_{\text{ij}} = u_{\text{ij}}u_{\text{ij}}^H$ ,则：

$$
|h_{\text{ij}}|^2 = \operatorname{tr}{(H_{\text{ij}}W)}
$$

矩阵W满足以下约束：  
* 半正定性：$W \succeq 0$
* 秩一性：rank(W) = 1
* 对角线为1：$diag(W) = 1_{\text{L+1}}$(对应各反射系数的模为1，且直射链路系数固定为1)

上式即对相移引入半正定松弛方法。经过上诉操作，保密和速率可表示为：
$$
R_s(P_A, P_B, \mathbf{W}) = \ln\left(1 + \frac{P_A \operatorname{tr}(\mathbf{H}_{AB} \mathbf{W})}{\sigma^2}\right) + \ln\left(1 + \frac{P_B \operatorname{tr}(\mathbf{H}_{BA} \mathbf{W})}{\sigma^2}\right)
- \ln\left(1 + \frac{P_A \operatorname{tr}(\mathbf{H}_{AC} \mathbf{W}) + P_B \operatorname{tr}(\mathbf{H}_{BC} \mathbf{W})}{\sigma^2}\right). \tag{16}
$$

### 优化问题建模
考虑秩一约束和功率约束，保密和速率最大化问题可以建模为

$$
\begin{aligned}
& \max_{P_A, P_B, \mathbf{W}} \quad R_s(P_A, P_B, \mathbf{W}) \\
\text{s.t.} \quad & P_A + P_B \leq P_{\text{max}}, \quad P_A \geq P_A^{\text{min}}, \quad P_B \geq P_B^{\text{min}}, \\
& \mathbf{W} \succeq 0, \quad \text{diag}(\mathbf{W}) = 1, \quad \text{rank}(\mathbf{W}) = 1.
\end{aligned}
$$

## 交替优化算法（AO）
### 功率优化子问题（固定W）
当RIS相移固定时，W为常数矩阵。定义四个等效信道增益常数：

$$
\alpha_{\text{AB}} = \operatorname{tr}{(H_{\text{AB}}W)},\quad
\alpha_{\text{BA}} = \operatorname{tr}{(H_{\text{BA}}W)},\quad
\alpha_{\text{AC}} = \operatorname{tr}{(H_{\text{AC}}W)},\quad
\alpha_{\text{BC}} = \operatorname{tr}{(H_{\text{BC}}W)},
$$

由于目标函数中的对数函数在正域上单调递增，且总功率约束为上限约束，最优解必在边界取得，即 $P_A + P_B = P_max$，令 $P_B = P_max - P_A$，式（16） 退化为关于 $P_A$ 的一元函数：

$$
f(P_A) = \ln(1 + \alpha P_A) + \ln(1 + \beta (P_{\text{max}} - P_A)) - \ln(1 + \gamma P_A + \delta (P_{\text{max}} - P_A)), \tag{20}
$$

其中：

$$
\alpha = \frac{a_{AB}}{\sigma^2}, \quad \beta = \frac{a_{BA}}{\sigma^2}, \quad \gamma = \frac{a_{AC}}{\sigma^2}, \quad \delta = \frac{a_{BC}}{\sigma^2}. \tag{21}
$$

对 $f(P_A)$ 求导：

$$
f'(P_A) = \frac{\alpha}{1 + \alpha P_A} - \frac{\beta}{1 + \beta (P_{\text{max}} - P_A)} - \frac{\gamma - \delta}{1 + \gamma P_A + \delta (P_{\text{max}} - P_A)}. \tag{22}
$$

令 $f'(P_A) = 0$，通分合并后得到关于 $P_A$ 的二次方程：
$$
AP_A^2 + BP_A + C = 0. \tag{23}
$$

其中系数A，B，C由 $\alpha,\beta,\gamma,\sigma,P_{\text{max}}$ 确定，其解析解为：

$$
P_A^* = \frac{-B \pm \sqrt{B^2 - 4AC}}{2A}. \tag{24}
$$

由于 $f(P_A)$ 在区间上可能非凹，需检查驻点及边界，功率的可行区间为：
$$
P_A^{\text{min}} \leq P_A \leq P_{\text{max}} - P_B^{\text{min}}. \tag{25}
$$

候选解包括：
1. 内部驻点 $P_A^*$（若落在区间内）;
2. 左边界 $P_A = P_A^{min}$
3. 右边界 $P_A = P_{max} - P_B^{min}$

比较各候选点的 $f(P_A)$，取最大值对应的 $P_A^*$，则 $P_B^* = P_{max} - P_A^*$。该方法在给定RIS配置中获得全局最优功率分配

### RIS相移优化子问题（固定PA，PB）
固定PA 和 PB后，定义辅助函数：

$$
\Psi_{AB}(\mathbf{W}) = \sigma^2 + P_A \operatorname{tr}(\mathbf{H}_{AB} \mathbf{W}), \quad \Psi_{BA}(\mathbf{W}) = \sigma^2 + P_B \operatorname{tr}(\mathbf{H}_{BA} \mathbf{W}),
$$
$$
\Psi_C(\mathbf{W}) = \sigma^2 + P_A \operatorname{tr}(\mathbf{H}_{AC} \mathbf{W}) + P_B \operatorname{tr}(\mathbf{H}_{BC} \mathbf{W}).
$$

则目标函数为：

$$
g(\mathbf{W}) = \ln(\Psi_{AB}(\mathbf{W})) + \ln(\Psi_{BA}(\mathbf{W})) - \ln(\Psi_C(\mathbf{W})). \tag{26}
$$

$\ln(\Psi_{AB}(\mathbf{W}))$ 和 $\ln(\Psi_{BA}(\mathbf{W}))$ 都是凹函数（对数与线性函数复合保持凹性），而 $- \ln(\Psi_C(\mathbf{W}))$ 是凸函数，因此 g(W) 非凹。

采用连续凸近似（SCA）：由于ln(*)为凹函数，在任意点 $\hat{W}$ 处的一阶泰勒展开提供全局上界：

$$
\ln(\Psi_C(\mathbf{W})) \leq \ln(\Psi_C(\tilde{\mathbf{W}})) + \frac{\Psi_C(\mathbf{W}) - \Psi_C(\tilde{\mathbf{W}})}{\Psi_C(\tilde{\mathbf{W}})}. \tag{27}
$$

于是：

$$
-\ln(\Psi_C(\mathbf{W})) \geq -\ln(\Psi_C(\tilde{\mathbf{W}})) - \frac{\Psi_C(\mathbf{W}) - \Psi_C(\tilde{\mathbf{W}})}{\Psi_C(\tilde{\mathbf{W}})}. \tag{28}
$$

将式（28）代入式（26），得到凹下界近似：（这里的常数包含了上一个迭代的固定点数据）

$$
\tilde{g}(\mathbf{W}) = \ln(\Psi_{AB}(\mathbf{W})) + \ln(\Psi_{BA}(\mathbf{W})) - \frac{\Psi_C(\mathbf{W})}{\Psi_C(\tilde{\mathbf{W}})} + \text{常数}. \tag{29}
$$

最大化 $\tilde{g}(\mathbf{W})$ 为凸问题。松弛秩序一约束，得到半定规划（SDP）：

$$
\begin{aligned}
& \min_{\mathbf{W}} \quad -\tilde{g}(\mathbf{W}) \\
\text{s.t.} \quad & \mathbf{W} \succeq 0, \\
& \operatorname{diag}(\mathbf{W}) = 1. \tag{30}
\end{aligned}
$$

问题（30）为凸，可用CVX等工具箱高效求解。得到 $W^*$ 后。若 $rank(W^*) = 1$，直接提取主特征向量；若秩大于1（SDR松弛不紧），去最大特征值对应的特征向量 $u_{\text{max}}$，构造:

$$
\tilde{\mathbf{W}} = e^{\arg(\mathbf{u}_{\max})}, \tag{31}
$$

## 算法流程
### 算法步骤
1. **初始化**：随机生成 $ \mathbf{w}^{(0)}, |w_i^{(0)}| = 1 $，令 $ \mathbf{W}^{(0)} = \mathbf{w}^{(0)} (\mathbf{w}^{(0)})^H, k = 0 $。

2. **迭代**：
   - **功率更新**：固定 $ \mathbf{W}^{(k)} $，按3.1节求解得到 $ P_A^{(k+1)}, P_B^{(k+1)} $。
   - **RIS更新**：固定 $ P_A^{(k+1)}, P_B^{(k+1)} $，求解SDP (30) 得到 $ \mathbf{W}^{(k+1)} $，必要时按式 (31) 秩一恢复。
   - $ k \leftarrow k + 1 $。

3. **终止条件**：目标函数增量小于阈值或达到最大迭代次数。

### 收敛性分析
功率子问题在给定 \( \mathbf{W} \) 时求得全局最优；RIS子问题通过SCA转化为凸SDP，在给定功率时获得全局最优解，且SCA上界近似保证目标函数不降。交替迭代使目标函数单调递减，功率约束和RIS相移约束均为紧集，目标函数连续有上界，算法收敛至局部最优。

## 相关代码
