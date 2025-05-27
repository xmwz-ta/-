clc; clear;

%% 参数设置

numRuns = 1000;
L = 0.15;  % 阵元间距

element_nom = [0, 0, 0;
               L, 0, 0;
               0, L, 0];

target_true = [30; 40; -100];
c0 = 1500;  % 声速 [m/s]

% 最大误差等级
max_c_err   = 10;         % 声速误差上限 ±10 m/s
max_tof_err = 0.00005;       % 测时误差上限 ±50 us
max_base_err = 0.02;      % 阵元间距误差 ±2%
max_phase_err = deg2rad(10);  % 相位误差上限 10°

numLevels = 100;
levels = linspace(0, 1, numLevels);

%% 灵敏度分析
meanErr_c     = zeros(1, numLevels);
meanErr_tof   = zeros(1, numLevels);
meanErr_base  = zeros(1, numLevels);
meanErr_phase = zeros(1, numLevels);

for i = 1:numLevels
    errs = zeros(numRuns,1);
    
    % 声速误差
    sig = levels(i) * max_c_err;
    for k = 1:numRuns
        p_est = estimatePosition(target_true, element_nom, ...
            0, 0, 0, sig, c0);
        errs(k) = norm(p_est - target_true);
    end
    meanErr_c(i) = mean(errs);

    % 测时误差
    sig = levels(i) * max_tof_err;
    for k = 1:numRuns
        p_est = estimatePosition(target_true, element_nom, ...
            0, sig, 0, 0, c0);
        errs(k) = norm(p_est - target_true);
    end
    meanErr_tof(i) = mean(errs);

    % 阵元基线误差
    sig = levels(i) * max_base_err;
    for k = 1:numRuns
        p_est = estimatePosition(target_true, element_nom, ...
            sig, 0, 0, 0, c0);
        errs(k) = norm(p_est - target_true);
    end
    meanErr_base(i) = mean(errs);

    % 相位误差
    sig = levels(i) * max_phase_err;
    for k = 1:numRuns
        p_est = estimatePosition(target_true, element_nom, ...
            0, 0, sig, 0, c0);
        errs(k) = norm(p_est - target_true);
    end
    meanErr_phase(i) = mean(errs);
end


subplot(2,2,1);
plot(levels * max_c_err, meanErr_c, '-o'); title('声速误差');
xlabel('声速误差 (m/s)'); ylabel('定位误差 (m)');

subplot(2,2,2);
plot(levels * max_tof_err * 1e6, meanErr_tof, '-o'); title('测时误差');
xlabel('测时误差 (\mus)'); ylabel('定位误差 (m)');

subplot(2,2,3);
plot(levels * max_base_err * 100, meanErr_base, '-o'); title('阵元间距误差');
xlabel('阵元误差 (%)'); ylabel('定位误差 (m)');

subplot(2,2,4);
plot(rad2deg(levels * max_phase_err), meanErr_phase, '-o'); title('相位误差');
xlabel('相位误差 (°)'); ylabel('定位误差 (m)');

function p_est = estimatePosition(target, element_nom, ...
        sigma_base, sigma_tof, sigma_phase, sigma_c, c_nom)

    freq = 30e3;
    c = c_nom + sigma_c * randn();
    lambda = c / freq;
    L = 0.15;

    R_true = norm(target);
    u_true = target / R_true;


    element_act = element_nom + sigma_base * L * randn(size(element_nom));

    % 使用实际阵元位置计算基线向量
    b12 = element_act(2,:) - element_act(1,:);
    b13 = element_act(3,:) - element_act(1,:);

    % 计算真实相位差
    delta_phi_true_x = (2*pi / lambda) * dot(b12, u_true);
    delta_phi_true_y = (2*pi / lambda) * dot(b13, u_true);

    % 模拟测量误差
    delta_phi_meas_x = delta_phi_true_x + sigma_phase * randn();
    delta_phi_meas_y = delta_phi_true_y + sigma_phase * randn();

    % 使用阵元位置反演
    b12_nom = element_nom(2,:) - element_nom(1,:);
    b13_nom = element_nom(3,:) - element_nom(1,:);

    Lx_nom = norm(b12_nom);
    Ly_nom = norm(b13_nom);

    u_x = (lambda / (2*pi)) * delta_phi_meas_x / Lx_nom;
    u_y = (lambda / (2*pi)) * delta_phi_meas_y / Ly_nom;


    sign_z = sign(target(3));  % 确保正负方向一致
    u_z = sign_z * sqrt(max(0, 1 - u_x^2 - u_y^2));

    u_est = [u_x; u_y; u_z];


    % 测距误差
    t_true = R_true / c_nom;
    t_meas = t_true + sigma_tof * randn();
    R_meas = c * t_meas;

    % 位置估计
    p_est = R_meas * u_est;
end

