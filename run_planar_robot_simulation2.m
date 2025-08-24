function run_planar_robot_simulation()
    % 이 스크립트는 flexible_ref.pdf 논문의 2-link planar manipulator
    % 시뮬레이션을 재현합니다. InvDynFlex와 ForDyn을 모두 포함하며,
    % 제공된 파라미터를 기반으로 자유 진동 시나리오를 실행하고 결과를 플롯합니다.
    clc; clear all; close all;

    % --- 1. 시뮬레이션 파라미터 및 초기 조건 설정 ---
    disp('로봇 파라미터를 설정합니다...');
    robot_params = setup_robot_parameters();

    % 논문에 명시된 초기 조건 (자유 진동)
    q0 = [0; 0; 0; 0; 0.1; 0.002]; % [theta1; delta1_1; delta1_2; theta2; delta2_1; delta2_2]
    q_dot0 = zeros(6, 1);          % 초기 속도
    t_span = [0 2]; 
    
    disp('ODE 솔버를 사용하여 시뮬레이션을 실행합니다...');
    ode_func = @(t, y) dynamics_system(t, y, robot_params);
    options = odeset('RelTol', 1e-4, 'AbsTol', 1e-5, 'Stats', 'on');
    [T, Y] = ode45(ode_func, t_span, [q0; q_dot0], options);

    disp('시뮬레이션 완료.');

    % --- 3. 결과 시각화 ---
    disp('결과를 플롯합니다...');
    plot_results(T, Y);
end

function dydt = dynamics_system(t, y, robot_params)
    % ode45 솔버에 의해 호출될 메인 동역학 시스템 함수
    n_coords = size(y,1) / 2;
    q_f = y(1:n_coords);
    q_f_dot = y(n_coords+1:end);
    
    tau_input = zeros(n_coords, 1);
    q_f_ddot_zero = zeros(n_coords, 1);
    
    % InvDynFlex는 h+Kq 계산 및 ForDyn에 필요한 행렬들을 반환
    [h_Kq, stored_matrices] = InvDynFlex(q_f, q_f_dot, q_f_ddot_zero, robot_params);
    
    b = tau_input - h_Kq;
    
    % ForDyn을 사용하여 실제 가속도 계산
    q_f_ddot = ForDyn(b, robot_params, stored_matrices);
    
    dydt = [q_f_dot; q_f_ddot];
end

function params = setup_robot_parameters()
    % 논문의 모든 파라미터를 구조체로 정리합니다.
    params.n = 2;
    params.n_f = [2; 2];
    params.nx = [0; 0]; params.ny = [2; 2]; params.nz = [0; 0]; params.na = [0; 0];

    params.m = [0.5; 0.5]; 
    params.L = [0; 1; 1];   
    params.I_o = {diag([0,0,0.1]), diag([0,0,0.0083])};
    params.c = {[0.5;0;0], [0.5;0;0]};
    params.p_a = {[0;0;0], [0;0;0]};
    params.e_x = {[1;0;0], [1;0;0]};
    
    params.payload.M = diag([0.1, 0.1, 0.1, 0, 0, 0.0005]);

    phi_e = [0.1859, 0.2151; 0.8833, -0.0693];
    phi_prime_e = [0.6571, -0.5604; 2.6413, -10.8526];

    params.Phi_matrices = {[0, 0; phi_e(1,1), phi_e(1,2); 0, 0], [0, 0; phi_e(2,1), phi_e(2,2); 0, 0]};
    params.Delta_matrices = {[0, 0; 0, 0; phi_prime_e(1,1), phi_prime_e(1,2)], [0, 0; 0, 0; phi_prime_e(2,1), phi_prime_e(2,2)]};
    
    params.K_dd = {diag([0.91, 12.74]), diag([18.73, 999.88])};
    
    v_data = [0.0066, 0.0131; 0.0333, 0.0544];
    w_data = [0.0024, 0.0045; 0.0122, 0.0156];
    z_data_1 = [0.0007, 0.0013; 0.0013, 0.0024];
    z_data_2 = [0.0185, 0.0205; 0.0205, 0.0406];

    params.integral_matrices = struct;
    params.integral_matrices(1).V_u = [0, 0; v_data(1,:); 0, 0];
    params.integral_matrices(1).W_u = [0, 0; w_data(1,:); 0, 0];
    params.integral_matrices(1).Z_uu = z_data_1;
    params.integral_matrices(2).V_u = [0, 0; v_data(2,:); 0, 0];
    params.integral_matrices(2).W_u = [0, 0; w_data(2,:); 0, 0];
    params.integral_matrices(2).Z_uu = z_data_2;

    J_S_matrix = [2 0 0; 0 1 0; 0 0 1];
    I_S_matrix = diag([0, 1, 1]);
    for i=1:params.n
        params.integral_matrices(i).V_psi = zeros(3, params.n_f(i));
        params.integral_matrices(i).Y_psi_psi = zeros(params.n_f(i));
        params.integral_matrices(i).Z_uu_bar = zeros(params.n_f(i));
        params.integral_matrices(i).Y_psi_psi_bar = zeros(params.n_f(i));
        params.integral_matrices(i).Z_uu_bar_total = params.integral_matrices(i).Z_uu;
        params.J_S{i} = J_S_matrix;
        params.I_S{i} = I_S_matrix;
    end
end

function [tau_f, stored_matrices] = InvDynFlex(q_f, q_f_dot, q_f_ddot, robot_params)
    n = robot_params.n;
    n_f = robot_params.n_f;
    t = cell(1, n + 2); t_dot = cell(1, n + 2); w = cell(1, n + 2);
    tau = cell(1, n); 
    A_f_storage = cell(1, n); A_f_dot_storage = cell(1, n);
    P_store = cell(1,n); M_store = cell(1,n); A_store = cell(1, n + 2);
    R_im1_i_store = cell(1, n + 1);
    t{1} = zeros(6, 1); t_dot{1} = zeros(6, 1);
    
    coord_idx = 1;
    for i = 1:(n + 1)
        if i > 1
            num_coords_prev = 1 + n_f(i-1);
            q_prev_start_idx = coord_idx - num_coords_prev;
            q_prev = q_f(q_prev_start_idx : q_prev_start_idx + num_coords_prev - 1);
            delta_prev = q_prev(2:end);
        else; delta_prev = []; end

        if i <= n
            num_coords_i = 1 + n_f(i);
            q_i = q_f(coord_idx : coord_idx + num_coords_i - 1);
            q_i_dot = q_f_dot(coord_idx : coord_idx + num_coords_i - 1);
            q_i_ddot = q_f_ddot(coord_idx : coord_idx + num_coords_i - 1);
            theta_i = q_i(1); delta_i = q_i(2:end);
        else; theta_i = 0; q_i_dot = []; q_i_ddot = []; end
            
        R_im1_i = get_total_rotation_matrix(theta_i, delta_prev, robot_params, i);
        R_im1_i_store{i} = R_im1_i;
        
        if i == 1; m_f1 = 6 + n_f(1); A_f_i_im1 = [eye(6); zeros(m_f1-6, 6)];
        else; A_f_i_im1 = A_f_storage{i-1}; end
        
        if i <= n; R_ext = blkdiag(R_im1_i', R_im1_i', eye(n_f(i)));
        else; R_ext = blkdiag(R_im1_i', R_im1_i'); end
        A_i_im1 = R_ext * A_f_i_im1;
        A_store{i+1} = A_i_im1;
        
        t_prev = [t{i}; zeros(size(A_i_im1, 2) - size(t{i}, 1), 1)];
        if i <= n; P_i = get_P_matrix_f(robot_params, i); P_store{i} = P_i;
        else; P_i = []; end
        
        if i <= n; t{i+1} = A_i_im1 * t_prev + P_i * q_i_dot;
        else; t{i+1} = A_i_im1 * t_prev; end
        
        omega_i = t{i+1}(4:6);
        omega_skew = [0,-omega_i(3),omega_i(2); omega_i(3),0,-omega_i(1); -omega_i(2),omega_i(1),0];
        if i <= n; Omega_i = blkdiag(omega_skew, omega_skew, zeros(n_f(i)));
        else; Omega_i = blkdiag(omega_skew, omega_skew); end
        
        if i == 1; A_f_i_im1_dot = zeros(size(A_f_i_im1));
        else; A_f_i_im1_dot = A_f_dot_storage{i-1}; end
        A_i_im1_dot = R_ext * A_f_i_im1_dot;
        
        if i <= n
            P_dot_i = Omega_i * P_i;
            t_dot{i+1} = A_i_im1 * t_dot{i} + A_i_im1_dot * t{i} + P_i * q_i_ddot + P_dot_i * q_i_dot;
        else
            t_dot{i+1} = A_i_im1 * t_dot{i} + A_i_im1_dot * t{i};
        end
            
        if i <= n
            delta_i_dot = q_i_dot(2:end);
            omega_im1 = t{i}(4:6);
            [A_f0, A_f1] = get_A_matrix_f_components(delta_i, robot_params, i);
            A_f_storage{i} = A_f0 + A_f1;
            [A_dot_f0, A_dot_f1] = get_A_dot_matrix_f_components(omega_i, omega_im1, delta_i, delta_i_dot, robot_params, i);
            A_f_dot_storage{i} = A_dot_f0 + A_dot_f1;
        end
        if i <= n; coord_idx = coord_idx + num_coords_i; end
    end
    
    M_payload = robot_params.payload.M;
    omega_payload = t{n+2}(4:6);
    gamma_payload = [zeros(3,1); cross(omega_payload, M_payload(4:6,4:6) * omega_payload)];
    w{n+2} = M_payload * t_dot{n+2} + gamma_payload;
    
    coord_idx_end = size(q_f,1);
    for i = n:-1:1
        num_coords_i = 1 + n_f(i);
        q_i_start_idx = coord_idx_end - num_coords_i + 1;
        q_i = q_f(q_i_start_idx:coord_idx_end);
        q_i_dot = q_f_dot(q_i_start_idx:coord_idx_end);
        
        M_f_i = get_M_f(q_i, q_i_dot, robot_params, i);
        M_store{i} = M_f_i;
        gamma_f_i = get_gamma_f(t{i+1}, q_i, q_i_dot, M_f_i, robot_params, i);
        W_i = M_f_i * t_dot{i+1} + gamma_f_i;
        
        A_ip1_i = A_store{i+2};
        w_ip1_resized = [w{i+2}; zeros(size(A_ip1_i,1)-size(w{i+2},1),1)];
        w{i+1} = W_i + A_ip1_i' * w_ip1_resized;
        
        P_i = P_store{i};
        tau{i} = P_i' * w{i+1};
        coord_idx_end = q_i_start_idx - 1;
    end

    tau_f = vertcat(tau{:});
    stored_matrices.A = A_store;
    stored_matrices.P = P_store;
    stored_matrices.M = M_store;
end

function q_ddot = ForDyn(b, robot_params, stored_matrices)
    n = robot_params.n; n_f = robot_params.n_f;
    A = stored_matrices.A; P = stored_matrices.P; M = stored_matrices.M;
    M_payload = robot_params.payload.M;
    M_bar = cell(1, n + 2); phi_bar = cell(1, n + 1); Delta = cell(1, n);
    phi = cell(1, n); eta = cell(1, n + 2); X_hat = cell(1, n);
    X_bar = cell(1, n); X = cell(1, n); mu = cell(1, n + 1);

    eta{n + 2} = zeros(size(M_payload, 1), 1);
    M_bar{n + 2} = M_payload;

    coord_idx_end = size(b, 1);
    for i = n:-1:1
        num_coords_i = 1 + n_f(i);
        b_i_start_idx = coord_idx_end - num_coords_i + 1;
        b_i = b(b_i_start_idx:coord_idx_end);
        
        A_ip1_i = A{i+2};
        eta_i_ip1 = A_ip1_i' * eta{i + 2};
        X_hat{i} = b_i - P{i}' * eta_i_ip1;
        
        if i == n; phi_ip1_term = zeros(size(M_bar{i+2}));
        else; phi_ip1_term = phi{i+1} * phi_bar{i+1}'; end
        M_bar{i+1} = M{i} + A_ip1_i' * (M_bar{i+2} - phi_ip1_term) * A_ip1_i;

        phi_bar{i} = M_bar{i+1} * P{i};
        Delta{i} = P{i}' * phi_bar{i};
        phi{i} = phi_bar{i} / Delta{i};
        
        eta{i+1} = phi{i} * X_hat{i} + eta_i_ip1;
        coord_idx_end = b_i_start_idx - 1;
    end
    
    for i = 1:n; X_bar{i} = Delta{i} \ X_hat{i}; end

    mu{1} = zeros(6, 1);
    for i = 1:n
        A_i_im1 = A{i+1};
        mu_i_im1 = A_i_im1 * mu{i};
        X{i} = X_bar{i} - phi{i}' * mu_i_im1;
        mu{i+1} = P{i} * X{i} + mu_i_im1;
    end
    q_ddot = vertcat(X{:});
end

function plot_results(T, Y)
    n_coords = size(Y,2)/2; q = Y(:, 1:n_coords);
    figure('Name', 'Planar Robot Simulation Results');
    subplot(2,2,1); plot(T, q(:,2), 'r-', T, q(:,3), 'b-'); title('Flexible coordinates \delta_1'); xlabel('Time [s]'); ylabel('\delta [-]'); legend('\delta_{1,1}', '\delta_{1,2}'); grid on;
    subplot(2,2,2); plot(T, q(:,5), 'g-', T, q(:,6), 'k-'); title('Flexible coordinates \delta_2'); xlabel('Time [s]'); ylabel('\delta [-]'); legend('\delta_{2,1}', '\delta_{2,2}'); grid on;
    subplot(2,2,3); plot(T, rad2deg(q(:,1)), 'm-', T, rad2deg(q(:,4)), 'c-'); title('Joint Angles \theta'); xlabel('Time [s]'); ylabel('\theta [deg]'); legend('\theta_1', '\theta_2'); grid on;
    subplot(2,2,4); title('Effector Position (Not Implemented)'); xlabel('Time [s]'); ylabel('Position [m]'); grid on;
end

% --- 이하 모든 Helper 함수들은 이전 답변들과 동일 ---
function R_total = get_total_rotation_matrix(theta, delta_prev, params, i)
    R_nominal = get_rotation_matrix_rigid(theta, params, i);
    if i == 1; R_flex_correction_matrix = eye(3);
    else; R_flex_correction_matrix = get_rotation_matrix_flexible(delta_prev, params, i-1); end
    R_total = R_flex_correction_matrix * R_nominal; 
end

function R_rigid = get_rotation_matrix_rigid(theta, params, i)
    if i <= params.n; R_rigid = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];
    else; R_rigid = eye(3); end
end

function R_flex = get_rotation_matrix_flexible(delta_k, params, k)
    if isempty(delta_k) || k < 1; R_flex = eye(3); return; end
    Delta = params.Delta_matrices{k};
    psi = Delta * delta_k;
    th_x=psi(1); th_y=psi(2); th_z=psi(3);
    Rx = [1 0 0; 0 cos(th_x) -sin(th_x); 0 sin(th_x) cos(th_x)];
    Ry = [cos(th_y) 0 sin(th_y); 0 1 0; -sin(th_y) 0 cos(th_y)];
    Rz = [cos(th_z) -sin(th_z) 0; sin(th_z) cos(th_z) 0; 0 0 1];
    R_flex = Rz * Ry * Rx;
end

function [A_f0, A_f1] = get_A_matrix_f_components(delta_i, params, i)
    p_r_i = [params.L(i); 0; 0];
    p_r_i_skew = [0,-p_r_i(3),p_r_i(2); p_r_i(3),0,-p_r_i(1); -p_r_i(2),p_r_i(1),0];
    Phi_i = params.Phi_matrices{i};
    Delta_i = params.Delta_matrices{i};
    A_f0_top = [eye(3),-p_r_i_skew,Phi_i; zeros(3),eye(3),Delta_i];
    if i < params.n; num_pad_rows = params.n_f(i+1);
    else; num_pad_rows = 0; end
    A_f0 = [A_f0_top; zeros(num_pad_rows, size(A_f0_top, 2))];
    u_i_delta = Phi_i * delta_i;
    u_i_delta_skew = [0,-u_i_delta(3),u_i_delta(2); u_i_delta(3),0,-u_i_delta(1);-u_i_delta(2),u_i_delta(1),0];
    R_psi_i = get_rotation_matrix_flexible(delta_i, params, i);
    A_f1_top = [zeros(3),-u_i_delta_skew,zeros(size(Phi_i)); zeros(3),zeros(3),(R_psi_i - eye(3))*Delta_i];
    A_f1 = [A_f1_top; zeros(num_pad_rows, size(A_f1_top, 2))];
end

function [A_dot_f0, A_dot_f1] = get_A_dot_matrix_f_components(omega_i, omega_im1, delta_i, delta_i_dot, params, i)
    p_r_im1 = [params.L(i); 0; 0];
    Phi_im1 = params.Phi_matrices{i};
    Delta_im1 = params.Delta_matrices{i};
    omega_im1_skew = [0,-omega_im1(3),omega_im1(2); omega_im1(3),0,-omega_im1(1); -omega_im1(2),omega_im1(1),0];
    top_right_block = [omega_im1_skew * Phi_im1; omega_im1_skew * Delta_im1];
    omega_cross_p = cross(omega_im1, -p_r_im1);
    omega_cross_p_skew = [0,-omega_cross_p(3),omega_cross_p(2); omega_cross_p(3),0,-omega_cross_p(1); -omega_cross_p(2),omega_cross_p(1),0];
    A_dot_f0_6x = [zeros(3), omega_cross_p_skew, top_right_block(1:3, :); zeros(3), zeros(3), top_right_block(4:6, :);];
    if i < params.n; num_pad_rows = params.n_f(i+1); else; num_pad_rows = 0; end
    A_dot_f0 = [A_dot_f0_6x; zeros(num_pad_rows, size(A_dot_f0_6x, 2))];

    [~, A_f1] = get_A_matrix_f_components(delta_i, params, i);
    d_dt_Phi_delta = Phi_im1 * delta_i_dot + omega_im1_skew * Phi_im1 * delta_i;
    d_dt_Phi_delta_skew = [0,-d_dt_Phi_delta(3),d_dt_Phi_delta(2); d_dt_Phi_delta(3),0,-d_dt_Phi_delta(1); -d_dt_Phi_delta(2),d_dt_Phi_delta(1),0];
    psi_i = Delta_im1 * delta_i; psi_i_dot = Delta_im1 * delta_i_dot;
    R_psi_i = get_rotation_matrix_flexible(delta_i, params, i);
    R_psi_i_dot = get_R_psi_dot(psi_i, psi_i_dot);
    d_dt_R_minus_I_Delta = (R_psi_i_dot + omega_im1_skew * (R_psi_i - eye(3))) * Delta_im1;
    A_ring_f1_6x = [zeros(3), -d_dt_Phi_delta_skew, zeros(size(Phi_im1)); zeros(3), zeros(3), d_dt_R_minus_I_Delta];
    A_ring_f1 = [A_ring_f1_6x; zeros(num_pad_rows, size(A_ring_f1_6x, 2))];

    omega_i_skew = [0,-omega_i(3),omega_i(2); omega_i(3),0,-omega_i(1); -omega_i(2),omega_i(1),0];
    if i < params.n; Omega_i_op = blkdiag(omega_i_skew, omega_i_skew, zeros(params.n_f(i+1))); else; Omega_i_op = blkdiag(omega_i_skew, omega_i_skew); end
    Omega_bar_im1_op = blkdiag(omega_im1_skew, omega_im1_skew, zeros(params.n_f(i)));
    A_dot_f1 = A_ring_f1 + Omega_i_op * A_f1 - A_f1 * Omega_bar_im1_op;
end

function R_psi_dot = get_R_psi_dot(psi, psi_dot)
    sy = sin(psi(2)); cy = cos(psi(2)); sz = sin(psi(3)); cz = cos(psi(3));
    psi_dot_y = psi_dot(2); psi_dot_z = psi_dot(3);
    Mat1 = [-cz*sy, 0, 0; -sz*sy, 0, 0; -cy, 0, 0];
    Mat2 = [-sz*cy, -cz, 0; cz*cy, -sz, 0; 0, 0, 0];
    R_psi_dot = Mat1 * psi_dot_y + Mat2 * psi_dot_z;
end

function P = get_P_matrix_f(params, i)
    if i > params.n; P = []; return; end
    z_i = [0; 0; 1];
    P_rigid = [zeros(3,1); z_i];
    P = [P_rigid, zeros(6, params.n_f(i)); zeros(params.n_f(i), 1), eye(params.n_f(i))];
end

function M_f = get_M_f(q_i, q_i_dot, params, i)
    M_f = eye(6 + params.n_f(i)); % Placeholder
end

function gamma_f = get_gamma_f(t_i, q_i, q_i_dot, M_f_i, params, i)
    gamma_f = zeros(6 + params.n_f(i), 1); % Placeholder
end
