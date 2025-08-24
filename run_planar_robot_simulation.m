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
    options = odeset('RelTol', 1e-20, 'AbsTol', 1e-20, 'Stats', 'on');
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








% =========================================================================
% --- Helper Functions ---
% =========================================================================
function [M_f, M1, M2_dot1, M2_dot2] = get_M_f(q_i, q_i_dot, params, i)
    delta = q_i(2:end);
    delta_dot = q_i_dot(2:end);
    m = params.m(i); c = params.c{i}; I_o = params.I_o{i};
    V_u = params.integral_matrices(i).V_u;
    % beam : pa = [0 0 0] ex = [1 0 0]
    p_a = params.p_a{i}; e_x = params.e_x{i}; W_u = params.integral_matrices(i).W_u;
    V_psi = params.integral_matrices(i).V_psi; J_S = params.J_S{i};
    Z_uu = params.integral_matrices(i).Z_uu; Z_uu_bar = params.integral_matrices(i).Z_uu_bar;
    Y_psi_psi = params.integral_matrices(i).Y_psi_psi; Y_psi_psi_bar = params.integral_matrices(i).Y_psi_psi_bar; I_S = params.I_S{i};
    Z_uu_bar_total = params.integral_matrices(i).Z_uu_bar_total;

    H_ru = cross_mat(p_a)*V_u + cross_mat(e_x)*W_u; % Eq (B.40)
    H_rpsi = J_S * V_psi; % Eq (B.43)
    
    % --- Constant part M_f^(0)
    M_vv_0 = m * eye(3);
    M_vw_0 = -[0 -c(3) c(2); c(3) 0 -c(1); -c(2) c(1) 0];
    M_ww_0 = I_o;
    M_dd_0 = Z_uu_bar_total; % Assuming Z_psi_psi_bar is negligible or included
    
    M0 = blkdiag(M_vv_0, M_ww_0, M_dd_0);
    M0(1:3, 4:6) = M_vw_0;
    M0(4:6, 7:8) = H_ru+H_rpsi;

    % --- Linear part M_f^(1)
    M_vw_1 = -cross_mat(V_u * delta);
    I_ru = -(cross_mat(p_a)*cross_mat(V_u*delta) + cross_mat(V_u*delta)*cross_mat(p_a)) ...
           -(cross_mat(e_x)*cross_mat(W_u*delta) + cross_mat(W_u*delta)*cross_mat(e_x)); % Eq (B.39)
    I_rpsi = cross_mat(V_psi*delta)*J_S - J_S*cross_mat(V_psi*delta); % Eq (B.42)
    H_uu = get_H_uu(delta, Z_uu); % Eq (B.31)
    H_psipsi = get_H_psipsi(delta, Y_psi_psi, I_S); % Eq (B.35)
    M_ww_1 = I_ru + I_rpsi;
    M_wd_1 = H_uu + H_psipsi;
    % M_vd_1 = V_u; % Eq (B.38)
    
    M1 = zeros(size(M0));
    M1(1:3, 4:6) = M_vw_1; M1(4:6, 1:3) = M_vw_1';
    M1(4:6, 4:6) = M_ww_1;
    M1(4:6, 7:end) = M_wd_1; M1(7:end, 4:6) = M_wd_1';
    
    % --- Quadratic part M_f^(2)
    % delta가 mode 2인 경우
    delta_u_mat = get_delta_u_mat(delta, Z_uu);
    delta_psi_mat = get_delta_psi_mat(delta, Y_psi_psi);
    
    I_uu = (delta' * Z_uu_bar * delta) * eye(3) - delta_u_mat' * Z_uu * delta_u_mat; % Eq (B.29)
    I_psipsi = ((delta'*Y_psi_psi_bar*delta)*eye(3) - delta_psi_mat'*Y_psi_psi*delta_psi_mat)*I_S + J_S*(delta_psi_mat'*Y_psi_psi*delta_psi_mat); % Eq (B.33)
    
    M_ww_2 = I_uu + I_psipsi;
    M2 = zeros(size(M0));
    M2(4:6, 4:6) = M_ww_2;

    % M2_dot (delta_dot, delta)
    delta_dot_u_mat = get_delta_u_mat(delta_dot, Z_uu);
    delta_dot_psi_mat = get_delta_psi_mat(delta_dot, Y_psi_psi);
    
    I_dot1_uu = (delta_dot' * Z_uu_bar * delta) * eye(3) - delta_dot_u_mat' * Z_uu * delta_dot_u_mat; % Eq (B.29)
    I_dot1_psipsi = ((delta_dot'*Y_psi_psi_bar*delta)*eye(3) - delta_dot_psi_mat'*Y_psi_psi*delta_dot_psi_mat)*I_S + J_S*(delta_dot_psi_mat'*Y_psi_psi*delta_dot_psi_mat); % Eq (B.33)
    
    M_dot1_ww_2 = I_dot1_uu + I_dot1_psipsi;
    M2_dot1 = zeros(size(M0));
    M2_dot1(4:6, 4:6) = M_dot1_ww_2;

    % M2_dot (delta, delta_dot)
    I_dot2_uu = (delta' * Z_uu_bar * delta_dot) * eye(3) - delta_dot_u_mat' * Z_uu * delta_dot_u_mat; % Eq (B.29)
    I_dot2_psipsi = ((delta'*Y_psi_psi_bar*delta_dot)*eye(3) - delta_dot_psi_mat'*Y_psi_psi*delta_dot_psi_mat)*I_S + J_S*(delta_dot_psi_mat'*Y_psi_psi*delta_dot_psi_mat); % Eq (B.33)
    
    M_dot2_ww_2 = I_dot2_uu + I_dot2_psipsi;
    M2_dot2 = zeros(size(M0));
    M2_dot2(4:6, 4:6) = M_dot2_ww_2;

    M_f = M0 + M1 + M2;
end

function M_f_dot = get_M_f_dot(q_i, q_i_dot, params, i)
    delta = q_i(2:end); 
    delta_dot = q_i_dot(2:end);
    % Eq (B.53)
    [~, M1_dot, ~, ~] = get_M_f(q_i_dot, q_i_dot, params, i); % M_f^(1)(delta_dot)
    
    [~, ~, M2_d1, ~] = get_M_f(q_i, q_i_dot, params, i); % M_f^(2)(delta_dot, delta_dot) components needed
    I_uu_d1 = M2_d1(4:6,4:6);
    
    [~, ~, ~, M2_d2] = get_M_f(q_i, q_i_dot, params, i);
    I_uu_d2 = M2_d2(4:6,4:6);

    M2_term1 = zeros(size(M1_dot));
    M2_term2 = zeros(size(M1_dot));
    M2_term1(4:6, 4:6) = I_uu_d1; % Simplified M_f^(2)(delta_dot, delta)
    M2_term2(4:6, 4:6) = I_uu_d2; % Simplified M_f^(2)(delta, delta_dot)

    M_f_dot = M1_dot + M2_d1 + M2_d2;
end

function M_f_tilde = get_M_f_tilde(q_i, q_i_dot, omega_i, params, i)
    % Eq (B.64) 및 하위 수식들을 기반으로 M_f_tilde 행렬을 계산합니다.
    
    n_f_i = params.n_f(i);
    delta_i = q_i(2:end);
    delta_i_dot = q_i_dot(2:end);

    % --- 각 하위 블록 계산 ---
    
    % 1. M_wv_tilde 계산
    [~, M1, ~, ~] = get_M_f(q_i, q_i_dot, params, i); % M_f^(1)
    M_wv = M1(1:3, 4:6);
    
    [~, M1_dot, ~, ~] = get_M_f(q_i_dot, q_i_dot, params, i); % M_f^(1)(delta_dot)
    M_wv_dot = M1_dot(1:3, 4:6);
    
    omega_skew = [0, -omega_i(3), omega_i(2); 
                  omega_i(3), 0, -omega_i(1); 
                 -omega_i(2), omega_i(1), 0];
                 
    M_wv_tilde = -(M_wv_dot + omega_skew * M_wv - M_wv * omega_skew);

    % 2. M_dv_tilde 계산
    M_dv = params.integral_matrices(i).V_u'; % M_vd' = (P_ru)' = V_u'
    M_dv_tilde = -(-M_dv * omega_skew);

    % 3. M_dw_tilde 계산 (I_tilde 항들의 합)
    I_tilde_ru = get_I_tilde_ru(omega_i, params, i);
    I_tilde_rpsi = get_I_tilde_rpsi(omega_i, params, i);
    I_tilde_uu = get_I_tilde_uu(delta_i, omega_i, params, i);
    I_tilde_psipsi = get_I_tilde_psipsi(delta_i, omega_i, params, i);
    
    M_dw_tilde = I_tilde_ru + I_tilde_rpsi + I_tilde_uu + I_tilde_psipsi;

    % 4. M_dd_tilde 계산 (H_tilde 항들의 합)
    H_tilde_uu = get_H_tilde_uu(omega_i, params, i);
    H_tilde_psipsi = get_H_tilde_psipsi(omega_i, params, i);
    
    M_dd_tilde = H_tilde_uu + H_tilde_psipsi;

    % --- 최종 M_f_tilde 행렬 조립 (Eq. B.64) ---
    M_f_tilde = zeros(6 + n_f_i);
    M_f_tilde(4:6, 1:3) = M_wv_tilde;
    M_f_tilde(7:end, 1:3) = M_dv_tilde;
    M_f_tilde(7:end, 4:6) = M_dw_tilde;
    M_f_tilde(7:end, 7:end) = M_dd_tilde;
end

% =========================================================================
% --- M_f_tilde 계산을 위한 추가 Helper 함수들 ---
% =========================================================================
% NOTE: 아래 함수들은 Appendix B.3의 상세 수식(B.57 ~ B.63)을 구현합니다.
function I_ru = get_I_tilde_ru(omega, params, i)
    V_u = params.integral_matrices(i).V_u;
    W_u = params.integral_matrices(i).W_u;
    p_a = params.p_a{i};
    e_x = params.e_x{i};
    omega_skew = [0,-omega(3),omega(2); omega(3),0,-omega(1); -omega(2),omega(1),0];
    p_a_skew = [0,-p_a(3),p_a(2); p_a(3),0,-p_a(1); -p_a(2),p_a(1),0];
    e_x_skew = [0,-e_x(3),e_x(2); e_x(3),0,-e_x(1); -e_x(2),e_x(1),0];
    
    I_ru = -V_u' * omega_skew * p_a_skew - W_u' * omega_skew * e_x_skew;
end

function I_rpsi = get_I_tilde_rpsi(omega, params, i)
    V_psi = params.integral_matrices(i).V_psi;
    J_S = params.J_S{i};
    omega_skew = [0,-omega(3),omega(2); omega(3),0,-omega(1); -omega(2),omega(1),0];

    I_rpsi = V_psi' * omega_skew * J_S;
end

function I_tilde_uu = get_I_tilde_uu(delta, omega, params, i)
    % 식 (B.59) 구현
    Z_uu = params.integral_matrices(i).Z_uu;
    Z_uu_bar = params.integral_matrices(i).Z_uu_bar;
    
    % [delta]_u 행렬 구성
    delta_u_mat = get_delta_u_mat(delta, i);
    
    % 식 내부 항 계산
    term1 = Z_uu * delta_u_mat * omega;
    term1_u_mat = get_delta_u_mat(term1, i); % [Z_uu [δ]_u ω_o]_u
    
    term2 = Z_uu_bar * delta * omega';
    
    % 최종 I_tilde_uu 계산
    I_tilde_uu = term1_u_mat - term2;
end

function I_tilde_psipsi = get_I_tilde_psipsi(delta, omega, params, i)
    % 식 (B.60) 구현
    Y_psi_psi = params.integral_matrices(i).Y_psi_psi;
    Y_psi_psi_bar = params.integral_matrices(i).Y_psi_psi_bar;
    I_S = params.I_S{i};
    J_S = params.J_S{i};

    % [omega]_Delta, [delta]_psi 행렬 구성
    omega_Delta_mat = get_omega_Delta_mat(omega, params, i);
    delta_psi_mat = get_delta_psi_mat(delta, i);
    
    % 식 내부 항 계산
    term1 = (omega_Delta_mat * Y_psi_psi * delta_psi_mat - Y_psi_psi_bar * delta * omega') * I_S
    omega
    J_S_omega = J_S * omega;
    J_S_omega_Delta_mat = get_omega_Delta_mat(J_S_omega, params, i);
    J_S_omega_Delta_mat
    Y_psi_psi
    delta_psi_mat
    term2 = J_S_omega_Delta_mat * Y_psi_psi * delta_psi_mat * omega
    
    % 최종 I_tilde_psipsi 계산
    I_tilde_psipsi = term1 - term2;
end

% =========================================================================
% --- 위 함수들을 위한 추가 Helper 함수들 ---
% =========================================================================

% function delta_u_mat = get_delta_u_mat(delta_vec, params, i)
%     % 식 (B.30)에 따른 [δ]_u 행렬 구성
%     % delta_vec = [delta_x; delta_y; delta_z; delta_alpha]
%     % planar case에서는 delta_x, delta_z, delta_alpha가 0이므로 delta_vec은 delta_y가 됩니다.
%     % 여기서는 일반적인 경우를 가정하여 구성합니다.
% 
%     nx = params.nx(i); % Number of traction modes
%     ny = params.ny(i); % Number of y-bending modes
%     nz = params.nz(i); % Number of z-bending modes
% 
%     delta_x = delta_vec(1:nx);
%     delta_y = delta_vec(nx+1 : nx+ny);
%     delta_z = delta_vec(nx+ny+1 : nx+ny+nz);
% 
%     delta_u_mat = blkdiag(diag(delta_x), diag(delta_y), diag(delta_z));
% end
% 
% function delta_psi_mat = get_delta_psi_mat(delta_vec, params, i)
%     % 식 (B.34)에 따른 [δ]_ψ 행렬 구성
%     nx = params.nx(i); ny = params.ny(i); nz = params.nz(i); na = params.na(i);
% 
%     delta_y = delta_vec(nx+1 : nx+ny);
%     delta_z = delta_vec(nx+ny+1 : nx+ny+nz);
%     delta_alpha = delta_vec(nx+ny+nz+1 : end);
% 
%     % 행렬의 비대각선 요소를 채우기 위한 로직
%     delta_psi_mat = zeros(nx+ny+nz+na);
%     delta_psi_mat(nx+ny+1:nx+ny+nz, nx+1:nx+ny) = diag(delta_z);
%     delta_psi_mat(nx+1:nx+ny, nx+ny+1:nx+ny+nz) = diag(delta_y);
%     delta_psi_mat(1:nx, nx+ny+nz+1:end) = diag(delta_alpha);
% end

function omega_Delta_mat = get_omega_Delta_mat(omega_vec, params, i)
    % 식 (B.60) 주변에 정의된 [ω_o]_Δ 연산자 구현
    nx = params.nx(i); ny = params.ny(i); nz = params.nz(i); na = params.na(i);
    
    omega_x = omega_vec(1);
    omega_y = omega_vec(2);
    omega_z = omega_vec(3);
    
    omega_Delta_mat = blkdiag(zeros(nx), omega_z*eye(ny), omega_y*eye(nz), omega_x*eye(na));
end

function H_tilde_uu = get_H_tilde_uu(omega, params, i)
    % % 식 (B.62) 구현
    % 
    % % 필요한 적분 행렬들을 params 구조체에서 가져옵니다.
    % % Z_ij는 ∫ρS φ_i' φ_j ds 형태의 행렬들입니다.
    % Z_xy = params.integral_matrices(i).Z_xy;
    % Z_xz = params.integral_matrices(i).Z_xz;
    % Z_yx = params.integral_matrices(i).Z_yx;
    % Z_yz = params.integral_matrices(i).Z_yz;
    % Z_zx = params.integral_matrices(i).Z_zx;
    % Z_zy = params.integral_matrices(i).Z_zy;
    % 
    % wx = omega(1); wy = omega(2); wz = omega(3);
    % 
    % nx = size(Z_xy, 2); % traction 모드 개수
    % ny = size(Z_xy, 1); % y-bending 모드 개수
    % nz = size(Z_xz, 1); % z-bending 모드 개수
    % na = 0; % Torsion 모드 개수 (이 행렬에서는 사용되지 않음)

    % 식 (B.62)의 블록 행렬 구조에 따라 조립 (planar 2x2)
    H_tilde_uu = zeros(2,2);
end

function H_tilde_psipsi = get_H_tilde_psipsi(omega, params, i)
    % % 식 (B.63) 구현
    % 
    % % 필요한 적분 행렬들을 params 구조체에서 가져옵니다.
    % % Y_ij는 ∫ρ Δ_i' J_s^-1 Δ_j ds 형태의 행렬들입니다.
    % Y_ay = params.integral_matrices(i).Y_ay;
    % Y_az = params.integral_matrices(i).Y_az;
    % Y_ya = params.integral_matrices(i).Y_ya;
    % Y_za = params.integral_matrices(i).Y_za;
    % 
    % % 식 (B.63) 아래의 정의에 따라 가중된 각속도 omega_tilde 계산
    % I_S = params.I_S{i};
    % omega_tilde = I_S * omega;
    % w_tilde_y = omega_tilde(2);
    % w_tilde_z = omega_tilde(3);
    % 
    % nx = 0; % Traction 모드 (이 행렬에서는 사용되지 않음)
    % ny = size(Y_ay, 2); % y-bending 모드 개수
    % nz = size(Y_az, 2); % z-bending 모드 개수
    % na = size(Y_ay, 1); % torsion 모드 개수
    % 
    % % 식 (B.63)의 블록 행렬 구조에 따라 조립
    % planar 2 x 2
    H_tilde_psipsi = zeros(2,2);
end

function gamma_f = get_gamma_f(t_i, q_i, q_i_dot, M_f_i, params, i)
    % Eq (5.53)
    omega_i = t_i(4:6);
    omega_skew = cross_mat(omega_i);
    Omega_i = blkdiag(omega_skew, omega_skew, zeros(params.n_f(i)));
    
    M_f_dot_i = get_M_f_dot(q_i, q_i_dot, params, i);
    M_f_tilde_i = get_M_f_tilde(q_i, q_i_dot, omega_i, params, i);

    K_delta_delta = params.K_dd{i};
    K_f_x_i = [zeros(6,1); K_delta_delta * q_i(2:end)];
    
    E_v = blkdiag(zeros(3), eye(3), eye(params.n_f(i)));
    
    gamma_f = (M_f_dot_i + Omega_i * M_f_i + M_f_tilde_i) * E_v * t_i + K_f_x_i;
end

function c_mat = cross_mat(vec)
    c_mat = [0, -vec(3), vec(2); vec(3), 0, -vec(1); -vec(2), vec(1), 0];
end

function H_uu = get_H_uu(delta, Z_uu)
    % Eq (B.31)
    % Assuming Z_uu is a struct with fields Zxy, Zxz, Zyx, ...
    H_uu = zeros(3, 2); % Placeholder
end
function H_psipsi = get_H_psipsi(delta, Y_psipsi, I_S)
    % Eq (B.35)
    H_psipsi = zeros(3, 2); % Placeholder
end
function d_u_mat = get_delta_u_mat(delta, Z_uu)
    % Eq (B.30)
    d_u_mat = [zeros(2,1), delta, zeros(2,1)]; % Placeholder
end
function d_psi_mat = get_delta_psi_mat(delta, Y_psipsi)
    % Eq (B.34)
    d_psi_mat = [zeros(2,1), zeros(2,1), delta]; % Placeholder
end


% (이전 답변의 Kinetics Loop Helper 함수들은 여기에 포함)
% --- Helper Functions ---
function R_total = get_total_rotation_matrix(theta, delta_prev, params, i)
    R_nominal = get_rotation_matrix_rigid(theta, params, i);
    % <<<<<<<<<<<<<<<< 오류 수정된 부분 >>>>>>>>>>>>>>>>
    if i == 1
        R_total = R_nominal + zeros(3,3);
    else
        % 이전 링크(i-1)의 변형으로 인한 보정
        R_flex_correction_matrix = get_rotation_matrix_flexible(delta_prev, params, i-1);
        % <<<<<<<<<<<<<<<< 여기까지 >>>>>>>>>>>>>>>>
        R_total = R_nominal + (R_flex_correction_matrix-eye(3)) * R_nominal; 
    end
    
end

function R_rigid = get_rotation_matrix_rigid(theta, params, i)
    if i <= params.n
        R_rigid = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1]; % z축 회전만
    else 
        R_rigid = eye(3); 
    end
end

function R_flex = get_rotation_matrix_flexible(delta_k, params, k)
% k = i-1
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
n_f = params.n_f(i)
% p => a
% 혹시 계산 안 맞으면 p, phi, Delta, delta, omega index 확인해보기
    p_r_i = [params.L(i); 0; 0];
    Phi_i = params.Phi_matrices{i};
    Delta_i = params.Delta_matrices{i};
    omega_i_skew = [0,-omega_i(3),omega_i(2); omega_i(3),0,-omega_i(1); -omega_i(2),omega_i(1),0];
    
    top_right_block = [omega_i_skew * Phi_i; omega_i_skew * Delta_i];
    omega_cross_p = cross(omega_i, p_r_i);
    omega_cross_p_skew = [0,-omega_cross_p(3),omega_cross_p(2); omega_cross_p(3),0,-omega_cross_p(1); -omega_cross_p(2),omega_cross_p(1),0];
    
    A_dot_f0_6x = [zeros(3), omega_cross_p_skew, top_right_block(1:3, :);
                   zeros(3), zeros(3),         top_right_block(4:6, :);];
    % <<<<<<<<<<<<<<<< 2. 인덱스 오류 수정: i=n일 때 예외 처리 >>>>>>>>>>>>>>>>
    if i < params.n
        num_pad_rows = params.n_f(i+1);
    else % i == n일 때, 다음은 페이로드이므로 유연좌표 없음
        num_pad_rows = 0;
    end
    A_dot_f0 = [A_dot_f0_6x; zeros(num_pad_rows, 3+3+params.n_f(i))];
    % <<<<<<<<<<<<<<<< 여기까지 >>>>>>>>>>>>>>>>

    d_dt_Phi_delta = Phi_i * delta_i_dot + omega_i_skew * Phi_i * delta_i;
    d_dt_Phi_delta_skew = [0,-d_dt_Phi_delta(3),d_dt_Phi_delta(2); d_dt_Phi_delta(3),0,-d_dt_Phi_delta(1); -d_dt_Phi_delta(2),d_dt_Phi_delta(1),0];
    
    psi_i = Delta_i * delta_i;
    psi_i_dot = Delta_i * delta_i_dot;
    R_psi_i = get_rotation_matrix_flexible(delta_i, params, i);
    R_psi_i_dot = get_R_psi_dot(psi_i, psi_i_dot);
    
    d_dt_R_minus_I_Delta = (R_psi_i_dot + omega_i_skew * (R_psi_i - eye(3))) * Delta_i;
    
    A_dot_f1_6x = [zeros(3), -d_dt_Phi_delta_skew, zeros(3, params.n_f(i));
                    zeros(3), zeros(3),             d_dt_R_minus_I_Delta];
    A_dot_f1 = [A_dot_f1_6x ; zeros(num_pad_rows, 3+3+params.n_f(i))];

    A_dot_f1 = A_dot_f0+A_dot_f1

end

function R_psi_dot = get_R_psi_dot(psi, psi_dot)
    sy = sin(psi(2)); cy = cos(psi(2));
    sz = sin(psi(3)); cz = cos(psi(3));
    psi_dot_y = psi_dot(2);
    psi_dot_z = psi_dot(3);
    Mat1 = [-cz*sy, 0, 0; -sz*sy, 0, 0; -cy, 0, 0];
    Mat2 = [-sz*cy, -cz, 0; cz*cy, -sz, 0; 0, 0, 0];
    R_psi_dot = Mat1 * psi_dot_y + Mat2 * psi_dot_z;
end

function P = get_P_matrix_f(params, i)
    if i > params.n; P = []; return; end
    z_i = [0; 0; 1]; 
    P_rigid = [zeros(3,1); z_i]; % revolute, read 96 page. Zi reads simply [0,−sαi,cαi]⊤ in Ri−1, or [0,0,1]⊤ in Ri
    P = [P_rigid, zeros(6, params.n_f(i)); zeros(params.n_f(i), 1), eye(params.n_f(i))];
end
% =========================================================================
% --- EXAMPLE USAGE ---
% =========================================================================
% % 1. 로봇 파라미터 정의
% robot_params.n = 2;
% robot_params.n_f = [2; 2];
% robot_params.L = [1; 1];
% robot_params.m = [1;1];
% robot_params.c = {[0.5;0;0], [0.5;0;0]};
% robot_params.I_o = {diag([0.1,0.1,0.1]), diag([0.1,0.1,0.1])};
% robot_params.p_a = {[0;0;0], [0;0;0]}; % Rigid part vector of segment
% robot_params.e_x = {[1;0;0], [1;0;0]}; % Slender part axis
% robot_params.J_S = {diag([0.1,0.1,0.1]), diag([0.1,0.1,0.1])};
% robot_params.I_S = {diag([0,0.1,0.1]), diag([0,0.1,0.1])};
% robot_params.K_dd = {diag([1,1]), diag([1,1])};
% 
% % Integral Matrices (Pre-computed from Appendix B.2)
% V_u1 = rand(3,2); W_u1 = rand(3,2); V_psi1 = rand(3,2);
% Z_uu1.Zxy = rand(1,1); % Example field
% Z_uu_bar1 = rand(2,2); Z_uu_bar_total1 = rand(2,2);
% Y_psi_psi1.Yzy = rand(1,1); Y_psi_psi_bar1 = rand(2,2);
% robot_params.integral_matrices(1) = struct('V_u', V_u1, 'W_u', W_u1, 'V_psi', V_psi1, 'Z_uu', Z_uu1, 'Z_uu_bar', Z_uu_bar1, 'Z_uu_bar_total', Z_uu_bar_total1, 'Y_psi_psi', Y_psi_psi1, 'Y_psi_psi_bar', Y_psi_psi_bar1);
% 
% V_u2 = rand(3,2); W_u2 = rand(3,2); V_psi2 = rand(3,2);
% % ... and so on for link 2
% robot_params.integral_matrices(2) = struct('V_u', V_u2, 'W_u', W_u2, 'V_psi', V_psi2, ...);
% 
% robot_params.Phi_matrices = {rand(3,2), rand(3,2)};
% robot_params.Delta_matrices = {rand(3,2), rand(3,2)};
% robot_params.payload.M = diag([1 1 1, 0.1 0.1 0.1]);
% 
% % 2. 상태 변수 정의
% n_total_coords = 2 + sum(robot_params.n_f);
% q_f = rand(n_total_coords, 1);
% q_f_dot = rand(n_total_coords, 1);
% q_f_ddot = rand(n_total_coords, 1);
% 
% % 3. 함수 호출
% tau_f = InvDynFlex_AppendixB_Complete(q_f, q_f_dot, q_f_ddot, robot_params);
% disp(tau_f);