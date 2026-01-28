clear; clc; close all;
rng(42);
[strain_rates, data_cell] = loadMultiRateExperimentalDataNoTime();

fitting_options = struct();
fitting_options.num_maxwell_elements = 3;  
fitting_options.ogden_order = 2;           
fitting_options.display_iter = true;       
fitting_options.strain_rates = strain_rates; 

[optimized_params, fitted_stress_cell, fitting_metrics] = fitMultiRateGeneralizedMaxwellModelNoTime(data_cell, fitting_options);
displayMultiRateFittingResults(optimized_params, fitting_metrics, fitting_options);
plotMultiRateFittingResultsNoTime(data_cell, fitted_stress_cell,optimized_params, fitting_metrics, fitting_options);
validateMultiRateModelNoTime(optimized_params, fitting_options);


%% 数据读取
function [strain_rates, data_cell] = loadMultiRateExperimentalDataNoTime()
    strain_rates = [0.10, 0.20, 0.30, 0.40, 0.50, 1.00, 2.00, 3.00, 4.00, 5.00];
    rate_names = {'10', '20', '30', '40', '50', '100', '200', '300', '400', '500'};
    data_cell = cell(length(strain_rates), 1);
    data_loaded = false(length(strain_rates), 1);
    
    for i = 1:length(strain_rates)
        filename = sprintf('strain_rate_%s_percent.csv', rate_names{i});
        if exist(filename, 'file')
            try
                data = readtable(filename);

                available_cols = lower(data.Properties.VariableNames);
                strain_idx = find(contains(available_cols, 'strain'), 1);
                stress_idx = find(contains(available_cols, 'stress'), 1);
                if ~isempty(strain_idx) && ~isempty(stress_idx)
                    strain = table2array(data(:, strain_idx));
                    stress = table2array(data(:, stress_idx));

                    [strain, stress] = preprocessDataNoTime(strain, stress, strain_rates(i));

                    data_cell{i} = struct('strain', strain,'stress', stress, 'strain_rate', strain_rates(i));
                    data_loaded(i) = true;
                else
                end
            catch ME
            end
        else
            fprintf('未找到文件: %s\n', filename);
        end
    end
end

%% 数据预处理
function [strain_clean, stress_clean] = preprocessDataNoTime(strain_raw, stress_raw, strain_rate)
    strain_raw = strain_raw(:);
    stress_raw = stress_raw(:);
    valid_idx = ~isnan(strain_raw) & ~isnan(stress_raw);
    strain_raw = strain_raw(valid_idx);
    stress_raw = stress_raw(valid_idx);
    
    [strain_sorted, sort_idx] = sort(strain_raw);
    stress_sorted = stress_raw(sort_idx);

    window_size = min(7, floor(length(strain_sorted)/20));
    if window_size > 1
        strain_smooth = movmean(strain_sorted, window_size);
        stress_smooth = movmean(stress_sorted, window_size);
    else
        strain_smooth = strain_sorted;
        stress_smooth = stress_sorted;
    end
    
    N_interp = min(200, length(strain_smooth)); 
    if length(strain_smooth) > 1
        if min(strain_smooth) > 0
            strain_min = 0;
        else
            strain_min = min(strain_smooth);
        end
        
        strain_max = max(strain_smooth);
        strain_interp = linspace(strain_min, strain_max, N_interp)';

        stress_interp = interp1(strain_smooth, stress_smooth, strain_interp, 'spline');
    else
        strain_interp = strain_smooth;
        stress_interp = stress_smooth;
    end

    strain_clean = strain_interp;
    stress_clean = stress_interp;
    
end
%% 
function strain = generateLoadingUnloadingStrain(strain_rate)
    N = 200; 
    N_half = floor(N/2);
    strain_loading = linspace(0, 1.0, N_half)';
    strain_unloading = linspace(1.0, 0, N_half)';
    strain = [strain_loading; strain_unloading];
end

%% maxwell stress
function stress = calculateGeneralizedMaxwellStressNoTime(strain, strain_rate, ogden_mu, ogden_alpha, maxwell_E, maxwell_tau, E_eq)
    N = length(strain);
    stress = zeros(N, 1);
    if N == 0
        return;
    end
    time = generateTimeFromStrain(strain, strain_rate);
    if N > 1
        dt = diff(time);
        if any(dt <= 0)
            dt = ones(N-1, 1) * mean(abs(diff(time)));
        end
    else
        dt = 1;
    end

    maxwell_stress = zeros(N, length(maxwell_E));  % 各Maxwell单元应力
    lambda = 1 + strain;
    stress_hyper = zeros(N, 1);
    for i = 1:length(ogden_mu)
        if ogden_alpha(i) == 0
            continue;
        end
        stress_hyper = stress_hyper + (2 * ogden_mu(i) / ogden_alpha(i)) *(lambda.^ogden_alpha(i) - lambda.^(-ogden_alpha(i)/2));
    end
    for i = 1:N
        if i == 1
            maxwell_stress(i, :) = zeros(1, length(maxwell_E));
        else
            if i-1 <= length(dt)
                dt_i = dt(i-1);
            else
                dt_i = mean(dt);
            end
            if i > 1 && strain(i) ~= strain(i-1)
                depsilon = log(lambda(i)/lambda(i-1));  
            else
                depsilon = 0;
            end

            for k = 1:length(maxwell_E)
                relaxation = exp(-dt_i / maxwell_tau(k));
                maxwell_stress(i, k) = maxwell_stress(i-1, k) * relaxation + maxwell_E(k) * depsilon;
            end
        end
        stress_visco = sum(maxwell_stress(i, :));
        stress_eq = E_eq * strain(i);
        stress(i) = stress_hyper(i) + stress_visco + stress_eq;
    end
end

%% 时间向量
function time = generateTimeFromStrain(strain, strain_rate)
    N = length(strain);
    time = zeros(N, 1);
    [max_strain, max_idx] = max(strain);
    if max_idx > 1
        for i = 1:max_idx
            time(i) = strain(i) / strain_rate;
        end
    end
    if max_idx < N
        for i = (max_idx+1):N
            strain_decrease = max_strain - strain(i);
            time(i) = max_strain / strain_rate + strain_decrease / strain_rate;
        end
    end
end

%% 
function [optimized_params, fitted_stress_cell, metrics] = fitMultiRateGeneralizedMaxwellModelNoTime(data_cell, options)
   
    n_rates = length(data_cell);
    n_maxwell = options.num_maxwell_elements;
    n_ogden = options.ogden_order;
    show_iter = options.display_iter;

    n_params = 2*n_ogden + 2*n_maxwell + 1;  % Ogden参数 + Maxwell参数 + 平衡模量
 
    initial_guess = estimateMultiRateInitialParametersNoTime(data_cell, n_ogden, n_maxwell);
    lb = zeros(1, n_params) + 0.01;  
    ub = ones(1, n_params) * 10;     
    lb(n_ogden+1:2*n_ogden) = -5;    % ogden_alpha
    ub(n_ogden+1:2*n_ogden) = 5;     
    
    lb(end-n_maxwell+1:end-1) = 0.1; % 松弛时间
    ub(end-n_maxwell+1:end-1) = 20;  

    optim_options = optimoptions('fmincon', ...
        'Display', 'iter', ...
        'MaxFunctionEvaluations', 10000, ...
        'MaxIterations', 2000, ...
        'TolFun', 1e-8, ...
        'TolX', 1e-8, ...
        'Algorithm', 'interior-point');
    
    if ~show_iter
        optim_options.Display = 'final';
    end
    
    objective_func = @(params) calculateMultiRateFittingErrorNoTime(params, ...
        data_cell, n_ogden, n_maxwell);
    tic;
    [optimized_params, fval, exitflag, output] = fmincon(objective_func, ...
        initial_guess, [], [], [], [], lb, ub, [], optim_options);
    elapsed_time = toc;

    fitted_stress_cell = cell(n_rates, 1);
    metrics_cell = cell(n_rates, 1);
    
    for i = 1:n_rates
        data = data_cell{i};
        fitted_stress = calculateModelStressNoTime(optimized_params, ...
            data.strain, data.strain_rate, n_ogden, n_maxwell);
        fitted_stress_cell{i} = fitted_stress;

        metrics_cell{i} = calculateGoodnessOfFit(data.stress, fitted_stress);
    end

    metrics = calculateMultiRateGoodnessOfFit(data_cell, fitted_stress_cell);
    metrics.elapsed_time = elapsed_time;
    metrics.exitflag = exitflag;
    metrics.iterations = output.iterations;
    metrics.fval = fval;
    metrics.individual_metrics = metrics_cell;
end

%%  
function initial_params = estimateMultiRateInitialParametersNoTime(data_cell, n_ogden, n_maxwell)
    
    [~, max_idx] = max(cellfun(@(x) x.strain_rate, data_cell));
    data = data_cell{max_idx};
    
    strain = data.strain;
    stress = data.stress;
    [~, max_idx] = max(strain);
    if max_idx < length(strain) - 10
        end_strain = strain(end-10:end);
        end_stress = stress(end-10:end);
        p = polyfit(end_strain, end_stress, 1);
        E_eq_est = max(p(1), 0.05);
    else
        E_eq_est = 0.1;
    end
    ogden_mu_est = linspace(0.5, 0.1, n_ogden);
    ogden_alpha_est = linspace(1.5, -0.5, n_ogden);
    maxwell_E_est = linspace(0.8, 0.2, n_maxwell);
    maxwell_tau_est = logspace(-0.5, 1, n_maxwell); 
    initial_params = [ogden_mu_est, ogden_alpha_est, ...
        maxwell_E_est, maxwell_tau_est, E_eq_est];
end

%% stress 计算
function stress = calculateModelStressNoTime(params, strain, strain_rate, n_ogden, n_maxwell)
    ogden_mu = params(1:n_ogden);
    ogden_alpha = params(n_ogden+1:2*n_ogden);
    maxwell_E = params(2*n_ogden+1:2*n_ogden+n_maxwell);
    maxwell_tau = params(2*n_ogden+n_maxwell+1:2*n_ogden+2*n_maxwell);
    E_eq = params(end);

    stress = calculateGeneralizedMaxwellStressNoTime(strain, strain_rate, ogden_mu, ogden_alpha, maxwell_E, maxwell_tau, E_eq);
end

%% 
function error = calculateMultiRateFittingErrorNoTime(params, data_cell, n_ogden, n_maxwell)   
    total_error = 0;
    n_rates = length(data_cell);
    
    for i = 1:n_rates
        data = data_cell{i};
        stress_pred = calculateModelStressNoTime(params, data.strain, data.strain_rate, n_ogden, n_maxwell);
        basic_error = sum((stress_pred - data.stress).^2);

        weight = length(data.stress) / 1000; % 归一化
        total_error = total_error + weight * basic_error;
    end
    regularization = 0.01 * sum(params.^2);
    maxwell_tau = params(2*n_ogden+n_maxwell+1:2*n_ogden+2*n_maxwell);
    tau_penalty = sum(max(0, -maxwell_tau).^2) * 100;
        error = total_error + regularization + tau_penalty;
end

%% R2
function metrics = calculateGoodnessOfFit(exp_data, fit_data)
    valid_idx = ~isnan(exp_data) & ~isnan(fit_data);
    exp_clean = exp_data(valid_idx);
    fit_clean = fit_data(valid_idx);
    
    if isempty(exp_clean) || isempty(fit_clean)
        metrics.R2 = 0;
        metrics.RMSE = Inf;
%         metrics.MAE = Inf;
%         metrics.MAPE = Inf;
        return;
    end
    
    n = length(exp_clean);

    SS_res = sum((exp_clean - fit_clean).^2);
    SS_tot = sum((exp_clean - mean(exp_clean)).^2);
    if SS_tot == 0
        R2 = 1;
    else
        R2 = max(0, 1 - SS_res/SS_tot);
    end
    
    RMSE = sqrt(mean((exp_clean - fit_clean).^2));
%     MAE = mean(abs(exp_clean - fit_clean));
%     MAPE = mean(abs((exp_clean - fit_clean) ./ (exp_clean + eps))) * 100;
    
    metrics = struct();
    metrics.R2 = R2;
    metrics.RMSE = RMSE;
%     metrics.MAE = MAE;
%     metrics.MAPE = MAPE;
    metrics.SS_res = SS_res;
    metrics.SS_tot = SS_tot;
    metrics.n_points = n;
end

%% 
function metrics = calculateMultiRateGoodnessOfFit(data_cell, fitted_stress_cell)
    n_rates = length(data_cell);
    
    all_exp = [];
    all_fit = [];
    
    for i = 1:n_rates
        data = data_cell{i};
        fit_stress = fitted_stress_cell{i};
        valid_idx = ~isnan(data.stress) & ~isnan(fit_stress);
        all_exp = [all_exp; data.stress(valid_idx)];
        all_fit = [all_fit; fit_stress(valid_idx)];
    end
    metrics = calculateGoodnessOfFit(all_exp, all_fit);
    metrics.n_rates = n_rates;
end

%% 
function displayMultiRateFittingResults(params, metrics, options)
    n_ogden = options.ogden_order;
    n_maxwell = options.num_maxwell_elements;
    strain_rates = options.strain_rates;    

    ogden_mu = params(1:n_ogden);
    ogden_alpha = params(n_ogden+1:2*n_ogden);
    maxwell_E = params(2*n_ogden+1:2*n_ogden+n_maxwell);
    maxwell_tau = params(2*n_ogden+n_maxwell+1:2*n_ogden+2*n_maxwell);
    E_eq = params(end); 
    E_total = sum(maxwell_E) + E_eq;
    
%     % 显示总体拟合优度
%     fprintf('\n[总体拟合优度]\n');
%     fprintf('  决定系数 R² = %.6f\n', metrics.R2);
%     fprintf('  均方根误差 RMSE = %.6f MPa\n', metrics.RMSE);
%     fprintf('  平均绝对误差 MAE = %.6f MPa\n', metrics.MAE);
%     fprintf('  平均绝对百分比误差 MAPE = %.2f%%\n', metrics.MAPE);
%     fprintf('  拟合应变率数量: %d\n', metrics.n_rates);
    
    % 显示各应变率的拟合优度
    if isfield(metrics, 'individual_metrics')
%         fprintf('\n[各应变率拟合优度]\n');
        for i = 1:length(strain_rates)
            rate_metrics = metrics.individual_metrics{i};
%             fprintf('  应变率 %.1f%%/s: R²=%.4f, RMSE=%.4f MPa, MAPE=%.1f%%\n', ...
%                 strain_rates(i)*100, rate_metrics.R2, rate_metrics.RMSE, rate_metrics.MAPE);
        end
    end

    saveMultiRateParametersToFile(params, metrics, options);
end

%% save
function saveMultiRateParametersToFile(params, metrics, options)
    n_ogden = options.ogden_order;
    n_maxwell = options.num_maxwell_elements;
    param_names = cell(1, length(params));
    param_values = zeros(length(params), 1);
    
    idx = 1;
    
    for i = 1:n_ogden
        param_names{idx} = sprintf('ogden_mu_%d', i);
        param_values(idx) = params(idx);
        idx = idx + 1;
    end
    
    for i = 1:n_ogden
        param_names{idx} = sprintf('ogden_alpha_%d', i);
        param_values(idx) = params(idx);
        idx = idx + 1;
    end

    for i = 1:n_maxwell
        param_names{idx} = sprintf('maxwell_E_%d', i);
        param_values(idx) = params(idx);
        idx = idx + 1;
    end
    
    for i = 1:n_maxwell
        param_names{idx} = sprintf('maxwell_tau_%d', i);
        param_values(idx) = params(idx);
        idx = idx + 1;
    end
    
    param_names{idx} = 'E_eq';
    param_values(idx) = params(idx);

    param_table = table(param_names', param_values, 'VariableNames', {'Parameter', 'Value'});
    writetable(param_table, 'multi_rate_fitted_parameters.csv');
    if isfield(metrics, 'individual_metrics')

        n_rates = length(options.strain_rates);
        rate_table = table('Size', [n_rates, 5], ...
            'VariableTypes', {'double', 'double', 'double', 'double', 'double'}, ...
            'VariableNames', {'StrainRate', 'R2', 'RMSE', 'MAE', 'MAPE'});
        
        for i = 1:n_rates
            rate_metrics = metrics.individual_metrics{i};
            rate_table.StrainRate(i) = options.strain_rates(i);
            rate_table.R2(i) = rate_metrics.R2;
            rate_table.RMSE(i) = rate_metrics.RMSE;
%             rate_table.MAE(i) = rate_metrics.MAE;
%             rate_table.MAPE(i) = rate_metrics.MAPE;
        end
        
        writetable(rate_table, 'multi_rate_fitting_metrics.csv');
    end
end

%% 画图
function plotMultiRateFittingResultsNoTime(data_cell,fitted_stress_cell,params, metrics, options)
    n_rates = length(data_cell);
    strain_rates = options.strain_rates;
    n_ogden = options.ogden_order;
    n_maxwell = options.num_maxwell_elements;

    ogden_mu = params(1:n_ogden);
    ogden_alpha = params(n_ogden+1:2*n_ogden);
    maxwell_E = params(2*n_ogden+1:2*n_ogden+n_maxwell);
    maxwell_tau = params(2*n_ogden+n_maxwell+1:2*n_ogden+2*n_maxwell);

    colors = jet(n_rates);
    
    figure (1);
    set(figure (1),'position',[200,0,1080,600]);
    hold on;
    
%     for i = 1:n_rates
   for i = 1:5:6
        data = data_cell{i};
        fit_stress = fitted_stress_cell{i};
        plot(data.strain*100, data.stress, 'o', 'Color', colors(i,:),'MarkerSize', 4, 'MarkerFaceColor', colors(i,:), 'DisplayName', sprintf('%.0f%%/s 实验', strain_rates(i)*100));
%         plot(data.strain*100, fit_stress, '-', 'Color', colors(i,:), 'LineWidth', 1.5, 'DisplayName', sprintf('%.0f%%/s 拟合', strain_rates(i)*100));
    end
    xlim([0,100]);
    ylim([-0.5,1]);
    xlabel('Strain (%)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Stress (MPa)', 'FontSize', 12, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 9, 'NumColumns', 2);
    grid on;
    box on;

    saveas(gcf, 'multi_rate_generalized_maxwell_fitting_results.png');
end

%% 
function validateMultiRateModelNoTime(params, options)

    n_ogden = options.ogden_order;
    n_maxwell = options.num_maxwell_elements;
    ogden_mu = params(1:n_ogden);
    ogden_alpha = params(n_ogden+1:2*n_ogden);
    maxwell_E = params(2*n_ogden+1:2*n_ogden+n_maxwell);
    maxwell_tau = params(2*n_ogden+n_maxwell+1:2*n_ogden+2*n_maxwell);
    E_eq = params(end);
%     test_strain_rates = [0.1;0.2;0.3;0.4;0.5;1;2;3;4;5];  
    test_strain_rates = [0.1;1];  
    colors = jet(length(test_strain_rates));
    
    figure (1);
%     set(figure (2),'position',[200,0,1080,600]);
    hold on;
    for i = 1:length(test_strain_rates)
        strain_rate = test_strain_rates(i);
        strain = generateLoadingUnloadingStrain(strain_rate);
        stress = calculateGeneralizedMaxwellStressNoTime(strain, strain_rate, ogden_mu, ogden_alpha, maxwell_E, maxwell_tau, E_eq);

        plot(strain*100, stress, '-', 'Color', colors(i,:),'LineWidth', 1.5, 'DisplayName', sprintf('%.3f s⁻¹', strain_rate));
    end
       xlabel('Strain (%)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Stress (MPa)', 'FontSize', 12, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 9, 'NumColumns', 2);
    grid on;
    box on;
    
    saveas(gcf, 'multi_rate_model_validation.png');

end