clear; clc; close all;

[fitting_params, options] = loadFittingParameters();
[time_data, strain_data] = loadStrainData();
strain_data = strain_data / 100;
% displayInputInfo(time_data, strain_data);
tic;
predicted_stress = predictStressResponse(time_data, strain_data, fitting_params, options);
elapsed_time = toc;

savePredictionResults(time_data, strain_data, predicted_stress);
plotPredictionResults(time_data, strain_data, predicted_stress, fitting_params, options);

% analyzeHysteresis(time_data, strain_data, predicted_stress);

%% 加载
function [params, options] = loadFittingParameters()
    param_files = {'fitted_parameters.csv'};    
    params_loaded = false;    
    for i = 1:length(param_files)
        if exist(param_files{i}, 'file')          
            if endsWith(param_files{i}, '.mat')
                try
                    data = load(param_files{i});
                    param_names = {'optimized_params', 'params', 'fitted_params', 'optimized_parameters', 'parameters'};                    
                    for j = 1:length(param_names)
                        if isfield(data, param_names{j})
                            params = data.(param_names{j});
                            params_loaded = true;
                            break;
                        end
                    end
                    opt_names = {'fitting_options', 'options', 'model_options'};
                    for j = 1:length(opt_names)
                        if isfield(data, opt_names{j})
                            options = data.(opt_names{j});
                            fprintf('加载选项: %s\n', opt_names{j});
                            break;
                        end
                    end
                    
                    if params_loaded
                        break;
                    end
                    
                catch ME
                    fprintf('加载MAT文件失败: %s\n', ME.message);
                end
                
            elseif endsWith(param_files{i}, '.csv')
                try
                    param_table = readtable(param_files{i});
                    if all(ismember({'Parameter', 'Value'}, param_table.Properties.VariableNames))
                        param_values = param_table.Value;
                        params = param_values';
                        params_loaded = true;
                    end
                    
                catch ME
                    fprintf('加载CSV文件失败: %s\n', ME.message);
                end
            end
        end
    end
    
    if ~params_loaded
        params = inputManualParameters();
        options = struct();
        options.num_maxwell_elements = 3;
        options.ogden_order = 2;
    else
        if ~exist('options', 'var')
            options = struct();
            options.num_maxwell_elements = 3;  
            options.ogden_order = 2;           
        end
    end
end

%% 若无参数可输入修改
function params = inputManualParameters()
    fprintf('\n--- 手动输入模型参数 ---\n');
    n_params = input('请输入参数总数: ');
    params = zeros(1, n_params);
    for i = 1:n_params
        param_name = input(sprintf('参数 %d 名称: ', i), 's');
        param_value = input(sprintf('参数 %d 值: ', i));
        params(i) = param_value;
        
        if i == 1
            param_names = {param_name};
        else
            param_names{end+1} = param_name;
        end
    end

    save('manual_parameters.mat', 'params', 'param_names');
end

%% 读取时间应变数据
function [time, strain] = loadStrainData()
    data_files = {'strain_data.csv','experimental_data.csv'};
    file_found = false;    
    for i = 1:length(data_files)
        if exist(data_files{i}, 'file')
            try
                data = readtable(data_files{i});
                available_cols = lower(data.Properties.VariableNames);
                time_cols = {'time', 't', '时间', 'times'};
                time_idx = [];
                for j = 1:length(time_cols)
                    col_match = find(contains(available_cols, lower(time_cols{j})), 1);
                    if ~isempty(col_match)
                        time_idx = col_match;
                        break;
                    end
                end
                strain_cols = {'strain', 'epsilon', '应变', 'e', 'strain_pct'};
                strain_idx = [];
                for j = 1:length(strain_cols)
                    col_match = find(contains(available_cols, lower(strain_cols{j})), 1);
                    if ~isempty(col_match)
                        strain_idx = col_match;
                        break;
                    end
                end
                
                if ~isempty(time_idx) && ~isempty(strain_idx)
                    time = table2array(data(:, time_idx));
                    strain = table2array(data(:, strain_idx));
                    time = time(:);
                    strain = strain(:);
                    
                    file_found = true;
                    break;
                else
                end
                
            catch ME
                fprintf('读取文件 %s 失败: %s\n', data_files{i}, ME.message);
            end
        end
    end
    [time, strain] = preprocessStrainData(time, strain);
end

%% 应变数据处理 平滑、插值
function [time_clean, strain_clean] = preprocessStrainData(time_raw, strain_raw)
    time_raw = time_raw(:);
    strain_raw = strain_raw(:);

    if length(time_raw) ~= length(strain_raw)
        error('时间和应变数据长度不一致: time=%d, strain=%d', ...
            length(time_raw), length(strain_raw));
    end
    valid_idx = ~isnan(time_raw) & ~isnan(strain_raw);
    time_raw = time_raw(valid_idx);
    strain_raw = strain_raw(valid_idx);
    [time_sorted, sort_idx] = sort(time_raw);
    strain_sorted = strain_raw(sort_idx);
    if any(diff(time_sorted) <= 0)
        [time_unique, unique_idx] = unique(time_sorted);
        strain_unique = strain_sorted(unique_idx);
        N = min(1000, length(time_unique));  % 不超过1000点
        time_interp = linspace(min(time_unique), max(time_unique), N)';
        strain_interp = interp1(time_unique, strain_unique, time_interp, 'pchip');
        
        time_clean = time_interp;
        strain_clean = strain_interp;
    else

        N_original = length(time_sorted);
        if N_original > 1000
            time_interp = linspace(min(time_sorted), max(time_sorted), 1000)';
            strain_interp = interp1(time_sorted, strain_sorted, time_interp, 'pchip');
            
            time_clean = time_interp;
            strain_clean = strain_interp;
        else
            time_clean = time_sorted;
            strain_clean = strain_sorted;
        end
    end
end

%% 显示信息啊方便修改
function displayInputInfo(time, strain)

    fprintf('数据点数: %d\n', length(time));
    fprintf('时间范围: %.2f 到 %.2f 秒\n', min(time), max(time));
    fprintf('总时长: %.2f 秒\n', max(time) - min(time));

    if length(time) > 1
        time_diffs = diff(time);
        fprintf('时间步长范围: %.4f 到 %.4f 秒\n', min(time_diffs), max(time_diffs));
        fprintf('平均时间步长: %.4f 秒\n', mean(time_diffs));
    end
    fprintf('应变范围: %.1f%% 到 %.1f%%\n', min(strain)*100, max(strain)*100);
    fprintf('平均应变: %.1f%%\n', mean(strain)*100);

    if length(time) > 1 && length(strain) > 1
        strain_rate = diff(strain) ./ diff(time);
        max_strain_rate = max(abs(strain_rate));
        fprintf('最大应变率: %.3f s⁻¹\n', max_strain_rate);
        fprintf('平均应变率: %.3f s⁻¹\n', mean(abs(strain_rate)));
    else
        fprintf('最大应变率: 数据点不足\n');
        strain_rate = [];
    end
    fprintf('\n=== 应变统计 ===\n');
    fprintf('应变标准差: %.2f%%\n', std(strain)*100);
    if mean(strain) ~= 0
        fprintf('应变变异系数: %.2f%%\n', (std(strain)/abs(mean(strain)))*100);
    else
        fprintf('应变变异系数: 无穷大（平均应变为0）\n');
    end
    
    if ~isempty(strain_rate)
        threshold = 0.001;  % 0.1%每秒
        
        loading_count = sum(strain_rate > threshold);
        unloading_count = sum(strain_rate < -threshold);
        holding_count = sum(abs(strain_rate) <= threshold);
        
        total_count = length(strain_rate);
        
        fprintf('加载阶段点数: %d (%.1f%%)\n', loading_count, 100*loading_count/total_count);
        fprintf('卸载阶段点数: %d (%.1f%%)\n', unloading_count, 100*unloading_count/total_count);
        fprintf('保持阶段点数: %d (%.1f%%)\n', holding_count, 100*holding_count/total_count);
    else
        fprintf('数据点不足，无法分析加载卸载\n');
    end
end

%% 应力预测主函数
function stress = predictStressResponse(time, strain, params, options)
    n_ogden = options.ogden_order;
    n_maxwell = options.num_maxwell_elements;
    
    ogden_mu = params(1:n_ogden);
    ogden_alpha = params(n_ogden+1:2*n_ogden);
    maxwell_E = params(2*n_ogden+1:2*n_ogden+n_maxwell);
    maxwell_tau = params(2*n_ogden+n_maxwell+1:2*n_ogden+2*n_maxwell);
    E_eq = params(end);
    
    stress = calculateGeneralizedMaxwellStress(strain, time,ogden_mu, ogden_alpha, maxwell_E, maxwell_tau, E_eq);

    if any(isnan(stress)) || any(isinf(stress))
        stress(isnan(stress) | isinf(stress)) = 0;
    end
end

%% 根据模型计算应力
function stress = calculateGeneralizedMaxwellStress(strain, time,ogden_mu, ogden_alpha, maxwell_E, maxwell_tau, E_eq)   
    N = length(strain);
    stress = zeros(N, 1);
    
    if N == 0
        return;
    end
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
        stress_hyper = stress_hyper +(2 * ogden_mu(i) / ogden_alpha(i)) *(lambda.^ogden_alpha(i) - lambda.^(-ogden_alpha(i)/2));
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
                maxwell_stress(i, k) = maxwell_stress(i-1, k) * relaxation + ...
                    maxwell_E(k) * depsilon;
            end
        end
        stress_visco = sum(maxwell_stress(i, :));
        stress_eq = E_eq * strain(i);
        stress(i) = stress_hyper(i) + stress_visco + stress_eq;
    end
end

%% 保存结果
function savePredictionResults(time, strain, stress)
    results_table = table(time, strain*100, stress,'VariableNames', {'time_s', 'strain_pct', 'stress_MPa'});
    csv_filename = 'stress_prediction_results.csv';
    writetable(results_table, csv_filename);
%     mat_filename = 'stress_prediction_results.mat';
%     save(mat_filename, 'time', 'strain', 'stress');

    stats = struct();
    stats.time_min = min(time);
    stats.time_max = max(time);
    stats.strain_min = min(strain)*100;
    stats.strain_max = max(strain)*100;
    stats.stress_min = min(stress);
    stats.stress_max = max(stress);
    stats.stress_mean = mean(stress);
    stats.stress_std = std(stress);
    fprintf('  - CSV文件: %s\n', csv_filename);
end

%% ----画图----
function plotPredictionResults(time, strain, stress, params, options)
    n_ogden = options.ogden_order;
    n_maxwell = options.num_maxwell_elements;

    ogden_mu = params(1:n_ogden);
    ogden_alpha = params(n_ogden+1:2*n_ogden);
    maxwell_E = params(2*n_ogden+1:2*n_ogden+n_maxwell);
    maxwell_tau = params(2*n_ogden+n_maxwell+1:2*n_ogden+2*n_maxwell);
    
    figure (1); %
    set(figure (1),'position',[200,0,1280,500]);
    yyaxis left
    plot(time, strain*100, 'b-', 'LineWidth', 2);hold on
    xlabel('Time (s)', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('Strain (%)', 'FontSize', 11, 'FontWeight', 'bold');
    box on;
    xlim([-0.1,0.7]);
    ylim([-10,120]);

    figure (2); %
    set(figure (2),'position',[200,0,1080,600]);
    if length(time) > 1
        strain_rate = diff(strain) ./ diff(time);
        time_rate = time(1:end-1) + diff(time)/2;
        plot(time_rate, strain_rate, 'g-', 'LineWidth', 1.5);hold on
        xlabel('Time (s)', 'FontSize', 11, 'FontWeight', 'bold');
        ylabel('Strain rate (s⁻¹)', 'FontSize', 11, 'FontWeight', 'bold');
        box on;
    end
    
    figure (1); %
    set(figure (1),'position',[200,0,1280,500]);
    yyaxis right
%     plot(time, stress, 'r-', 'LineWidth', 2);hold on
t1=328;t2=655;ernum=5;er=ernum./0.1;
    time=[time(1:t1);time(t1+1:t2)];
% %     stress=[stress(1:t1)-ernum*time(1:t1).^0.5;  stress(t1+1:t2)-ernum*time(1:t2-t1).^0.5-(stress(t1+er+1)-stress(t1))];
%     stress=[stress(1:t1)+ernum*(time(1:t1)-0.5*time(t1)).^2-ernum*(0.5*time(t1)).^2;  stress(t1+1:t2)];
    stress(stress<0)=0;
    plot(time, stress, 'r-', 'LineWidth', 2);
    hold on
    xlabel('Time (s)', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('Stress (MPa)', 'FontSize', 11, 'FontWeight', 'bold');
    box on; 
    ylim([-0.5,3.5]);
    saveas(gcf, 'stress_prediction_results.png');
end
