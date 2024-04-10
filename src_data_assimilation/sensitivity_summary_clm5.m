clc;
clear;
%%
model_name = 'cesm2_clm5_cen_vr_v2'; % 'cesm2_clm5_mic_vr_v22';
nn_exp_name = 'exp_pc_cesm2_23';
time_domain = 'whole_time';

% server
% data_dir_output = '/GFPS8p/cess11/taof/datahub/ensemble/output_data/';
% data_dir_input = '/GFPS8p/cess11/taof/datahub/ensemble/input_data/';

% mac
data_dir_output = '/Users/ft254/DATAHUB/ENSEMBLE/OUTPUT_DATA/';
data_dir_input = '/Users/ft254/DATAHUB/ENSEMBLE/INPUT_DATA/';

%%
bootstrap_num = 200;

stat_matric_index = 3; % for NSE

component_name = {'H_NN', 'A_A', 'F_I', 'D_K', 'E_V', 'B_Xi', 'C_NPP', 'Ga_Homo', 'default'};

process_list =  {'A', 'K', 'Xi', 'V', 'I', 'NPP'};
process_label = {'A_A', 'D_K', 'B_Xi', 'E_V', 'F_I', 'C_NPP'};

process_description = {'CUE', 'Baseline Decomposition', 'Environmental Impacts', 'Vertical Transport', 'Input Allocation', 'Carbon Input'};
climate_list = {'A', 'B', 'C', 'D', 'E_all'};

manage_proposal = -0.1:0.02:0.1;
depth_name = {'C_30cm', 'B_100cm', 'A_200cm', 'D_full_depth'};

%% component control
control_summary = nan(length(depth_name)*length(component_name), 8);
colnames = {'var_soc', 'upper_soc', 'lower_soc', 'var_stat', 'upper_stat', 'lower_stat', 'depth', 'component'};

iquantile = 4;
ibootstrap = 1;
icomponent = 2;
counter = 1;
for iquantile = 1:4
    for icomponent = 1:length(component_name)
        control_matric_total_stock = nan(bootstrap_num, 1);
        control_matric_stat = nan(bootstrap_num, 1);

        for ibootstrap =1:bootstrap_num
            disp(['processing quantile ', num2str(iquantile), ' component ', num2str(icomponent), ' bootstrap  ', num2str(ibootstrap)]);
            % soc global stock
            % control_total_stock = load([data_dir_output, 'world_simulation_analyses/compontent_control_soc_total_stock_', model_name, '_', nn_exp_name, ...
            %     '_cross_valid_0_', num2str(ibootstrap), '_control_test.mat']);
            control_total_stock = load([data_dir_output, 'world_simulation_analyses/compontent_control_soc_total_stock_GCB_', model_name, '_', nn_exp_name, ...
                '_bootstrap_', num2str(ibootstrap), '_control_test.mat']);

            control_total_stock = control_total_stock.var_data_middle;
            % record
            control_matric_total_stock(ibootstrap) = control_total_stock(iquantile, icomponent);

            control_test = load([data_dir_output, 'world_simulation_analyses/compontent_control_spatial_variation_', model_name, '_', time_domain,...
                '_', nn_exp_name, '_bootstrap_', num2str(ibootstrap), '.mat']);

            control_test = control_test.var_data_middle; % [iquantile, icomponent, imatric]

            if icomponent == 1
                control_matric_stat(ibootstrap) = control_test(4, icomponent, stat_matric_index);
                valid_loc = find(control_matric_stat > 0);
            else
                control_matric_stat(ibootstrap) = control_test(4, icomponent, stat_matric_index) - control_test(4, 1, stat_matric_index);
            end
        end
        % summary
        control_summary(counter, 1) = median(control_matric_total_stock(valid_loc), 'omitnan');
        control_summary(counter, 2) = quantile(control_matric_total_stock(valid_loc), 0.86);
        control_summary(counter, 3) = quantile(control_matric_total_stock(valid_loc), 0.14);
        control_summary(counter, 4) = median(control_matric_stat(valid_loc), 'omitnan');
        control_summary(counter, 5) = quantile(control_matric_stat(valid_loc), 0.86);
        control_summary(counter, 6) = quantile(control_matric_stat(valid_loc), 0.14);


        control_summary(counter, 7) = iquantile;
        control_summary(counter, 8) = icomponent;

        counter = counter + 1;

    end
end

save([data_dir_output, 'world_simulation_analyses/compontent_control_summary_bootstrap_GCB_', model_name, '.mat'], 'control_summary');
