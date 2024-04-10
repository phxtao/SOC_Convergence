% clear;
% clc;

%% cesm2 Settings
% parallel setting
% delete(gcp('nocreate'));
% parpool(100);


cesm2_case_name = 'sasu_f05_g16_checked_step4';
start_year = 661;
end_year = 680;

model_name = 'cesm2_clm5_mic_vr_v22';

nn_exp_name = 'exp_pc_cesm2_23';

time_domain = 'whole_time'; % 'whole_time', 'before_1985', 'after_1985'

ibootstrap = 7;
% parfor ibootstrap = start_id:end_id % 1:200

nn_exp_name = ['exp_pc_cesm2_23_cross_valid_0_', num2str(ibootstrap)];
% nn_exp_name = ['exp_pc_cesm2_23_bootstrap_', num2str(ibootstrap)];

disp(['Processing ibootstrap: ', num2str(ibootstrap)]);

warning('off');
format long e;

% paths
% mac
% data_path = '/Users/ft254/DATAHUB/ENSEMBLE/';
% cd(['/Users/ft254/GitHub/ENSEMBLE/SRC_DA/src_', model_name, '/']);

% server
data_path = '/glade/work/taofeng/datahub/ensemble/'; % '/GFPS8p/cess11/taof/datahub/ensemble/';
cd(['/glade/u/home/taofeng/ensemble/src_da/src_', model_name, '/']);

%% set vertical soil pools
month_num = 12;
soil_cpool_num = 7;
soil_decom_num = 20;

%% load wosis data
env_info = load([data_path, 'input_data/wosis_2019_snap_shot/wosis_2019_snapshot_hugelius_mishra_env_info.mat'], 'EnvInfo'); % NPP in this file will be used
env_info = env_info.EnvInfo;
% layer_info: "profile_id, date, upper_depth, lower_depth, node_depth, soc_layer_weight, soc_stock, bulk_denstiy, is_pedo"
wosis_profile_info = ncread([data_path, 'input_data/wosis_2019_snap_shot/soc_profile_wosis_2019_snapshot_hugelius_mishra.nc'], 'soc_profile_info'); % wosis profile info
wosis_soc_info = ncread([data_path, 'input_data/wosis_2019_snap_shot/soc_data_integrate_wosis_2019_snapshot_hugelius_mishra.nc'], 'data_soc_integrate'); % wosis SOC info

%% soil depth information
% width between two interfaces
dz = [2.000000000000000E-002, 4.000000000000000E-002, 6.000000000000000E-002, ...
    8.000000000000000E-002, 0.120000000000000, 0.160000000000000, ...
    0.200000000000000, 0.240000000000000, 0.280000000000000, ...
    0.320000000000000, 0.360000000000000, 0.400000000000000, ...
    0.440000000000000, 0.540000000000000, 0.640000000000000, ...
    0.740000000000000, 0.840000000000000, 0.940000000000000, ...
    1.04000000000000, 1.14000000000000, 2.39000000000000, ...
    4.67553390593274, 7.63519052838329, 11.1400000000000, ...
    15.1154248593737]';

% depth of the interface
zisoi = [2.000000000000000E-002, 6.000000000000000E-002, ...
    0.120000000000000, 0.200000000000000, 0.320000000000000, ...
    0.480000000000000, 0.680000000000000, 0.920000000000000, ...
    1.20000000000000, 1.52000000000000, 1.88000000000000, ...
    2.28000000000000, 2.72000000000000, 3.26000000000000, ...
    3.90000000000000, 4.64000000000000, 5.48000000000000, ...
    6.42000000000000, 7.46000000000000, 8.60000000000000, ...
    10.9900000000000, 15.6655339059327, 23.3007244343160, ...
    34.4407244343160, 49.5561492936897]';

% depth of the node
zsoi = [1.000000000000000E-002, 4.000000000000000E-002, 9.000000000000000E-002, ...
    0.160000000000000, 0.260000000000000, 0.400000000000000, ...
    0.580000000000000, 0.800000000000000, 1.06000000000000, ...
    1.36000000000000, 1.70000000000000, 2.08000000000000, ...
    2.50000000000000, 2.99000000000000, 3.58000000000000, ...
    4.27000000000000, 5.06000000000000, 5.95000000000000, ...
    6.94000000000000, 8.03000000000000, 9.79500000000000, ...
    13.3277669529664, 19.4831291701244, 28.8707244343160, ...
    41.9984368640029]';

% depth between two node
dz_node = zsoi - [0; zsoi(1:end-1)];

%% input from cesm2
% cesm2 resolution
cesm2_resolution_lat = 180/384;
cesm2_resolution_lon = 360/576;
lon_grid = (-180 + cesm2_resolution_lon/2 : cesm2_resolution_lon : 180 - cesm2_resolution_lon/2)';
lat_grid = (90 - cesm2_resolution_lat/2 : -cesm2_resolution_lat : -90 + cesm2_resolution_lat/2)';

% nbedrock
cesm2_simu_nbedrock = load([data_path, 'input_data/cesm2_simu/spinup_ss/',...
    cesm2_case_name, '_cesm2_ss_4da_', num2str(start_year), '_', num2str(end_year),...
    '_nbedrock.mat'], 'var_record_monthly_mean');
cesm2_simu_nbedrock = cesm2_simu_nbedrock.var_record_monthly_mean;
% ALTMAX
cesm2_simu_altmax = load([data_path, 'input_data/cesm2_simu/spinup_ss/',...
    cesm2_case_name, '_cesm2_ss_4da_', num2str(start_year), '_', num2str(end_year),...
    '_ALTMAX.mat'], 'var_record_monthly_mean');
cesm2_simu_altmax = cesm2_simu_altmax.var_record_monthly_mean;
% ALTMAX_LASTYEAR
cesm2_simu_altmax_last_year = load([data_path, 'input_data/cesm2_simu/spinup_ss/',...
    cesm2_case_name, '_cesm2_ss_4da_', num2str(start_year), '_', num2str(end_year),...
    '_ALTMAX_LASTYEAR.mat'], 'var_record_monthly_mean');
cesm2_simu_altmax_last_year = cesm2_simu_altmax_last_year.var_record_monthly_mean;
% CELLSAND
cesm2_simu_cellsand = load([data_path, 'input_data/cesm2_simu/spinup_ss/',...
    cesm2_case_name, '_cesm2_ss_4da_', num2str(start_year), '_', num2str(end_year),...
    '_CELLSAND.mat'], 'var_record_monthly_mean');
cesm2_simu_cellsand = cesm2_simu_cellsand.var_record_monthly_mean;
% NPP
cesm2_simu_npp = load([data_path, 'input_data/cesm2_simu/spinup_ss/',...
    cesm2_case_name, '_cesm2_ss_4da_', num2str(start_year), '_', num2str(end_year),...
    '_NPP.mat'], 'var_record_monthly_mean');
cesm2_simu_npp = cesm2_simu_npp.var_record_monthly_mean;
% SOILPSI
cesm2_simu_soil_water_potnetial = load([data_path, 'input_data/cesm2_simu/spinup_ss/',...
    cesm2_case_name, '_cesm2_ss_4da_', num2str(start_year), '_', num2str(end_year),...
    '_SOILPSI.mat'], 'var_record_monthly_mean');
cesm2_simu_soil_water_potnetial = cesm2_simu_soil_water_potnetial.var_record_monthly_mean;
% TSOI
cesm2_simu_soil_temperature = load([data_path, 'input_data/cesm2_simu/spinup_ss/',...
    cesm2_case_name, '_cesm2_ss_4da_', num2str(start_year), '_', num2str(end_year),...
    '_TSOI.mat'], 'var_record_monthly_mean');
cesm2_simu_soil_temperature = cesm2_simu_soil_temperature.var_record_monthly_mean;
% W_SCALAR
cesm2_simu_w_scalar = load([data_path, 'input_data/cesm2_simu/spinup_ss/',...
    cesm2_case_name, '_cesm2_ss_4da_', num2str(start_year), '_', num2str(end_year),...
    '_W_SCALAR.mat'], 'var_record_monthly_mean');
cesm2_simu_w_scalar = cesm2_simu_w_scalar.var_record_monthly_mean;
% T_SCALAR
cesm2_simu_t_scalar = load([data_path, 'input_data/cesm2_simu/spinup_ss/',...
    cesm2_case_name, '_cesm2_ss_4da_', num2str(start_year), '_', num2str(end_year),...
    '_T_SCALAR.mat'], 'var_record_monthly_mean');
cesm2_simu_t_scalar = cesm2_simu_t_scalar.var_record_monthly_mean;
% O_SCALAR
cesm2_simu_o_scalar = load([data_path, 'input_data/cesm2_simu/spinup_ss/',...
    cesm2_case_name, '_cesm2_ss_4da_', num2str(start_year), '_', num2str(end_year),...
    '_O_SCALAR.mat'], 'var_record_monthly_mean');
cesm2_simu_o_scalar = cesm2_simu_o_scalar.var_record_monthly_mean;
% FPI_vr
cesm2_simu_n_scalar = load([data_path, 'input_data/cesm2_simu/spinup_ss/',...
    cesm2_case_name, '_cesm2_ss_4da_', num2str(start_year), '_', num2str(end_year),...
    '_FPI_vr.mat'], 'var_record_monthly_mean');
cesm2_simu_n_scalar = cesm2_simu_n_scalar.var_record_monthly_mean;
% LITR1_INPUT_ACC_VECTOR
cesm2_simu_input_vector_litter1 = load([data_path, 'input_data/cesm2_simu/spinup_ss/',...
    cesm2_case_name, '_cesm2_ss_4da_', num2str(start_year), '_', num2str(end_year),...
    '_LITR1_INPUT_ACC_VECTOR.mat'], 'var_record_monthly_mean');
cesm2_simu_input_vector_litter1 = cesm2_simu_input_vector_litter1.var_record_monthly_mean;
% LITR2_INPUT_ACC_VECTOR
cesm2_simu_input_vector_litter2 = load([data_path, 'input_data/cesm2_simu/spinup_ss/',...
    cesm2_case_name, '_cesm2_ss_4da_', num2str(start_year), '_', num2str(end_year),...
    '_LITR2_INPUT_ACC_VECTOR.mat'], 'var_record_monthly_mean');
cesm2_simu_input_vector_litter2 = cesm2_simu_input_vector_litter2.var_record_monthly_mean;
% LITR3_INPUT_ACC_VECTOR
cesm2_simu_input_vector_litter3 = load([data_path, 'input_data/cesm2_simu/spinup_ss/',...
    cesm2_case_name, '_cesm2_ss_4da_', num2str(start_year), '_', num2str(end_year),...
    '_LITR3_INPUT_ACC_VECTOR.mat'], 'var_record_monthly_mean');
cesm2_simu_input_vector_litter3 = cesm2_simu_input_vector_litter3.var_record_monthly_mean;
% CWD_INPUT_ACC_VECTOR
cesm2_simu_input_vector_cwd = load([data_path, 'input_data/cesm2_simu/spinup_ss/',...
    cesm2_case_name, '_cesm2_ss_4da_', num2str(start_year), '_', num2str(end_year),...
    '_CWD_INPUT_ACC_VECTOR.mat'], 'var_record_monthly_mean');
cesm2_simu_input_vector_cwd = cesm2_simu_input_vector_cwd.var_record_monthly_mean;
% TOTSOMC
cesm2_simu_soc_stock = load([data_path, 'input_data/cesm2_simu/spinup_ss/',...
    cesm2_case_name, '_cesm2_ss_4da_', num2str(start_year), '_', num2str(end_year),...
    '_TOTSOMC.mat'], 'var_record_monthly_mean');
cesm2_simu_soc_stock = cesm2_simu_soc_stock.var_record_monthly_mean;


for ilayer = 1:soil_decom_num
    cesm2_simu_input_vector_litter1(:, :, ilayer, :) = cesm2_simu_input_vector_litter1(:, :, ilayer, :).*dz(ilayer);
    cesm2_simu_input_vector_litter2(:, :, ilayer, :) = cesm2_simu_input_vector_litter2(:, :, ilayer, :).*dz(ilayer);
    cesm2_simu_input_vector_litter3(:, :, ilayer, :) = cesm2_simu_input_vector_litter3(:, :, ilayer, :).*dz(ilayer);
    cesm2_simu_input_vector_cwd(:, :, ilayer, :) = cesm2_simu_input_vector_cwd(:, :, ilayer, :).*dz(ilayer);
end

cesm2_simu_input_sum_litter1 = reshape(sum(cesm2_simu_input_vector_litter1, 3), [384, 576, 12]);
cesm2_simu_input_sum_litter2 = reshape(sum(cesm2_simu_input_vector_litter2, 3), [384, 576, 12]);
cesm2_simu_input_sum_litter3 = reshape(sum(cesm2_simu_input_vector_litter3, 3), [384, 576, 12]);
cesm2_simu_input_sum_cwd = reshape(sum(cesm2_simu_input_vector_cwd, 3), [384, 576, 12]);


% clearvars cesm2_simu_input_vector_litter1 cesm2_simu_input_vector_litter2 cesm2_simu_input_vector_litter3 cesm2_simu_input_vector_cwd


%% Load NN Prediction Results
GlobalGrid = load([data_path, 'input_data/data4nn/world_grid_envinfo_present.mat']);
GlobalGrid = GlobalGrid.EnvInfo;

grid_para_result = csvread([data_path, 'output_data/neural_networking/grid_para_result_', model_name, '_', time_domain, '_', nn_exp_name, '.csv']);
valid_grid_loc = csvread([data_path, 'output_data/neural_networking/valid_grid_loc_', model_name, '_', time_domain, '_', nn_exp_name, '.csv']);

selected_grid_loc = randsample(length(valid_grid_loc), 5000, false);
grid_para_result = grid_para_result(selected_grid_loc, :);
%% load MCMC MLE results
% grid info
GlobalGrid = GlobalGrid(valid_grid_loc(selected_grid_loc), :);
GridInfo = GlobalGrid(:, [1:2, 12]);

ValidGrid = length(GridInfo(:, 1));
save([data_path, '/output_data/', model_name, '_recover/grid_info.mat', 'GridInfo']);
%%
% cesm2 resolution
cesm2_resolution_lat = 180/384;
cesm2_resolution_lon = 360/576;
LonGrid_Range = (-180 + cesm2_resolution_lon/2 : cesm2_resolution_lon : 180 - cesm2_resolution_lon/2)';
LatGrid_Range = (90 - cesm2_resolution_lat/2 : -cesm2_resolution_lat : -90 + cesm2_resolution_lat/2)';

parfor iGrid = 1:length(GridInfo)

    warning off
    disp(['Processing Grid: ', num2str(iGrid)]);

    Grid_Lon = GridInfo(iGrid,1);
    Grid_Lat = GridInfo(iGrid,2);

    lat_loc = find(abs(Grid_Lat - LatGrid_Range) == min(abs(Grid_Lat - LatGrid_Range)));
    lon_loc = find(abs(Grid_Lon - LonGrid_Range) == min(abs(Grid_Lon - LonGrid_Range)));

    if length(lon_loc) > 1
        lon_loc = lon_loc(1);
    end

    if length(lat_loc) > 1
        lat_loc = lat_loc(1);
    end

    % input vector
    input_vector_cwd = reshape(cesm2_simu_input_sum_cwd(lat_loc, lon_loc, :), [1, month_num]);
    input_vector_litter1 = reshape(cesm2_simu_input_sum_litter1(lat_loc, lon_loc, :), [1, month_num]);
    input_vector_litter2 = reshape(cesm2_simu_input_sum_litter2(lat_loc, lon_loc, :), [1, month_num]);
    input_vector_litter3 = reshape(cesm2_simu_input_sum_litter3(lat_loc, lon_loc, :), [1, month_num]);


    % no input information
    if isnan(input_vector_cwd(1)) == 1 || isnan(input_vector_litter1(1)) == 1 ...
            || isnan(input_vector_litter2(1)) == 1 || isnan(input_vector_litter3(1)) == 1
        disp(['No input info ', ' Grid ', num2str(iGrid)]);
        continue
    end

    % npp from CESM2 simulatoin
    npp_mean = reshape(cesm2_simu_npp(lat_loc, lon_loc, :), [month_num, 1]);
    % NPP info of studied profile from soc_modIS NPP mean (year 2001-2016)
    % npp_mean = env_info(iprofile, 15);

    % altmax current and last year
    altmax_current_profile = reshape(cesm2_simu_altmax(lat_loc, lon_loc, :), [month_num, 1]);
    altmax_lastyear_profile = reshape(cesm2_simu_altmax_last_year(lat_loc, lon_loc, :), [month_num, 1]);

    % nbedrock
    nbedrock = reshape(cesm2_simu_nbedrock(lat_loc, lon_loc, :), [month_num, 1]);

    % oxygen scalar
    xio = reshape(cesm2_simu_o_scalar(lat_loc, lon_loc, 1:soil_decom_num, :), [soil_decom_num, month_num, 1]);
    % nitrogen
    xin = reshape(cesm2_simu_n_scalar(lat_loc, lon_loc, 1:soil_decom_num, :), [soil_decom_num, month_num, 1]);

    % sand content
    sand = reshape(cesm2_simu_cellsand(lat_loc, lon_loc, 1:soil_decom_num, :), [soil_decom_num, month_num, 1]);
    sand_vector = sand;

    % soil temperature and water potential
    soil_temp_profile = reshape(cesm2_simu_soil_temperature(lat_loc, lon_loc, 1:soil_decom_num, :), [soil_decom_num, month_num, 1]);
    % soil_water_profile = reshape(cesm2_simu_soil_water_potnetial(lat_loc, lon_loc, 1:soil_decom_num, :), [soil_decom_num, month_num, 1]);
    % soil_temp_profile = reshape(cesm2_simu_t_scalar(lat_loc, lon_loc, 1:soil_decom_num, :), [soil_decom_num, month_num, 1]);
    soil_water_profile = reshape(cesm2_simu_w_scalar(lat_loc, lon_loc, 1:soil_decom_num, :), [soil_decom_num, month_num, 1]);

    % Parameter names and their initial values in MCMC
    para_name = {'bio', 'cryo', 'q10', 'w_scaling', ...
        'tau4cwd', 'tau4l1', 'tau4l2', 'tau4s1', 'tau4s2_death', 'tau4s2_enz', 'tau4s3', 'tau4s4', ...
        'mm_const_assim', 'mm_const_decom', ...
        'fcwdl2', ...
        'pl1s1', 'pl2s1', 'pl3s4', 'l1_cue', 'l2l3_cue', ...
        'mic_cue', 'pdeath2soc', ...
        'beta'};

    is_default = 0;
    Para = grid_para_result(iGrid, :);

    Para(Para > 1) = 1;
    Para(Para <= 0) = 10^(-4);

    [~, ~, pseudo_obs] = ...
        matrix_fun(Para, nbedrock, sand_vector, npp_mean, ...
        input_vector_cwd, input_vector_litter1, input_vector_litter2, input_vector_litter3, ...
        altmax_current_profile, altmax_lastyear_profile, soil_temp_profile, soil_water_profile, xio, xin);
    pseudo_obs = pseudo_obs(:, 5);

    if isempty(find(pseudo_obs < 0, 1)) == 0 || isnan(sum(pseudo_obs)) == 1
        continue
    end

    %% MCMC
    npara = length(para_name);
    para_min = zeros(npara, 1);
    para_max = ones(npara, 1);
    % info of the node depth of profile, and change unit from cm to m
    wosis_layer_depth = zsoi(1:20);
    % observed C info (gC/m3)
    wosis_layer_obs = pseudo_obs;
    % specify the number of layers in studied profiles
    num_layers = length(wosis_layer_obs);

    % layer_weight([1, end], 1) = 10;
    [~, soc_stock_rank] = sort(abs(wosis_layer_obs - mean(wosis_layer_obs)), 'descend');
    [~, soc_stock_rank] = sort(soc_stock_rank);
    [~, soc_depth_rank] = sort(abs(wosis_layer_depth - mean(wosis_layer_depth)), 'descend');
    [~, soc_depth_rank] = sort(soc_depth_rank);
    soc_rank = (soc_stock_rank + soc_depth_rank)/2;

    % plot(wosis_layer_depth, wosis_layer_obs)

    layer_weight = exp(-(soc_rank - 1))*10;
    layer_weight(layer_weight < 1) = 1;

    % profiles
    optimize_profile_soc = -999*ones(20, 1);


    %% MCMC core part
    % try
        is_second_time = 0;
        % Predefine the size of some coefficients
        nsimu1 = 2000;
        nsimu2 = 5000;
        nparallel = 3;

        parameters_rec1 = nan(nparallel, npara, nsimu1);
        parameters_rec2 = nan(nparallel, npara, nsimu2);
        parameters_keep1 = nan(nparallel, npara, nsimu1);
        parameters_keep2 = nan(nparallel, npara, nsimu2);

        cost_record1 = nan(nsimu1, nparallel);
        cost_record2 = nan(nsimu2, nparallel);
        cost_keep1 = nan(nsimu1, nparallel);
        cost_keep2 = nan(nsimu2, nparallel);

        upgrade1 = nan(nparallel, 1);
        upgrade2 = nan(nparallel, 1);

        soc_mod_trace = nan(nparallel, nsimu2, length(wosis_layer_depth));
        %calculating prior cost function value
        cost_old1 = cost_fun(layer_weight, optimize_profile_soc, wosis_layer_obs);

        sd_controling_factor = 2.4; %2.4; % default to be 2.4^2

        %% Part1: MCMC run to calculate parameter covariances
        % iparallel = 1;
        parallel_loc = 1;
        counter2 = 0;
        for iparallel = 1 : 6
            parameters_rec1(parallel_loc, :, :) = nan;
            parameters_rec2(parallel_loc, :, :) = nan;
            parameters_keep1(parallel_loc, :, :) = nan;
            parameters_keep2(parallel_loc, :, :) = nan;

            cost_record1(:, parallel_loc) = nan;
            cost_record2(:, parallel_loc) = nan;
            cost_keep1(:, parallel_loc) = nan;
            cost_keep2(:, parallel_loc) = nan;

            soc_mod_trace(parallel_loc, :, :) = nan;

            cost_old = cost_old1;
            cost_control = 1;
            allow = 5*ones(npara, 1);
            % allow(strcmp(para_name, 'tau4p')) = 90;
            % give original parameter values to par_old
            par_old = 0.3*ones(npara, 1);
            % interval width
            diff = para_max - para_min;
            upgrade1(parallel_loc, 1) = 0;
            simu = 1;
            counter1 = 0;
            ifbreak1 = 0;
            while simu <= nsimu1
                while (true)
                    % generate new parameter values in the interval of
                    % [-0.5, 0.5]*diff/allow + par_old
                    par_new = par_old+(rand(npara,1)-0.5).*diff./allow;
                    % all the new valued parameters should be in the preassumed
                    % intervals
                    if isempty(find((par_new - para_min) <=0, 1)) == 1 && ...
                            isempty(find((par_new - para_max) >=0, 1)) == 1
                        break;
                    end
                end


                % soc_model simulation

                % clear warning info
                lastwarn('');
                [~, ~, soc_mod] = ...
                    matrix_fun(par_new, nbedrock, sand_vector, npp_mean, ...
                    input_vector_cwd, input_vector_litter1, input_vector_litter2, input_vector_litter3, ...
                    altmax_current_profile, altmax_lastyear_profile, soil_temp_profile, soil_water_profile, xio, xin);

                % original soc_modelled data
                [~, msgid] = lastwarn;
                if strcmp(msgid, 'MATLAB:singularMatrix') || strcmp(msgid, 'MATLAB:nearlySingularMatrix') ...
                        || isnan(sum(soc_mod(:, 5))) == 1
                    cost_new = cost_old + 100;
                else
                    % write in the newly soc_modeled data
                    optimize_profile_soc = soc_mod(:, 5);
                    % new cost function
                    cost_new = cost_fun(layer_weight, optimize_profile_soc, wosis_layer_obs);

                end

                delta_cost = cost_new-cost_old;
                % to decide whether to accept new parameter values
                if min(1, exp(-delta_cost*log(simu + 2))) > rand
                    % update the number of simulation
                    upgrade1(parallel_loc, 1) = upgrade1(parallel_loc, 1) + 1;

                    % disp(['upgrade1_parallel_', num2str(parallel_loc), ':', num2str(upgrade1(parallel_loc, 1))]);

                    % record accepted parameter values
                    parameters_keep1(parallel_loc, :, upgrade1(parallel_loc, 1)) = par_new;
                    % record accepted values of cost function
                    cost_keep1(upgrade1(parallel_loc, 1), parallel_loc)=cost_new;
                    % update the value of parameter values
                    par_old = par_new;
                    % update old cost function value
                    cost_old = cost_new;
                end
                % give parameter values to parameters_rec for covarance matrix,
                % rec: record
                parameters_rec1(parallel_loc, :, simu)=par_old;
                cost_record1(simu, parallel_loc)=cost_old;
                simu = simu + 1;

                if simu == nsimu1 && is_second_time == 0
                    if upgrade1(parallel_loc, 1) < 50
                        counter1 = counter1 + 1;
                        simu = 1;
                        upgrade1(parallel_loc, 1) = 0;
                        parameters_keep1(parallel_loc, :, :) = nan;
                        parameters_rec1(parallel_loc, :, :) = nan;
                        cost_record1(:, parallel_loc) = nan;
                        cost_keep1(:, parallel_loc) = nan;
                        allow = allow*1.1;
                    end
                end
                if counter1 == 3
                    ifbreak1 = 1;
                    is_second_time = 1;
                    break
                end
            end
            if ifbreak1 == 1
                continue
            end

            %% Part 2: MCMC run with updating covariances
            % sd_controling_factor = 2; % default to be 2.4^2
            sd = sd_controling_factor/npara;
            epsilon = 0;
            covars=sd*cov(reshape(parameters_rec1(parallel_loc, :, 1:nsimu1), [npara, length(1:nsimu1)])') + sd*epsilon*diag(ones(npara, 1));
            % covars=sd*cov(reshape(parameters_keep1(iparallel, :, 1:upgrade1(iparallel, 1)), [npara, upgrade1(iparallel, 1)])') + sd*epsilon*diag(ones(npara, 1));
            % update the value of upgrade and simu for new task
            upgrade2(parallel_loc, 1) = 0;
            % the very first cost functino value in simulation
            cost_old = cost_old1;
            simu = 1;
            ifbreak2 = 0;
            ifbreak3 = 0;
            while simu <= nsimu2
                monitor_time_2 = tic;
                while (true)
                    % mointor the consumed time at the initial stage (if too long i.e. 10min, then back to the first proposal step)
                    if upgrade2(parallel_loc, 1) == 0
                        consume_time_2 = toc(monitor_time_2);
                        if consume_time_2 > 1*60
                            ifbreak2 = 1;
                        end
                    end
                    if ifbreak2 == 1
                        break
                    end
                    try
                        par_new = mvnrnd(par_old, covars)';
                    catch
                        ifbreak3 = 1;
                        break
                    end
                    if isempty(find((par_new - para_min) <=0, 1)) == 1 && ...
                            isempty(find((par_new - para_max) >=0, 1)) == 1
                        break;
                    end
                end
                if ifbreak2 == 1 || ifbreak3 == 1
                    break
                end

                lastwarn('');
                % soc_model simulation
                [~, ~, soc_mod] = ...
                    matrix_fun(par_new, nbedrock, sand_vector, npp_mean, ...
                    input_vector_cwd, input_vector_litter1, input_vector_litter2, input_vector_litter3, ...
                    altmax_current_profile, altmax_lastyear_profile, soil_temp_profile, soil_water_profile, xio, xin);

                % original soc_modelled data
                [~, msgid] = lastwarn;
                if strcmp(msgid, 'MATLAB:singularMatrix') || strcmp(msgid, 'MATLAB:nearlySingularMatrix') ...
                        || isnan(sum(soc_mod(:, 5))) == 1
                    cost_new = cost_old + 100;
                else
                    % write in the newly soc_modeled data
                    optimize_profile_soc = soc_mod(:, 5);
                    % new cost function
                    cost_new = cost_fun(layer_weight, optimize_profile_soc, wosis_layer_obs);
                end

                delta_cost = cost_new-cost_old;

                % to determie whether to accept the newly soc_modelled data
                if min(1, exp(-delta_cost*log(simu + 2))) > rand
                    upgrade2(parallel_loc, 1) = upgrade2(parallel_loc, 1) + 1;

                    % disp(['upgrade2_parallel_', num2str(parallel_loc), ':', num2str(upgrade2(parallel_loc, 1))]);

                    % record the accepted parameter value
                    parameters_keep2(parallel_loc, :,upgrade2(parallel_loc, 1)) = par_new;
                    cost_keep2(upgrade2(parallel_loc, 1), parallel_loc) = cost_new;
                    par_old = par_new;
                    cost_old = cost_new;
                    coef = upgrade2(parallel_loc, 1)/simu;
                    soc_mod_trace(parallel_loc, upgrade2(parallel_loc, 1), :) = optimize_profile_soc;
                end
                parameters_rec2(parallel_loc, :, simu) = par_old;
                cost_record2(simu, parallel_loc) = cost_old;
                if simu > 4000
                    % covars=sd*(cov(reshape(parameters_keep2(iparallel, :, 1:upgrade2(iparallel, 1)), [npara, upgrade2(iparallel, 1)])')) + sd*epsilon*diag(ones(npara, 1));
                    covars=sd*(cov(reshape(parameters_rec2(parallel_loc, :, 1:simu), [npara, length(1:simu)])')) + sd*epsilon*diag(ones(npara, 1));
                end
                simu = simu + 1;
                % to test if the acceptance rate is in a resonable level
                if simu == nsimu2
                    if coef < 0.1 || coef > 0.5
                        ifbreak2 = 1;
                        counter2 = counter2 + 1;

                        if coef < 0.1
                            sd_controling_factor = sd_controling_factor*(1 - 0.3);
                        end

                        if coef > 0.5
                            sd_controling_factor = sd_controling_factor*2;
                        end

                        break
                    end
                end
            end
            if ifbreak2 == 0 && ifbreak3 == 0
                % iparallel = iparallel + 1;
                parallel_loc = parallel_loc + 1;
            end
            if counter2 > 5
                error('Error')
            end
        end
        % disp('MCMC has been finished');


        %% Figures (only used in single profile test
        % % autocorrelation plot
        % plot_simu = 1;
        % for ipara = 1:length(para_name)
        %     autocorr_log = 200;
        %     subplot(5, 5, ipara)
        %     autocorr(reshape(parameters_keep2(plot_simu, ipara, round(upgrade2(plot_simu)/2):upgrade2(plot_simu)), [length(round(upgrade2(plot_simu)/2):upgrade2(plot_simu)), 1]), 'NumLags', autocorr_log);
        %     title(para_name{ipara})
        % end
        %
        % close
        %
        % % obs v.s.modelled soc
        % scatter(wosis_layer_depth, wosis_layer_obs, 100, '*')
        % hold on
        % plot_mod_mean = reshape(mean(soc_mod_trace(plot_simu, round(upgrade2(plot_simu)/2):upgrade2(plot_simu), :), 2), [length(wosis_layer_depth), 1]);
        % plot_mod_std =  reshape(std(soc_mod_trace(plot_simu, round(upgrade2(plot_simu)/2):upgrade2(plot_simu), :), 0, 2), [length(wosis_layer_depth), 1]);
        % errorbar(wosis_layer_depth, plot_mod_mean, plot_mod_std, '*', 'MarkerSize', 10)
        %
        % close
        %
        %
        % % histogram plot
        % for ipara = 1:length(para_name)
        %     hist_bin_num = 30;
        %     subplot(5, 5, ipara)
        %     histogram(parameters_keep2(1, ipara, round(upgrade2(1)/2):upgrade2(1)), 'BinWidth', 0.03 );
        %     hold on
        %     % histogram(parameters_keep2(2, ipara, round(upgrade2(2)/2):upgrade2(2)), 'BinWidth', 0.03 );
        %     % hold on
        %     % histogram(parameters_keep2(3, ipara, round(upgrade2(3)/2):upgrade2(3)), 'BinWidth', 0.03 );
        %     % hold on
        %     xline(Para(ipara), 'LineWidth', 3, 'Color', 'red')
        %     title(para_name{ipara})
        %     xlim([0 1])
        % end
        %
        % close

        %% Test of convergence
        % set the first half number as burn-in
        valid_num = nan(nparallel, 1);
        soc_mod_parallel = nan(length(wosis_layer_obs), nparallel);
        soc_std_parallel = nan(length(wosis_layer_obs), nparallel);
        % col1: correlation, col2: RMSE
        stat_parallel = nan(nparallel, 2);
        for iparallel = 1 : nparallel
            parallel_loc = iparallel;
            para_keep_middle = reshape(parameters_keep2(parallel_loc, 1, :), [nsimu2, 1]);
            % find valid soc_modeled number
            valid_num(parallel_loc, 1) = length(para_keep_middle(~isnan(para_keep_middle), 1));

            soc_mod_parallel(:, parallel_loc) = mean(reshape(soc_mod_trace(parallel_loc, round(upgrade2(parallel_loc, 1)/2) : upgrade2(parallel_loc, 1), :),...
                [length(round(upgrade2(parallel_loc, 1)/2) : upgrade2(parallel_loc, 1)), length(wosis_layer_obs)]));

            soc_std_parallel(:, parallel_loc) = std(reshape(soc_mod_trace(parallel_loc, round(upgrade2(parallel_loc, 1)/2) : upgrade2(parallel_loc, 1), :),...
                [length(round(upgrade2(parallel_loc, 1)/2) : upgrade2(parallel_loc, 1)), length(wosis_layer_obs)]));

            SStot = sum((soc_mod_parallel(:, parallel_loc) - mean(soc_mod_parallel(:, parallel_loc))).^2);
            SSres = sum((soc_mod_parallel(:, parallel_loc) - wosis_layer_obs).^2);
            % R squared
            stat_parallel(parallel_loc, 1) = 1 - SSres/SStot;
            % RMSE
            stat_parallel(parallel_loc, 2) = sqrt(sum((wosis_layer_obs - soc_mod_parallel(:, parallel_loc)).^2)/num_layers);
        end

        % set the parallel has the lowest RMSE as the optimal
        opt_parallel = find(stat_parallel(:, 2) == min(stat_parallel(:, 2)));
        if length(opt_parallel) ~= 1
            opt_parallel = find(stat_parallel(opt_parallel, 1) == para_max(stat_parallel(opt_parallel, 1)));
        end

        para_after_burnin = ...
            reshape(parameters_keep2(opt_parallel, :, round(valid_num(opt_parallel, 1)/2)+1:valid_num(opt_parallel, 1)),...
            [npara, length(round(valid_num(opt_parallel, 1)/2)+1:valid_num(opt_parallel, 1))]);

        para_summary_mean = mean(para_after_burnin, 2, 'omitnan');
        para_summary_std = std(para_after_burnin, 0, 2, 'omitnan');


        soc_mod_opt = soc_mod_parallel(:, opt_parallel);
        soc_std_opt = soc_std_parallel(:, opt_parallel);


        %% Para MLE
        para_es = nan(npara, 1);
        para_mle = nan(npara, 3);
        for i = 1 : npara
            % para_mle = mle(Parameters_keep(i, upgrade-2000:(upgrade-1)),'distribution','gev');
            % [para_mle_middle, ~] = mle(ParaAfterBurnin(i, :), 'distribution', 'gev');
            para_mle_middle = fitdist(para_after_burnin(i, :)', 'gev');
            para_mle(i, :) = [para_mle_middle.k, para_mle_middle.sigma, para_mle_middle.mu];
            % the mode of gev distribution
            if para_mle_middle.k == 0
                para_es(i) = para_mle_middle.mu;
            else
                para_es(i) = para_mle_middle.mu + ...
                    para_mle_middle.sigma*((1+para_mle_middle.k)^(-para_mle_middle.k)-1)/para_mle_middle.k;
            end
        end

        %% Gelman-Rubin test
        GR = nan(npara, 1);
        between_var = nan(npara, 1);
        within_var = nan(npara, 1);
        for ipara = 1 : npara
            para_mean = nan(nparallel, 1);
            para_sum = nan(nparallel, 1);
            for iparallel = 1 : nparallel
                parallel_loc = iparallel;
                para_keep_middle = reshape(parameters_keep2(parallel_loc, 1, :), [nsimu2, 1]);
                valid_num = length(para_keep_middle(~isnan(para_keep_middle), 1));
                para_mean(parallel_loc, 1) = mean(parameters_keep2(parallel_loc, ipara, round(valid_num/2)+1:valid_num), 'omitnan');
                para_sum(parallel_loc, 1) = sum((parameters_keep2(parallel_loc, ipara, round(valid_num/2)+1:valid_num) - para_mean(parallel_loc, 1)).^2);
            end
            % calculate between variance and within variance
            between_var = ...
                valid_num/(nparallel-1)*nansum((para_mean - mean(para_mean)).^2);
            within_var = 1/nparallel/(valid_num-1)*sum(para_sum);
            GR(ipara) = sqrt(((within_var*(valid_num-1))/valid_num + between_var/valid_num)/within_var)';
        end

        %% Saving
        MatFileName = 'none';
        fun_save_da_recover(model_name, iGrid,...
            Para, wosis_layer_obs, soc_mod_opt, soc_std_opt, ...
            para_summary_mean, para_summary_std,...
            para_es,...
            GR, parameters_keep2, stat_parallel,...
            upgrade2);

        disp([datestr(now), ' Profile ', num2str(iGrid), ' has been finished']);

    % catch
    %     disp(['MCMC_Error_in_Grid_', num2str(iGrid)]);
    % end
end

disp([datestr(now), ' Program finished']);
