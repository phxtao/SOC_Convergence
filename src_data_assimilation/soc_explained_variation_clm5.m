clear;
clc;
start_id = 1;
end_id = 1;
%%
% parallel setting
% delete(gcp('nocreate'));
% parpool(10);

cesm2_case_name = 'sasu_f05_g16_checked_step4';
start_year = 661;
end_year = 680;

model_name = 'cesm2_clm5_cen_vr_v2'; % 'cesm2_clm5_cen_vr_v2'; % 'cesm2_clm5_cen_vr_v2';

time_domain = 'whole_time'; % 'whole_time', 'before_1985', 'after_1985', 'random_half', 'random_half_2'


for ibootstrap = 1:10
    disp([datestr(now), ' Processing ibootstrap ', num2str(ibootstrap)]);
    exp_name = ['exp_pc_cesm2_23_cross_valid_0_', num2str(ibootstrap)];
    % exp_name = ['exp_pc_cesm2_23_bootstrap_', num2str(ibootstrap)];
    
    
    %% paths
    % mac
    data_path = '/Users/ft254/DATAHUB/ENSEMBLE/';
    cd(['/Users/ft254/Github/ENSEMBLE/SRC_DA/src_', model_name, '/']);
    
    % server
    % data_path = '/GFPS8p/cess11/taof/datahub/ensemble/';
    % cd(['/GFPS8p/cess11/taof/ensemble/src_da/src_', model_name, '/']);
    
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
    
    
    cesm2_simu_input_sum_litter1 = reshape(sum(cesm2_simu_input_vector_litter1, 4), [384, 576, soil_decom_num]);
    cesm2_simu_input_sum_litter2 = reshape(sum(cesm2_simu_input_vector_litter2, 4), [384, 576, soil_decom_num]);
    cesm2_simu_input_sum_litter3 = reshape(sum(cesm2_simu_input_vector_litter3, 4), [384, 576, soil_decom_num]);
    cesm2_simu_input_sum_cwd = reshape(sum(cesm2_simu_input_vector_cwd, 4), [384, 576, soil_decom_num]);
    
    % clearvars cesm2_simu_input_vector_litter1 cesm2_simu_input_vector_litter2 cesm2_simu_input_vector_litter3 cesm2_simu_input_vector_cwd

    %% load results from neural networking
    nn_predict = csvread([data_path, 'output_data/neural_networking/nn_para_result_full_set_', model_name, '_', time_domain, '_', exp_name, '.csv']);
    nn_site_loc = csvread([data_path, 'output_data/neural_networking/nn_site_loc_full_set_', model_name, '_', time_domain, '_', exp_name, '.csv']);
    
    %% Individual Simulation
    indi_map_project_mod = nan(length(nn_site_loc), 150, 10);
    indi_map_project_obs = nan(length(nn_site_loc), 150, 10);
    indi_map_project_depth = nan(length(nn_site_loc), 150, 10);
    
    site_by_site_mod = nan(length(nn_site_loc), soil_decom_num);
    
    for iprofile_hat = 1:length(nn_site_loc)
        
        iprofile = nn_site_loc(iprofile_hat);
        warning('off')
        
        disp([datestr(now), ' Processing profile ', num2str(iprofile_hat)]);
        
        profile_id = wosis_profile_info(iprofile, 1);
        % find currently using profile
        loc_profile = find(wosis_soc_info(:, 1) == profile_id);
        % info of the node depth of profile, and change unit from cm to m
        wosis_layer_depth = wosis_soc_info(loc_profile, 5)/100;
        % observed C info (gC/m3)
        wosis_layer_obs = wosis_soc_info(loc_profile, 7);
        % specify the number of layers in studied profiles
        num_layers = length(wosis_layer_obs);
        
        % find the lon and lat info of soil profile
        lon_profile = wosis_profile_info(iprofile, 4);
        lat_profile = wosis_profile_info(iprofile, 5);
        lat_loc = find(abs(lat_profile - lat_grid) == min(abs(lat_profile - lat_grid)));
        lon_loc = find(abs(lon_profile - lon_grid) == min(abs(lon_profile - lon_grid)));
        
        if length(lon_loc) > 1
            lon_loc = lon_loc(1);
        end
        
        if length(lat_loc) > 1
            lat_loc = lat_loc(1);
        end
        
        % input vector
        % input_vector_cwd = reshape(cesm2_simu_input_vector_cwd(lat_loc, lon_loc, :, :), [soil_decom_num, month_num]);
        % input_vector_litter1 = reshape(cesm2_simu_input_vector_litter1(lat_loc, lon_loc, :, :), [soil_decom_num, month_num]);
        % input_vector_litter2 = reshape(cesm2_simu_input_vector_litter2(lat_loc, lon_loc, :, :), [soil_decom_num, month_num]);
        % input_vector_litter3 = reshape(cesm2_simu_input_vector_litter3(lat_loc, lon_loc, :, :), [soil_decom_num, month_num]);
        
        % input_vector_cwd = reshape(cesm2_simu_input_sum_cwd(lat_loc, lon_loc, :), [1, month_num]);
        % input_vector_litter1 = reshape(cesm2_simu_input_sum_litter1(lat_loc, lon_loc, :), [1, month_num]);
        % input_vector_litter2 = reshape(cesm2_simu_input_sum_litter2(lat_loc, lon_loc, :), [1, month_num]);
        % input_vector_litter3 = reshape(cesm2_simu_input_sum_litter3(lat_loc, lon_loc, :), [1, month_num]);
        
        input_vector_cwd = reshape(cesm2_simu_input_sum_cwd(lat_loc, lon_loc, :), [soil_decom_num, 1]);
        input_vector_litter1 = reshape(cesm2_simu_input_sum_litter1(lat_loc, lon_loc, :), [soil_decom_num, 1]);
        input_vector_litter2 = reshape(cesm2_simu_input_sum_litter2(lat_loc, lon_loc, :), [soil_decom_num, 1]);
        input_vector_litter3 = reshape(cesm2_simu_input_sum_litter3(lat_loc, lon_loc, :), [soil_decom_num, 1]);
        
        % no input information
        if isnan(input_vector_cwd(1)) == 1 || isnan(input_vector_litter1(1)) == 1 ...
                || isnan(input_vector_litter2(1)) == 1 || isnan(input_vector_litter3(1)) == 1
            disp(['No input info ', ' Profile ', num2str(iprofile)]);
            continue
        end
        
        % npp from CESM2 simulatoin
        npp_mean = reshape(cesm2_simu_npp(lat_loc, lon_loc, :), [month_num, 1]);
        % NPP info of studied profile from soc_modIS NPP mean (year 2001-2016)
        % modis_npp_mean = env_info(iprofile, 15);
        
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
        
        %% data cleansing
        % deal with zero value in mcmc
        % if zero, cost function will always be inf
        if length(find(wosis_layer_obs == 0)) == length(wosis_layer_obs)
            disp(['All zero value ', ' Profile ', num2str(iprofile)]);
            continue
        elseif isempty(find(wosis_layer_obs == 0, 1)) == 0
            wosis_layer_depth = wosis_layer_depth(wosis_layer_obs > 0);
            wosis_layer_obs = wosis_layer_obs(wosis_layer_obs > 0);
        end
        
        % deal with nan values in obs or depth
        layer_valid_loc = find(isnan(wosis_layer_depth) == 0 & isnan(wosis_layer_obs) == 0);
        
        if isempty(layer_valid_loc) == 0
            wosis_layer_depth = wosis_layer_depth(layer_valid_loc);
            wosis_layer_obs = wosis_layer_obs(layer_valid_loc);
        else
            disp(['No valid obsrvations in ', ' Profile ', num2str(iprofile)]);
            continue
        end
        
        % eliminate profiles having single layer
        if length(wosis_layer_obs) == 1
            disp(['Single layer ', ' Profile ', num2str(iprofile)]);
            continue
        end
        
        %% Parameterization
        para_name = {'diffus'; 'cryo';...
            'q10';...
            'efolding';...
            'taucwd'; 'taul1'; 'taul2';...
            'tau4s1'; 'tau4s2'; 'tau4s3';...
            'fl1s1'; 'fl2s1'; 'fl3s2'; 'fs1s2'; 'fs1s3'; 'fs2s1'; 'fs2s3'; 'fs3s1'; 'fcwdl2';...
            'w-scaling'; 'beta'};
                
        %% Fully Varying (NN)
        % if using default parameterization
        is_default = 0;
        % parameters values
        para = nn_predict(iprofile_hat, :);
        
        % clear warning info
        lastwarn('');
        [~, ~, ~, soc_mod, ~, ~, ~, ~, ~] ...
            = fun_var_decom(para, is_default, nbedrock, sand_vector, npp_mean, ...
            input_vector_cwd, input_vector_litter1, input_vector_litter2, input_vector_litter3, ...
            altmax_current_profile, altmax_lastyear_profile, soil_temp_profile, soil_water_profile, xio, xin);
        
        [~, msgid] = lastwarn;
        if strcmp(msgid, 'MATLAB:singularMatrix') || strcmp(msgid, 'MATLAB:nearlySingularMatrix')
            disp(['Assilimation_failed_in_Profile_', num2str(iprofile)]);
            % msgid = [];
            continue
        end
        
        if isnan(soc_mod(1)) == 1
            disp(['nan value ', ' Profile ', num2str(iprofile)]);
            continue
        end
        
        % profiles
        if length(wosis_layer_obs) > 1 && isnan(soc_mod(1)) == 0
            optimize_profile_soc = interp1(zsoi(1:soil_decom_num, 1), soc_mod, wosis_layer_depth, 'pchip');
            % scatter(current_obsC(:,1), current_obsC(:,2));
            site_by_site_mod(iprofile_hat, :) = soc_mod;
            indi_map_project_mod(iprofile_hat, 1:length(wosis_layer_obs), ibootstrap) = optimize_profile_soc;
            indi_map_project_obs(iprofile_hat, 1:length(wosis_layer_obs), ibootstrap) = wosis_layer_obs;
            indi_map_project_depth(iprofile_hat, 1:length(wosis_layer_obs), ibootstrap) = wosis_layer_depth;
        end
    end
end

%% save simu results
proda_explained_variation.soc_proda_mean = mean(indi_map_project_mod, 3, 'omitnan');
proda_explained_variation.soc_obs = mean(indi_map_project_obs, 3, 'omitnan');
proda_explained_variation.soc_obs_depth = mean(indi_map_project_depth, 3, 'omitnan');

save([data_path, 'output_data/converged_soc/proda_explained_variation_', model_name, '.mat'], 'proda_explained_variation')

disp('Programme Finished');
