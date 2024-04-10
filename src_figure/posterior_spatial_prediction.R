## Packages
library(R.matlab)
library(maps)
library(ggplot2)
library(R.matlab)
library(cowplot)
library(viridis)
library(cowplot)
library(scales)

# dev.off()
##
rm(list = ls())

setwd('/Users/ft254/Google_Drive/R')



#############################################################################
# Data Path
#############################################################################
data_dir_output = '/Users/ft254/DATAHUB/ENSEMBLE/OUTPUT_DATA/'
data_dir_input = '/Users/ft254/DATAHUB/ENSEMBLE/INPUT_DATA/'

#################################################################################
# profile info
#################################################################################
grid_var_names = c('Lon', 'Lat', 'Date', 
                   'Rmean', 'Rmax', 'Rmin', 
                   'ESA_Land_Cover', 
                   'ET',
                   'IGBP', 'Climate', 'Soil_Type', 'NPPmean', 'NPPmax', 'NPPmin',
                   'Veg_Cover', 
                   'Annual Mean Temperature', 'Mean Diurnal Range', 'Isothermality', 'Temperature Seasonality', 'Max Temperature of Warmest Month', 'Min Temperature of Coldest Month', 'Temperature Annual Range', 'Mean Temperature of Wettest Quarter', 'Mean Temperature of Driest Quarter', 'Mean Temperature of Warmest Quarter', 'Mean Temperature of Coldest Quarter', 'Annual Precipitation', 'Precipitation of Wettest Month', 'Precipitation of Driest Month', 'Precipitation Seasonality', 'Precipitation of Wettest Quarter', 'Precipitation of Driest Quarter', 'Precipitation of Warmest Quarter', 'Precipitation of Coldest Quarter', 
                   'Abs_Depth_to_Bedrock',
                   'Bulk_Density_0cm', 'Bulk_Density_30cm', 'Bulk_Density_100cm',
                   'CEC_0cm', 'CEC_30cm', 'CEC_100cm',
                   'Clay_Content_0cm', 'Clay_Content_30cm', 'Clay_Content_100cm',
                   'Coarse_Fragments_v_0cm', 'Coarse_Fragments_v_30cm', 'Coarse_Fragments_v_100cm', 
                   'Depth_Bedrock_R', 
                   'Garde_Acid', 
                   'Occurrence_R_Horizon', 
                   'pH_Water_0cm', 'pH_Water_30cm', 'pH_Water_100cm', 
                   'Sand_Content_0cm', 'Sand_Content_30cm', 'Sand_Content_100cm',
                   'Silt_Content_0cm', 'Silt_Content_30cm', 'Silt_Content_100cm', 
                   'SWC_v_Wilting_Point_0cm', 'SWC_v_Wilting_Point_30cm', 'SWC_v_Wilting_Point_100cm', 
                   'Texture_USDA_0cm', 'Texture_USDA_30cm', 'Texture_USDA_100cm', 
                   'USDA_Suborder', 
                   'WRB_Subgroup', 
                   'Drought',
                   'Elevation',
                   'Max_Depth', 
                   'Koppen_Climate_2018', 
                   'cesm2_npp', 'cesm2_npp_std',
                   'cesm2_gpp', 'cesm2_gpp_std',
                   'cesm2_vegc',
                   'nbedrock')

env_info_ss = readMat(paste(data_dir_input, 'wosis_2019_snap_shot/wosis_2019_snapshot_hugelius_mishra_env_info.mat', sep = ''))
env_info_ss = env_info_ss$EnvInfo[ , 4:80]
colnames(env_info_ss) = grid_var_names

# load climate info
global_climate = env_info_ss[ , 'Koppen_Climate_2018']
var_climate = array(NA, dim = c(length(global_climate), 1))
# koppen climatte class 2
var_climate = global_climate
var_climate[which(global_climate == 1)] = 101 # Af
var_climate[which(global_climate == 2)] = 102 # Am
var_climate[which(global_climate == 3)] = 103 # Aw
var_climate[which(global_climate >= 4 & global_climate <= 5)] = 104 # BwX
var_climate[which(global_climate >= 6 & global_climate <= 7)] = 105 # BsW
var_climate[which(global_climate >= 8 & global_climate <= 10)] = 106 # CsX
var_climate[which(global_climate >= 11 & global_climate <= 13)] = 107 # CwX
var_climate[which(global_climate >= 14 & global_climate <= 16)] = 108 # CfX
var_climate[which(global_climate >= 17 & global_climate <= 20)] = 109 # DsX
var_climate[which(global_climate >= 21 & global_climate <= 24)] = 110 # DwX
var_climate[which(global_climate >= 25 & global_climate <= 28)] = 111 # DfX
var_climate[which(global_climate >= 29 & global_climate <= 30)] = 112 # E

# load soil order
global_soilorder = env_info_ss[ , 'USDA_Suborder']
var_soilorder = array(NA, dim = c(length(global_soilorder), 1))

var_soilorder[] = 113 # others
var_soilorder[which(global_soilorder >= 5 & global_soilorder <= 7)] = 101 # Gelisols
var_soilorder[which(global_soilorder >= 10 & global_soilorder <= 13)] = 102 # Histosols
var_soilorder[which(global_soilorder >= 15 & global_soilorder <= 19)] = 103 # Spodosols
var_soilorder[which(global_soilorder >= 20 & global_soilorder <= 27)] = 104 # Andisols
var_soilorder[which(global_soilorder >= 30 & global_soilorder <= 34)] = 105 # Oxisols
var_soilorder[which(global_soilorder >= 40 & global_soilorder <= 45)] = 106 # Vertisols
var_soilorder[which(global_soilorder >= 50 & global_soilorder <= 56)] = 107 # Aridisols
var_soilorder[which(global_soilorder >= 60 & global_soilorder <= 64)] = 108 # Ultisols
var_soilorder[which(global_soilorder >= 69 & global_soilorder <= 77)] = 109 # Mollisols
var_soilorder[which(global_soilorder >= 80 & global_soilorder <= 84)] = 110 # Alfisols
var_soilorder[which(global_soilorder >= 85 & global_soilorder <= 86)] = 111 # Inceptisols
var_soilorder[which(global_soilorder >= 89 & global_soilorder <= 94)] = 111 # Inceptisols
var_soilorder[which(global_soilorder >= 95 & global_soilorder <= 99)] = 112 # Entisols

# load ESA land cover
global_landcover = env_info_ss[ , 'ESA_Land_Cover']
var_landcover = array(NA, dim = c(length(global_landcover), 1))

var_landcover[which(global_landcover >= 1 & global_landcover <= 2)] = 101 # agriculture
var_landcover[which(global_landcover >= 3 & global_landcover <= 4)] = 102 # mosaic agriculture
var_landcover[which(global_landcover >= 5 & global_landcover <= 6)] = 103 # broadleaved forest 
var_landcover[which(global_landcover >= 7 & global_landcover <= 8)] = 104 # needleleaved forest 
var_landcover[which(global_landcover >= 9 & global_landcover <= 9)] = 105 # mixed forest 
var_landcover[which(global_landcover >= 10 & global_landcover <= 11)] = 106 # Mosaic tree and shrub
var_landcover[which(global_landcover >= 12 & global_landcover <= 12)] = 107 # shrub
var_landcover[which(global_landcover >= 13 & global_landcover <= 13)] = 108 # grassland
var_landcover[which(global_landcover >= 14 & global_landcover <= 15)] = 109 # Lichen & mosses, sparse vegatation
var_landcover[which(global_landcover >= 16 & global_landcover <= 18)] = 110 # wetland
var_landcover[which(global_landcover >=19)] = 111 # urban and other


#################################################################################
# valid loc
#################################################################################
# COMPAS model
model_name = 'cesm2_clm5_mic_vr_v22'
nn_exp_name = 'exp_pc_cesm2_23'
time_domain = 'whole_time'

da_summary = readMat(paste(data_dir_output, 'mcmc_summary_', model_name, '/', model_name, '_da_summary_mic.mat', sep = ''))
da_summary = da_summary$da.summary.mic
colnames(da_summary) = c('profile_num', 'mic_cue', 
                         'soc_0_30cm', 'soc_30_100cm', 'soc_100_200cm', 'mic_0_30cm', 
                         'mic_stock', 'enz_stock', 'doc_stock', 'poc_stock', 'soc_total', 
                         'bulk_A', 'bulk_I', 'bulk_K', 'bulk_V', 'bulk_Xi', 'bulk_E', 'bulk_NPP')
valid_loc_compas = which(is.na(da_summary[ , 'profile_num']) == 0)

# CLM5
model_name = 'cesm2_clm5_cen_vr_v2'
nn_exp_name = 'exp_pc_cesm2_23'
time_domain = 'whole_time'

eligible_loc = readMat(paste(data_dir_input, 'data4nn/eligible_profile_loc_0_cesm2_clm5_cen_vr_v2_whole_time.mat', sep = ''))
eligible_loc = eligible_loc$eligible.loc.0

mcmc_results_r2 = readMat(paste('/Users/ft254/DATAHUB/ENSEMBLE/OUTPUT_DATA/mcmc_summary_', model_name, '/', model_name, '_stat_r2.mat', sep = ''))
mcmc_results_r2 = apply(mcmc_results_r2$stat.r2, 1, max, na.rm = TRUE)

mcmc_results_gr = readMat(paste('/Users/ft254/DATAHUB/ENSEMBLE/OUTPUT_DATA/mcmc_summary_', model_name, '/', model_name, '_para_gr.mat', sep = ''))
mcmc_results_gr = apply(mcmc_results_gr$para.gr, 1, mean, na.rm = TRUE)

valid_loc_clm5 = intersect(eligible_loc, which(mcmc_results_r2 > 0.0 & mcmc_results_gr < 1.05))

#################################################################################
#  PRODA site prediction
#################################################################################
proda_clm5 = readMat(paste(data_dir_output, '/converged_soc/proda_explained_variation_cesm2_clm5_cen_vr_v2.mat', sep = ''))
proda_clm5 = proda_clm5$proda.explained.variation

proda_compas = readMat(paste(data_dir_output, '/converged_soc/proda_explained_variation_cesm2_clm5_mic_vr_v22.mat', sep = ''))
proda_compas = proda_compas$proda.explained.variation

wosis_obs = proda_clm5[[2]]
proda_clm5 = proda_clm5[[1]]
proda_compas = proda_compas[[1]]

current_data_clm5 = cbind(as.vector(wosis_obs[valid_loc_clm5, ]), as.vector(proda_clm5[valid_loc_clm5, ]))
current_data_compas = cbind(as.vector(wosis_obs[valid_loc_compas, ]), as.vector(proda_compas[valid_loc_compas, ]))

current_data_clm5 = current_data_clm5[is.na(apply(current_data_clm5, 1, mean)) == 0, ]
current_data_compas = current_data_compas[is.na(apply(current_data_compas, 1, mean)) == 0, ]

colnames(current_data_clm5) = c('obs', 'proda')
colnames(current_data_compas) = c('obs', 'proda')
# NSE calculation

nse_clm5 = 1 - sum((current_data_clm5[ , 'obs'] - current_data_clm5[ , 'proda'])**2)/sum((current_data_clm5[ , 'obs'] - mean(current_data_clm5[ , 'obs']))**2)
nse_compas = 1 - sum((current_data_compas[ , 'obs'] - current_data_compas[ , 'proda'])**2)/sum((current_data_compas[ , 'obs'] - mean(current_data_compas[ , 'obs']))**2)

#################################################################################
# plot figure
#################################################################################

# COMPAS
current_data_compas = data.frame(current_data_compas)/1000
text_data = data.frame('x_axis' = c(0.01), 
                       'y_axis' = c(1000), 
                       'nse' = paste('Explained variation = ', round(nse_compas*100, 2), '%', sep = ''))
p_compas =
  ggplot(data = current_data_compas) + 
    stat_bin_hex(aes(x = obs, y = proda), bins = 100) +
    geom_abline(intercept = 0, slope = 1, size = 1, color = 'black') + 
    scale_fill_gradientn(name = 'Count', colors = viridis(7), trans = 'identity', limits = c(1, 100), oob = scales::squish) +
    scale_x_continuous(limits = c(0.01, 250), trans = 'identity') + 
    scale_y_continuous(limits = c(0.01, 250), trans = 'identity') + 
    scale_x_continuous(limits = c(0.01, 1000), trans = 'log10', labels = trans_format('log10', math_format(10^.x))) +
    scale_y_continuous(limits = c(0.01, 1000), trans = 'log10', labels = trans_format('log10', math_format(10^.x))) +
  geom_text(data = text_data, aes(x = x_axis, y = y_axis, label = nse), vjust = 1, hjust = 0, size = 9, show.legend = FALSE) + 
    theme_classic() + 
    # add title
    labs(title = 'COMPAS', x = expression(paste('Observed SOC (kg C m'^'-3', ')', sep = '')), y = expression(paste('PRODA prediction (kg C m'^'-3', ')', sep = ''))) + 
    # change the legend properties
    guides(fill = guide_colorbar(direction = 'horizontal', barwidth = 15, barheight = 2.5, title.position = 'left', title.hjust = 0, title.vjust = 0.8, label.hjust = 0.5, frame.linewidth = 0), reverse = FALSE) +
    theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25))  +
    theme(legend.justification = c(0, 1), legend.position = 'none', legend.background = element_rect(fill = NA)) + 
    # modify the position of title
    theme(plot.title = element_text(hjust = 0.5, size = 50)) + 
    # modify the font size
    # modify the margin
    # theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
    theme(plot.margin = unit(c(0., 0.2, 0.2, 0.2), 'inch')) +
    theme(axis.text=element_text(size = 30, color = 'black'), axis.title = element_text(size = 30), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1, color = 'black'), axis.ticks.length = unit(0.12, 'inch')) 
  


# CLM5
current_data_clm5 = data.frame(current_data_clm5)/1000
text_data = data.frame('x_axis' = c(0.01), 
                       'y_axis' = c(1000), 
                       'nse' = paste('Explained variation = ', round(nse_clm5*100, 2), '%', sep = ''))
p_clm5 =
  ggplot(data = current_data_clm5) + 
  stat_bin_hex(aes(x = obs, y = proda), bins = 100) +
  geom_abline(intercept = 0, slope = 1, size = 1, color = 'black') + 
  scale_fill_gradientn(name = 'Count', colors = viridis(7), trans = 'identity', limits = c(1, 100), oob = scales::squish) +
  scale_x_continuous(limits = c(0.01, 1000), trans = 'log10', labels = trans_format('log10', math_format(10^.x))) +
  scale_y_continuous(limits = c(0.01, 1000), trans = 'log10', labels = trans_format('log10', math_format(10^.x))) +
  geom_text(data = text_data, aes(x = x_axis, y = y_axis, label = nse), vjust = 1, hjust = 0, size = 9, show.legend = FALSE) + 
  theme_classic() + 
  # add title
  labs(title = 'CLM5', x = expression(paste('Observed SOC (kg C m'^'-3', ')', sep = '')), y = expression(paste('PRODA prediction (kg C m'^'-3', ')', sep = ''))) + 
  # change the legend properties
  guides(fill = guide_colorbar(direction = 'horizontal', barwidth = 15, barheight = 2.5, title.position = 'left', title.hjust = 0, title.vjust = 0.3, label.position = 'top', label.hjust = 0.5, frame.linewidth = 0), reverse = FALSE) +
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25))  +
  theme(legend.justification = c(1, 0), legend.position = c(1, 0), legend.background = element_rect(fill = NA)) + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, size = 50)) + 
  # modify the font size
  # modify the margin
  # theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0., 0.2, 0.2, 0.2), 'inch')) +
  theme(axis.text=element_text(size = 30, color = 'black'), axis.title = element_text(size = 30), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1, color = 'black'), axis.ticks.length = unit(0.12, 'inch')) 

jpeg(paste('./Converged_SOC/posterior_model_prediction.jpeg', sep = ''), width = 16, height = 8, units = 'in', res = 300)

plot_grid(p_clm5, p_compas, 
          labels = c('a', 'b'),
          rel_widths = c(1, 1),
          label_size = 45,
          label_x = 0.02, label_y = 1.03,
          label_fontfamily = 'Arial',
          nrow = 1)
dev.off()

