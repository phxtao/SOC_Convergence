
library(ncdf4)
library(R.matlab)
library(ggplot2)
# library(qlcMatrix)
library(jcolors)
library(cowplot)

# dev.off()
rm(list = ls())

setwd('/Users/phoenix/Google_Drive/R/Converged_SOC')

#######################################
# microbial model
#######################################
model_name = 'cesm2_clm5_mic_vr_v22'

ncfname = '/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/ENSEMBLE/INPUT_DATA/wosis_2019_snap_shot/soc_profile_wosis_2019_snapshot_hugelius_mishra.nc'
profile_info = nc_open(ncfname)
profile_info = ncvar_get(profile_info)

max_depth = readMat('/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/ENSEMBLE/INPUT_DATA/wosis_2019_snap_shot/wosis_2019_snapshot_hugelius_mishra_env_info.mat')
max_depth = max_depth$EnvInfo[ , 73]

mcmc_results_r2 = readMat(paste('/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/ENSEMBLE/OUTPUT_DATA/mcmc_summary_', model_name, '/', model_name, '_stat_r2.mat', sep = ''))
mcmc_results_r2 = mcmc_results_r2$stat.r2

mcmc_results_gr = readMat(paste('/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/ENSEMBLE/OUTPUT_DATA/mcmc_summary_', model_name, '/', model_name, '_para_gr.mat', sep = ''))
mcmc_results_gr = mcmc_results_gr$para.gr

mcmc_results_std = readMat(paste('/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/ENSEMBLE/OUTPUT_DATA/mcmc_summary_', model_name, '/', model_name, '_para_var.mat', sep = ''))
mcmc_results_std = mcmc_results_std$para.var


current_data_mic = data.frame(cbind(profile_info[ , c(1, 4, 5, 6)], apply(mcmc_results_r2, 1, max, na.rm = TRUE), apply(mcmc_results_gr, 1, mean, na.rm = TRUE), max_depth, NA))
colnames(current_data_mic) = c('id', 'lon', 'lat', 'layer_num', 'r2', 'gr', 'max_depth', 'source')

#######################################
# CLM5
#######################################
model_name = 'cesm2_clm5_cen_vr_v2'

ncfname = '/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/ENSEMBLE/INPUT_DATA/wosis_2019_snap_shot/soc_profile_wosis_2019_snapshot_hugelius_mishra.nc'
profile_info = nc_open(ncfname)
profile_info = ncvar_get(profile_info)

max_depth = readMat('/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/ENSEMBLE/INPUT_DATA/wosis_2019_snap_shot/wosis_2019_snapshot_hugelius_mishra_env_info.mat')
max_depth = max_depth$EnvInfo[ , 73]

mcmc_results_r2 = readMat(paste('/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/ENSEMBLE/OUTPUT_DATA/mcmc_summary_', model_name, '/', model_name, '_stat_r2.mat', sep = ''))
mcmc_results_r2 = mcmc_results_r2$stat.r2

mcmc_results_gr = readMat(paste('/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/ENSEMBLE/OUTPUT_DATA/mcmc_summary_', model_name, '/', model_name, '_para_gr.mat', sep = ''))
mcmc_results_gr = mcmc_results_gr$para.gr

mcmc_results_std = readMat(paste('/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/ENSEMBLE/OUTPUT_DATA/mcmc_summary_', model_name, '/', model_name, '_para_var.mat', sep = ''))
mcmc_results_std = mcmc_results_std$para.var


current_data_clm5 = data.frame(cbind(profile_info[ , c(1, 4, 5, 6)], apply(mcmc_results_r2, 1, max, na.rm = TRUE), apply(mcmc_results_gr, 1, mean, na.rm = TRUE), max_depth, NA))
colnames(current_data_clm5) = c('id', 'lon', 'lat', 'layer_num', 'r2', 'gr', 'max_depth', 'source')



## Jet colorbar function
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
# can be changed to state or world to have US and world map
world_coastline = rgdal::readOGR(dsn='/Users/phoenix/Google_Drive/Tsinghua_Luo/World_Vector_Shape/ne110m/ne_110m_land.shp',layer = 'ne_110m_land')
world_coastline <- fortify(world_coastline)
Map.Using = world_coastline


valid_site_wosis = which((current_data_mic$id > 10000 & current_data_mic$r2 > 0.0 & current_data_mic$max_depth > 50 & current_data_mic$gr < 1.05) | 
                           (current_data_clm5$id > 10000 & current_data_clm5$r2 > 0.0 & current_data_clm5$max_depth > 50 & current_data_clm5$gr < 1.05))


valid_site_mishra = which((current_data_mic$id < 10000 & current_data_mic$r2 > 0.0 & current_data_mic$max_depth > 50 & current_data_mic$gr < 1.05) | 
                            (current_data_clm5$id < 10000 & current_data_clm5$r2 > 0.0 & current_data_clm5$max_depth > 50 & current_data_clm5$gr < 1.05))

invalid_site = which((current_data_mic$r2 <= 0.0 | current_data_mic$max_depth <= 50 | current_data_mic$gr >= 1.05) |
                       (current_data_clm5$r2 <= 0.0 | current_data_clm5$max_depth <= 50 | current_data_clm5$gr >= 1.05))

current_data_mic$source[valid_site_wosis] = 'A_WoSIS'
current_data_mic$source[valid_site_mishra] = 'B_Mishra & Hugelius 2020'
current_data_mic$source[invalid_site] = 'C_Invalid'


current_data = current_data_mic[order(current_data_mic$source, decreasing = TRUE), ]


jpeg(paste('./profile_loc.jpeg', sep = ''), width = 10, height = 5, units = 'in', res = 300)
# p_wosis = 
  ggplot(data = current_data) +
  geom_point(aes(x = lon, y = lat, color = source), shape = 16, size = 0.2, alpha = 1) + 
  scale_color_manual(name = '', values = c('#005AB5', '#DC3220', 'grey'), labels = c('WoSIS 2019 Snapshot', 'Mishra & Hugelius 2020', 'Invalid')) +
  geom_polygon(data = Map.Using, aes(x = long, y = lat, group = group), fill = NA, color = 'black', size = 0.3) +
  ylim(c(-56, 80)) +
  # change the background to black and white
  theme_bw() +
  # change the legend properties
  # theme(legend.position = 'none') +
  theme(legend.justification = c(0, 0), legend.position = c(0, 0), legend.background = element_rect(fill = NA), legend.text.align = 0) +
  theme(legend.text = element_text(size = 35), legend.title = element_text(size = 35))  +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  theme(legend.text = element_text(size = 15), legend.title = element_text(size = 20)) +
  # add title
  labs(x = '', y = '') + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, size = 40)) + 
  # modify the font size
  theme(axis.title = element_text(size = 20)) + 
  # modify the margin
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'inch')) +
  theme(axis.text=element_text(size = 30))

dev.off()
