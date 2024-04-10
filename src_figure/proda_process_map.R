## Packages
library(R.matlab)
library(ggplot2)
library(cowplot)
# library(jcolors)
library(gridExtra)
library(viridis)
library(sf)
library(sp)
library(GGally)
library(raster)
library(proj4)
##
rm(list = ls())

setwd('/Users/ft254/Google_Drive/R')


## Jet colorbar function
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
diff.colors <- colorRampPalette(c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "white", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"))

#############################################################################
# Data Path
#############################################################################
data_dir_output = '/Users/ft254/DATAHUB/ENSEMBLE/OUTPUT_DATA/'
data_dir_input = '/Users/ft254/DATAHUB/ENSEMBLE/INPUT_DATA/'

#############################################################################
# function to increase vertical spacing between legend keys
#############################################################################
# @clauswilke
draw_key_polygon3 <- function(data, params, size) {
  lwd <- min(data$size, min(size) / 4)
  
  grid::rectGrob(
    width = grid::unit(0.6, "npc"),
    height = grid::unit(0.6, "npc"),
    gp = grid::gpar(
      col = data$colour,
      fill = alpha(data$fill, data$alpha),
      lty = data$linetype,
      lwd = lwd * .pt,
      linejoin = "mitre"
    ))
}

# register new key drawing function, 
# the effect is global & persistent throughout the R session
GeomBar$draw_key = draw_key_polygon3


#################################################################################
# soc stock microbial model
#################################################################################
model_name = 'cesm2_clm5_mic_vr_v22'
nn_exp_name = 'exp_pc_cesm2_23'
time_domain = 'whole_time'

valid_grid_loc_mic = read.csv(paste(data_dir_output, 'neural_networking/valid_grid_loc_', model_name, '_', time_domain, '_', nn_exp_name, '_cross_valid_0_', as.character(1), '.csv', sep = ''), header = FALSE)
valid_grid_loc_mic = valid_grid_loc_mic$V1

global_lat_lon_mic = readMat(paste(data_dir_output, 'converged_soc/soc_simu_grid_info_', model_name, '_', nn_exp_name, '_cross_valid_0_', as.character(1), '.mat', sep = ''))
global_lat_lon_mic = global_lat_lon_mic$var.data.middle[ , 1:2]
colnames(global_lat_lon_mic) = c('lon', 'lat')

bulk_process_mic = array(NA, dim = c(nrow(global_lat_lon_mic), 7, 10))

icross_valid = 2
for (icross_valid in 1:10) {
  global_simu = readMat(paste(data_dir_output, 'converged_soc/bulk_process_summary_', model_name, '_', nn_exp_name, '_cross_valid_0_', as.character(icross_valid), '.mat', sep = ''))
  bulk_process_mic[ , , icross_valid] = global_simu$var.data.middle[ , c(1:7)]
}
bulk_process_mean_mic = apply(bulk_process_mic, c(1, 2), median, na.rm = TRUE)

#################################################################################
# soc stock clm5
#################################################################################
model_name = 'cesm2_clm5_cen_vr_v2'
nn_exp_name = 'exp_pc_cesm2_23'
time_domain = 'whole_time'

valid_grid_loc_clm = read.csv(paste(data_dir_output, 'neural_networking/valid_grid_loc_', model_name, '_', time_domain, '_', nn_exp_name, '_cross_valid_0_', as.character(1), '.csv', sep = ''), header = FALSE)
valid_grid_loc_clm = valid_grid_loc_clm$V1

global_lat_lon_clm = readMat(paste(data_dir_output, 'converged_soc/soc_simu_grid_info_', model_name, '_', nn_exp_name, '_cross_valid_0_', as.character(1), '.mat', sep = ''))
global_lat_lon_clm = global_lat_lon_clm$var.data.middle[ , 1:2]
colnames(global_lat_lon_clm) = c('lon', 'lat')

bulk_process_clm = array(NA, dim = c(nrow(global_lat_lon_clm), 7, 10))

icross_valid = 2
for (icross_valid in 1:10) {
  global_simu = readMat(paste(data_dir_output, 'converged_soc/bulk_process_summary_', model_name, '_', nn_exp_name, '_cross_valid_0_', as.character(icross_valid), '.mat', sep = ''))
  bulk_process_clm[ , , icross_valid] = global_simu$var.data.middle[ , c(1:7)]
}
bulk_process_mean_clm = apply(bulk_process_clm, c(1, 2), median, na.rm = TRUE)

#################################################################################
# plot figures
#################################################################################
#-------------------------------------soc stock and Residence Time map
world_coastline = st_read('/Users/ft254/Google_Drive/Tsinghua_Luo/World_Vector_Shape/ne110m/ne_110m_land.shp', layer = 'ne_110m_land')
world_coastline <- st_transform(world_coastline, CRS('+proj=robin'))

ocean_left = cbind(rep(-180, 100), seq(from = 80, to = -56, by = -(80 + 56)/(100 -1)))
ocean_right = cbind(rep(180, 100), seq(from = -56, to = 80, by = (80 + 56)/(100 -1)))
ocean_top = cbind(seq(from = 180, to = -180, by = -(360)/(100 -1)), rep(80, 100))
ocean_bottom = cbind(seq(from = -180, to = 180, by = (360)/(100 -1)), rep(-56, 100))

world_ocean = rbind(ocean_left, ocean_bottom, ocean_right, ocean_top)
world_ocean = as.matrix(world_ocean)

world_ocean <- project(xy = world_ocean, proj = '+proj=robin')

world_ocean = data.frame(world_ocean)
colnames(world_ocean) = c('long', 'lat')

# process map mic
current_data_mic = data.frame(cbind(global_lat_lon_mic, bulk_process_mean_mic))
colnames(current_data_mic) = c('lon', 'lat', 'A', 'I', 'K', 'V', 'Xi', 'NPP', 'F')
lon_lat_transfer = project(xy = as.matrix(current_data_mic[ , c('lon', 'lat')]), proj = '+proj=robin') 
current_data_mic[ , c('lon', 'lat')] = lon_lat_transfer
# process map clm
current_data_clm = data.frame(cbind(global_lat_lon_clm, bulk_process_mean_clm))
colnames(current_data_clm) = c('lon', 'lat', 'A', 'I', 'K', 'V', 'Xi', 'NPP', 'F')
lon_lat_transfer = project(xy = as.matrix(current_data_clm[ , c('lon', 'lat')]), proj = '+proj=robin') 
current_data_clm[ , c('lon', 'lat')] = lon_lat_transfer

lat_limits = rbind(c(0, -56), c(0, 80))
lat_limits_robin = project(xy = as.matrix(lat_limits), proj = '+proj=robin') 

####################
process_scale_option = c('identity', 'identity', 'log10', 'identity', 'identity', 'identity', 'identity', 'identity')

process_name =  c('Carbon transfer', 
                  'Carbon input allocation', 
                  'Baseline decomposition', 
                  'Vertical transport rate',
                  'Environmental modifer', 
                  'Plant carbon inputs', 
                  'Litter to mineral soil fraction')
process_unit = c('unitless',
                 'unitless', 
                 expression(paste('yr'^'-1', sep = '')),
                 expression(paste('yr'^'-1', sep = '')),
                 'unitless', 
                 expression(paste('g C m'^'-2', ' yr'^'-1', sep = '')), 
                 'unitless')


corr_process_summary = array(NA, dim = c(length(process_name), 2))

ipara = 3
## Projected para 
for (ipara in 1:length(process_name)){
  
  middle_data_mic = current_data_mic[ , c(1, 2, (ipara+2))]
  middle_data_clm = current_data_clm[ , c(1, 2, (ipara+2))]
  colnames(middle_data_mic) = c('lon', 'lat', 'process')
  colnames(middle_data_clm) = c('lon', 'lat', 'process')
  
  legend_lower_mic = apply(current_data_mic[ , c(3:9)], 2, quantile, prob = 0.05 , na.rm = TRUE)
  legend_upper_mic = apply(current_data_mic[ , c(3:9)], 2, quantile, prob = 0.95, na.rm = TRUE)
  
  legend_lower_clm = apply(current_data_clm[ , c(3:9)], 2, quantile, prob = 0.05, na.rm = TRUE)
  legend_upper_clm = apply(current_data_clm[ , c(3:9)], 2, quantile, prob = 0.95, na.rm = TRUE)
  
  legend_lower_clm = apply(rbind(legend_lower_clm, legend_lower_mic), 2, min)
  legend_upper_clm = apply(rbind(legend_upper_clm, legend_upper_mic), 2, max)
  
  legend_lower_mic = legend_lower_clm
  legend_upper_mic = legend_upper_clm
  
  p_mic =
    ggplot() +
    geom_tile(data = middle_data_mic, aes(x = lon, y = lat, fill = process), height = 60000, width = 60000, na.rm = TRUE) +
    scale_fill_gradientn(name = process_unit[ipara], colours = rev(viridis(15)), na.value="transparent", limits = c(legend_lower_mic[ipara], legend_upper_mic[ipara]), trans = process_scale_option[ipara], oob = scales::squish) +
    geom_sf(data = world_coastline, fill = NA, color = 'black', linewidth = 1) + 
    geom_polygon(data = world_ocean, aes(x = long, y = lat), fill = NA, color = 'black', size = 2) +
    # change the background to black and white
    # coord_equal() +
    # theme_map() +
    ylim(lat_limits_robin[ , 2]) +
    # change the legend properties
    # theme(legend.position = 'none') +
    theme(legend.justification = c(0, 0), legend.position = c(0.03, 0.02), legend.background = element_rect(fill = NA), legend.text.align = 0) +
    # theme(legend.justification = c(0.5, 0), legend.position = c(0.5, 0), legend.background = element_rect(fill = NA), legend.direction = 'horizontal') +
    # change the size of colorbar
    guides(fill = guide_colorbar(direction = 'vertical', barwidth = 2.5, barheight = 14, title.position = 'top', title.hjust = 0, label.hjust = 0, frame.linewidth = 0), reverse = FALSE) +
    theme(legend.text = element_text(size = 35, ), legend.title = element_text(size = 40)) +
    # add title
    labs(title = paste('COMPAS: ', process_name[ipara], sep = ''), x = '', y = '') + 
    # modify the position of title
    theme(plot.title = element_text(hjust = 0.5, vjust = -1, size = 40)) + 
    # modify the font size
    theme(axis.title = element_text(size = 20)) + 
    theme(panel.background = element_rect(fill = NA, colour = NA)) +
    # modify the margin
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
    theme(plot.margin = unit(c(0, 0, 0, 0), 'inch')) +
    theme(axis.text=element_text(size = 35, color = 'black'))
  
  
  p_clm =
    ggplot() +
    geom_tile(data = middle_data_clm, aes(x = lon, y = lat, fill = process), height = 60000, width = 60000, na.rm = TRUE) +
    scale_fill_gradientn(name = process_unit[ipara], colours = rev(viridis(15)), na.value="transparent", limits = c(legend_lower_clm[ipara], legend_upper_clm[ipara]), trans = process_scale_option[ipara], oob = scales::squish) +
    geom_sf(data = world_coastline, fill = NA, color = 'black', linewidth = 1) + 
    geom_polygon(data = world_ocean, aes(x = long, y = lat), fill = NA, color = 'black', size = 2) +
    # change the background to black and white
    # coord_equal() +
    # theme_map() +
    ylim(lat_limits_robin[ , 2]) +
    # change the legend properties
    # theme(legend.position = 'none') +
    theme(legend.justification = c(0, 0), legend.position = c(0.03, 0.02), legend.background = element_rect(fill = NA), legend.text.align = 0) +
    # theme(legend.justification = c(0.5, 0), legend.position = c(0.5, 0), legend.background = element_rect(fill = NA), legend.direction = 'horizontal') +
    # change the size of colorbar
    guides(fill = guide_colorbar(direction = 'vertical', barwidth = 2.5, barheight = 14, title.position = 'top', title.hjust = 0, label.hjust = 0, frame.linewidth = 0), reverse = FALSE) +
    theme(legend.text = element_text(size = 35, ), legend.title = element_text(size = 40)) +
    # add title
    labs(title = paste('CLM5: ', process_name[ipara], sep = ''), x = '', y = '') + 
    # modify the position of title
    theme(plot.title = element_text(hjust = 0.5, vjust = -1, size = 40)) + 
    # modify the font size
    theme(axis.title = element_text(size = 20)) + 
    theme(panel.background = element_rect(fill = NA, colour = NA)) +
    # modify the margin
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
    theme(plot.margin = unit(c(0, 0, 0, 0), 'inch')) +
    theme(axis.text=element_text(size = 35, color = 'black'))
  
  # correlation
  middle_data_corr = cbind(middle_data_clm$process, middle_data_mic$process)
  middle_data_corr = data.frame(middle_data_corr)
  colnames(middle_data_corr) = c('clm', 'mic')
  
  if (ipara == 3) {
    limit_clm = c(0.001, 1)
    limit_mic = c(0.001, 1)
  } else {
    limit_clm = quantile(middle_data_clm$process, probs = c(0, 1), na.rm = TRUE)
    limit_mic = quantile(middle_data_mic$process, probs = c(0, 1), na.rm = TRUE)
  }
  limit_clm = c(min(limit_clm[1], limit_mic[1]), max(limit_clm[2], limit_mic[2]))
  limit_mic = limit_clm
  
  corr_process_middle = cor.test(middle_data_corr$clm, middle_data_corr$mic, na.rm = TRUE)
  corr_process_summary[ipara, ] = c(corr_process_middle$estimate, corr_process_middle$p.value)
  p_corr =
    ggplot() + 
    stat_bin_hex(data = middle_data_corr, aes(x = clm, y = mic), bins = 100) +
    scale_fill_gradientn(name = 'Count', colors = viridis(7), trans = 'identity', oob = scales::squish) +
    geom_abline(slope = 1, intercept = 0, size = 2, color = 'black') +
    scale_x_continuous(trans = process_scale_option[ipara], limits = limit_clm) +
    scale_y_continuous(trans = process_scale_option[ipara], limits = limit_mic) +
    theme_classic() + 
    # add title
    labs(title = '', x = 'CLM5', y = 'COMPAS') + 
    # change the legend properties
    guides(fill = guide_colorbar(direction = 'horizontal', barwidth = 15, barheight = 2.5, title.position = 'right', title.hjust = 0, title.vjust = 0.8, label.hjust = 0.5, frame.linewidth = 0), reverse = FALSE) +
    theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25))  +
    theme(legend.justification = c(1, 0), legend.position = 'None', legend.background = element_rect(fill = NA)) + 
    # modify the position of title
    theme(plot.title = element_text(hjust = 0.5, size = 50)) + 
    # modify the font size
    # modify the margin
    # theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
    theme(plot.margin = unit(c(0., 0.2, 0.2, 0.2), 'inch')) +
    theme(axis.text=element_text(size = 30, color = 'black'), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1, color = 'black'), axis.ticks.length = unit(0.12, 'inch')) 
  
  eval(parse(text = paste('p_corr', ipara, ' = p_corr', sep = '')))
  eval(parse(text = paste('p_clm', ipara, ' = p_clm', sep = '')))
  eval(parse(text = paste('p_mic', ipara, ' = p_mic', sep = '')))
}


jpeg(paste('./Converged_SOC/bulk_process.jpeg', sep = ''), width = 30, height = 35, units = 'in', res = 300)
plot_grid(p_clm1, p_mic1, p_corr1, NULL, 
          p_clm3, p_mic3, p_corr3, NULL, 
          p_clm5, p_mic5, p_corr5, NULL, 
          p_clm2, p_mic2, p_corr2, NULL, 
          p_clm4, p_mic4, p_corr4, NULL, 
          p_clm6, p_mic6, p_corr6, NULL, 
          nrow = 6, ncol = 4 ,
          rel_widths = c(3, 3, 2, 0.10),
          labels = c('a', 'b', 'c', ' ',
                     'd', 'e', 'f', ' ',
                     'g', 'h', 'i', ' ',
                     'j', 'k', 'l', ' ',
                     'm', 'n', 'o', ' ',
                     'p', 'q', 'r', ' '),
          label_size = 70,
          label_x = 0.05, label_y = 1.05,
          label_fontfamily = 'Arial',
          label_fontface = 'bold'
)
dev.off()

jpeg(paste('./Converged_SOC/litter_to_soil_fraction.jpeg', sep = ''), width = 30, height = 6, units = 'in', res = 300)
plot_grid(p_clm7, p_mic7, p_corr7, NULL, 
          nrow = 1, ncol = 4 ,
          rel_widths = c(3, 3, 2, 0.10),
          labels = c('a', 'b', 'c', ' '),
          label_size = 70,
          label_x = 0.05, label_y = 1.05,
          label_fontfamily = 'Arial',
          label_fontface = 'bold'
)
dev.off()




