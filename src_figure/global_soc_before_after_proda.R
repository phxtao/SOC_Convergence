## Packages
library(R.matlab)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(viridis)

library(GGally)
library(raster)

library(sf)
library(sp)
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

global_lat_lon_mic = readMat(paste(data_dir_output, 'converged_soc/soc_simu_grid_info_', model_name, '_', nn_exp_name, '_cross_valid_0_', as.character(1), '.mat', sep = ''))
global_lat_lon_mic = global_lat_lon_mic$var.data.middle[ , 1:2]
colnames(global_lat_lon_mic) = c('lon', 'lat')

soc_stock_mic = array(NA, dim = c(nrow(global_lat_lon_mic), 2, 10))

icross_valid = 2
for (icross_valid in 1:10) {
  global_simu = readMat(paste(data_dir_output, 'converged_soc/soc_simu_100cm_', model_name, '_', nn_exp_name, '_cross_valid_0_', as.character(icross_valid), '.mat', sep = ''))
  soc_stock_mic[ , , icross_valid] = global_simu$var.data.middle[ , c(1, 9)]/1000
}
soc_stock_mean_mic = apply(soc_stock_mic, c(1, 2), median, na.rm = TRUE)

#################################################################################
# soc stock clm5
#################################################################################
model_name = 'cesm2_clm5_cen_vr_v2'
nn_exp_name = 'exp_pc_cesm2_23'
time_domain = 'whole_time'

is_clm5_default_para = 1


global_lat_lon_clm = readMat(paste(data_dir_output, 'converged_soc/soc_simu_grid_info_', model_name, '_', nn_exp_name, '_cross_valid_0_', as.character(1), '.mat', sep = ''))
global_lat_lon_clm = global_lat_lon_clm$var.data.middle[ , 1:2]
colnames(global_lat_lon_clm) = c('lon', 'lat')

soc_stock_clm = array(NA, dim = c(nrow(global_lat_lon_clm), 2, 10))

icross_valid = 2
for (icross_valid in 1:10) {
  global_simu = readMat(paste(data_dir_output, 'converged_soc/soc_simu_100cm_', model_name, '_', nn_exp_name, '_cross_valid_0_', as.character(icross_valid), '.mat', sep = ''))
  if (is_clm5_default_para == 1) {
    soc_stock_clm[ , , icross_valid] = global_simu$var.data.middle[ , c(1, 3)]/1000
  } else {
    soc_stock_clm[ , , icross_valid] = global_simu$var.data.middle[ , c(1, 2)]/1000
  }
}
soc_stock_mean_clm = apply(soc_stock_clm, c(1, 2), median, na.rm = TRUE)

#################################################################################
# latitudinal soc stock
#################################################################################
resolution = 0.5
lat_seq = seq(from = 90 - resolution/2, to = -90 + resolution/2, by = -resolution)
lon_seq = seq(from = -180 + resolution/2, to = 180 - resolution/2, by = resolution)

col_name = c('lat', 'stock', 'lower', 'upper', 'model')
proda_lat_soc_stock = array(NA, dim = c(length(lat_seq)*2, length(col_name)))
adhoc_lat_soc_stock = array(NA, dim = c(length(lat_seq)*2, length(col_name)))
colnames(proda_lat_soc_stock) = col_name
colnames(adhoc_lat_soc_stock) = col_name

ilat = 100
for (ilat in 1:length(lat_seq)) {
  # mic
  proda_lat_soc_stock[ilat, 'lat'] = lat_seq[ilat]
  adhoc_lat_soc_stock[ilat, 'lat'] = lat_seq[ilat]
  profile_loc = which(global_lat_lon_mic[ , 'lat'] >= (lat_seq[ilat] - resolution/2) & global_lat_lon_mic[ , 'lat'] <= (lat_seq[ilat] + resolution/2))
  if (length(profile_loc) > 1) {
    # proda
    proda_lat_soc_stock[ilat, 'stock'] = median(soc_stock_mean_mic[profile_loc, 1], na.rm = TRUE)
    proda_lat_soc_stock[ilat, 'lower'] = quantile(soc_stock_mean_mic[profile_loc, 1], 0.16, na.rm = TRUE)
    proda_lat_soc_stock[ilat, 'upper'] = quantile(soc_stock_mean_mic[profile_loc, 1], 0.84, na.rm = TRUE)
    proda_lat_soc_stock[ilat, 'model'] = 1
    # ad hoc
    adhoc_lat_soc_stock[ilat, 'stock'] = median(soc_stock_mean_mic[profile_loc, 2], na.rm = TRUE)
    adhoc_lat_soc_stock[ilat, 'lower'] = quantile(soc_stock_mean_mic[profile_loc, 2], 0.16, na.rm = TRUE)
    adhoc_lat_soc_stock[ilat, 'upper'] = quantile(soc_stock_mean_mic[profile_loc, 2], 0.84, na.rm = TRUE)
    adhoc_lat_soc_stock[ilat, 'model'] = 1
  }
  
  # clm
  proda_lat_soc_stock[ilat + length(lat_seq), 'lat'] = lat_seq[ilat]
  adhoc_lat_soc_stock[ilat + length(lat_seq), 'lat'] = lat_seq[ilat]
  profile_loc = which(global_lat_lon_clm[ , 'lat'] >= (lat_seq[ilat] - resolution/2) & global_lat_lon_clm[ , 'lat'] <= (lat_seq[ilat] + resolution/2))
  if (length(profile_loc) > 1) {
    # proda
    proda_lat_soc_stock[ilat + length(lat_seq), 'stock'] = median(soc_stock_mean_clm[profile_loc, 1], na.rm = TRUE)
    proda_lat_soc_stock[ilat + length(lat_seq), 'lower'] = quantile(soc_stock_mean_clm[profile_loc, 1], 0.16, na.rm = TRUE)
    proda_lat_soc_stock[ilat + length(lat_seq), 'upper'] = quantile(soc_stock_mean_clm[profile_loc, 1], 0.84, na.rm = TRUE)
    proda_lat_soc_stock[ilat + length(lat_seq), 'model'] = 2
    
    # ad hoc
    adhoc_lat_soc_stock[ilat + length(lat_seq), 'stock'] = median(soc_stock_mean_clm[profile_loc, 2], na.rm = TRUE)
    adhoc_lat_soc_stock[ilat + length(lat_seq), 'lower'] = quantile(soc_stock_mean_clm[profile_loc, 2], 0.16, na.rm = TRUE)
    adhoc_lat_soc_stock[ilat + length(lat_seq), 'upper'] = quantile(soc_stock_mean_clm[profile_loc, 2], 0.84, na.rm = TRUE)
    adhoc_lat_soc_stock[ilat + length(lat_seq), 'model'] = 2
  }
}


proda_lat_soc_stock = data.frame(proda_lat_soc_stock)
adhoc_lat_soc_stock = data.frame(adhoc_lat_soc_stock)

proda_lat_soc_stock = proda_lat_soc_stock[is.na(proda_lat_soc_stock$model) == 0, ]
adhoc_lat_soc_stock = adhoc_lat_soc_stock[is.na(adhoc_lat_soc_stock$model) == 0, ]


### total soc storage
resolution = 0.5
lat_seq = global_lat_lon_clm[ , 2]
# area of grid 
radius = 6371008.8
length_top = (2*pi*radius*cos(abs(lat_seq+resolution/2)/180*pi)/360)*resolution
length_down = (2*pi*radius*cos(abs(lat_seq-resolution/2)/180*pi)/360)*resolution
height = (pi*radius/180)*resolution
lat_grid_area = (length_top + length_down)*height/2

# ad hoc
sum(soc_stock_mean_clm[ , 2]*lat_grid_area, na.rm = TRUE)/10**12 # unit Pg
sum(soc_stock_mean_mic[ , 2]*lat_grid_area, na.rm = TRUE)/10**12 # unit Pg
# proda
sum(soc_stock_mean_clm[ , 1]*lat_grid_area, na.rm = TRUE)/10**12 # unit Pg
sum(soc_stock_mean_mic[ , 1]*lat_grid_area, na.rm = TRUE)/10**12 # unit Pg

# ad hoc
cor.test(log10(soc_stock_mean_clm[ , 2]), log10(soc_stock_mean_mic[ , 2]))
# proda
cor.test(log10(soc_stock_mean_clm[ , 1]), log10(soc_stock_mean_mic[ , 1]))

#################################################################################
# plot figures
#################################################################################
line_color = c('#994F00', '#006CD1')
line_name = c('COMPAS', 'CLM5')
p_lat_soc_adhoc =
  ggplot(data = adhoc_lat_soc_stock) +
  geom_ribbon(aes(x = lat, ymin = lower, ymax = upper, group = as.factor(model), fill = as.factor(model)), alpha = 0.3) +
  geom_line(aes(x = lat, y = stock, group = as.factor(model), color = as.factor(model)), alpha = 1, size = 2) +
  scale_color_manual(name = '', labels = line_name, values = line_color) +
  scale_fill_manual(name = '', labels = line_name, values = line_color) +
  scale_x_continuous(n.breaks = 7) + 
  scale_y_continuous(trans = 'identity') +
    coord_flip(ylim = c(0, 50)) +
  # change the background to black and white
  theme_classic() +
  # add title
  labs(title = '', x = 'Latitude (º)', y = expression(paste('SOC stock (kg C m'^'-2', ')'))) +
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, size = 50)) + 
  theme(legend.justification = c(0, 1), legend.position = c(0.05, 1.05), legend.background = element_rect(fill = NA), legend.text.align = 0) +
  theme(legend.text = element_text(size = 30), legend.title = element_text(size = 30))  +
  theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.5, 'inch')) +
  theme(legend.direction = 'horizontal') + 
  # modify the margin
  # theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0., 0.2, 0.2, 0.2), 'inch')) +
  theme(axis.text=element_text(size = 30, color = 'black'), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1, color = 'black'), axis.ticks.length = unit(0.12, 'inch')) 

p_lat_soc_proda =
  ggplot(data = proda_lat_soc_stock) +
  geom_ribbon(aes(x = lat, ymin = lower, ymax = upper, group = as.factor(model), fill = as.factor(model)), alpha = 0.3) +
  geom_line(aes(x = lat, y = stock, group = as.factor(model), color = as.factor(model)), alpha = 1, size = 2) +
  scale_color_manual(name = '', labels = line_name, values = line_color) +
  scale_fill_manual(name = '', labels = line_name, values = line_color) +
  scale_x_continuous(n.breaks = 7) + 
  scale_y_continuous(trans = 'identity') +
  coord_flip(ylim = c(0, 50)) +
    # change the background to black and white
    theme_classic() +
    # add title
    labs(title = '', x = 'Latitude (º)', y = expression(paste('SOC stock (kg C m'^'-2', ')'))) +
    # modify the position of title
    theme(plot.title = element_text(hjust = 0.5, size = 50)) + 
    theme(legend.justification = c(0, 1), legend.position = c(0.05, 1.05), legend.background = element_rect(fill = NA), legend.text.align = 0) +
    theme(legend.text = element_text(size = 30), legend.title = element_text(size = 30))  +
    theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.5, 'inch')) +
    theme(legend.direction = 'horizontal') + 
    # modify the margin
    # theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
    theme(plot.margin = unit(c(0., 0.2, 0.2, 0.2), 'inch')) +
    theme(axis.text=element_text(size = 30, color = 'black'), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1, color = 'black'), axis.ticks.length = unit(0.12, 'inch')) 
  
  

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

# soc map mic
current_data_mic = data.frame(cbind(global_lat_lon_mic, soc_stock_mean_mic))
colnames(current_data_mic) = c('lon', 'lat', 'soc_proda', 'soc_adhoc')

lon_lat_transfer = project(xy = as.matrix(current_data_mic[ , c('lon', 'lat')]), proj = '+proj=robin') 
current_data_mic[ , c('lon', 'lat')] = lon_lat_transfer

lat_limits = rbind(c(0, -56), c(0, 80))
lat_limits_robin = project(xy = as.matrix(lat_limits), proj = '+proj=robin') 

p_soc_mic_adhoc =
  ggplot() +
  geom_tile(data = current_data_mic, aes(x = lon, y = lat, fill = soc_adhoc), height = 60000, width = 60000, na.rm = TRUE) +
  scale_fill_gradientn(name = expression(paste('kg C m'^'-2', sep = '')), colours = rev(viridis(15)), na.value="transparent", limits = c(5, 50), breaks = c(5, 15, 50), trans = 'log10', oob = scales::squish) +
  geom_sf(data = world_coastline, fill = NA, color = 'black', linewidth = 1) + 
  geom_polygon(data = world_ocean, aes(x = long, y = lat), fill = NA, color = 'black', size = 2) +
  # change the background to black and white
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
  labs(title = 'COMPAS', x = '', y = '') + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, vjust = -1, size = 50)) + 
  # modify the font size
  theme(axis.title = element_text(size = 20)) + 
  theme(panel.background = element_rect(fill = NA, colour = NA)) +
  # modify the margin
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0, 0, 0, 0), 'inch')) +
  theme(axis.text=element_text(size = 35, color = 'black'))


p_soc_mic_proda =
  ggplot() +
  geom_tile(data = current_data_mic, aes(x = lon, y = lat, fill = soc_proda), height = 60000, width = 60000, na.rm = TRUE) +
  scale_fill_gradientn(name = expression(paste('kg C m'^'-2', sep = '')), colours = rev(viridis(15)), na.value="transparent", limits = c(5, 50), breaks = c(5, 15, 50), trans = 'log10', oob = scales::squish) +
  geom_sf(data = world_coastline, fill = NA, color = 'black', linewidth = 1) + 
  geom_polygon(data = world_ocean, aes(x = long, y = lat), fill = NA, color = 'black', size = 2) +
  # change the background to black and white
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
  labs(title = 'COMPAS', x = '', y = '') + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, vjust = -1, size = 50)) + 
  # modify the font size
  theme(axis.title = element_text(size = 20)) + 
  theme(panel.background = element_rect(fill = NA, colour = NA)) +
  # modify the margin
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0, 0, 0, 0), 'inch')) +
  theme(axis.text=element_text(size = 35, color = 'black'))



# soc map clm
current_data_clm = data.frame(cbind(global_lat_lon_clm, soc_stock_mean_clm))
colnames(current_data_clm) = c('lon', 'lat', 'soc_proda', 'soc_adhoc')

lon_lat_transfer = project(xy = as.matrix(current_data_clm[ , c('lon', 'lat')]), proj = '+proj=robin') 
current_data_clm[ , c('lon', 'lat')] = lon_lat_transfer

lat_limits = rbind(c(0, -56), c(0, 80))
lat_limits_robin = project(xy = as.matrix(lat_limits), proj = '+proj=robin') 

p_soc_clm_adhoc =
  ggplot() +
  geom_tile(data = current_data_clm, aes(x = lon, y = lat, fill = soc_adhoc), height = 60000, width = 60000, na.rm = TRUE) +
  scale_fill_gradientn(name = expression(paste('kg C m'^'-2', sep = '')), colours = rev(viridis(15)), na.value="transparent", limits = c(5, 50), breaks = c(5, 15, 50), trans = 'log10', oob = scales::squish) +
  geom_sf(data = world_coastline, fill = NA, color = 'black', linewidth = 1) + 
  geom_polygon(data = world_ocean, aes(x = long, y = lat), fill = NA, color = 'black', size = 2) +
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
  labs(title = 'CLM5', x = '', y = '') + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, vjust = -1, size = 50)) + 
  # modify the font size
  theme(axis.title = element_text(size = 20)) + 
  theme(panel.background = element_rect(fill = NA, colour = NA)) +
  # modify the margin
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0, 0, 0, 0), 'inch')) +
  theme(axis.text=element_text(size = 35, color = 'black'))


p_soc_clm_proda =
  ggplot() +
  geom_tile(data = current_data_clm, aes(x = lon, y = lat, fill = soc_proda), height = 60000, width = 60000, na.rm = TRUE) +
  scale_fill_gradientn(name = expression(paste('kg C m'^'-2', sep = '')), colours = rev(viridis(15)), na.value="transparent", limits = c(5, 50), breaks = c(5, 15, 50), trans = 'log10', oob = scales::squish) +
  geom_sf(data = world_coastline, fill = NA, color = 'black', linewidth = 1) + 
  geom_polygon(data = world_ocean, aes(x = long, y = lat), fill = NA, color = 'black', size = 2) +
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
  labs(title = 'CLM5', x = '', y = '') + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, vjust = -1, size = 50)) + 
  # modify the font size
  theme(axis.title = element_text(size = 20)) + 
  theme(panel.background = element_rect(fill = NA, colour = NA)) +
  # modify the margin
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0, 0, 0, 0), 'inch')) +
  theme(axis.text=element_text(size = 35, color = 'black'))

#######################################
# ad hoc
p_map_adhoc = plot_grid(p_soc_clm_adhoc, p_soc_mic_adhoc,
          labels = c('a', 'b'),
          nrow = 2,
          label_size = 50,
          label_x = 0.02, label_y = 1.0,
          label_fontface = 'bold'
          )

jpeg(paste('./Converged_SOC/ad_hoc_soc_maps_clm5_default_', is_clm5_default_para,  '.jpeg', sep = ''), width = 20, height = 12, units = 'in', res = 300)

plot_grid(p_map_adhoc, p_lat_soc_adhoc, NULL,
          labels = c('', 'c'),
          nrow = 1,
          rel_widths = c(2, 1, 0.2),
          label_size = 50,
          label_x = 0.02, label_y = 1.0,
          label_fontfamily = 'Arial',
          label_fontface = 'bold'
          )
dev.off()


# proda
p_map_proda = plot_grid(p_soc_clm_proda, p_soc_mic_proda,
                        labels = c('a', 'b'),
                        nrow = 2,
                        label_size = 50,
                        label_x = 0.02, label_y = 1.0,
                        label_fontface = 'bold'
)

jpeg(paste('./Converged_SOC/proda_soc_maps.jpeg', sep = ''), width = 20, height = 12, units = 'in', res = 300)

plot_grid(p_map_proda, p_lat_soc_proda, NULL,
          labels = c('', 'c'),
          nrow = 1,
          rel_widths = c(2, 1, 0.2),
          label_size = 50,
          label_x = 0.02, label_y = 1.0,
          label_fontfamily = 'Arial',
          label_fontface = 'bold'
)
dev.off()
