## Packages
library(R.matlab)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(viridis)
library(ncdf4)
library(GGally)

library(raster)

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
# site level data assimilation summary
#################################################################################
mic_model_name = 'cesm2_clm5_mic_vr_v22'
clm_model_name = 'cesm2_clm5_cen_vr_v2'
nn_exp_name = 'exp_pc_cesm2_23'
time_domain = 'whole_time'

# para names
para_list_clm =  c('diffus', 'cryo', 
                   'q10', 'efolding', 
                   'taucwd', 'taul1', 'taul2', 'tau4s1', 'tau4s2', 'tau4s3', 
                   'fs1l1', 'fs1l2', 'fs2l3', 'fs2s1', 'fs3s1', 'fs1s2', 'fs3s2', 'fs1s3', 'fl2cwd', 
                   'w-scaling', 'b')
default_clm5 = c(0.0001, 0.0005, 
                 1.5, NA, 
                 3.33, 0.0541, 0.2041, 0.1370, 5, 222.222, 
                 0.42, 0.5, 0.5, NA, NA, 0.42, 0.03, 0.45, 0.786, 
                 1, NA)
prior_lower_clm5 = c(0.00003, 0.00003, 
                     1.2, 0, 
                     1, 0, 0.1, 0, 1, 200, 
                     0.1, 0.2, 0.2, 0, 0, 0.1, 0, 0, 0.5, 
                     0, 0.5)


prior_upper_clm5 = c(0.0005, 0.0016, 
                     3, 1, 
                     6, 0.11, 0.3, 1, 50, 1000, 
                     0.8, 0.8, 0.8, 0.4, 0.05, 0.74, 0.1, 0.9, 1, 
                     5, 1)

scaled_default_clm5 = (default_clm5 - prior_lower_clm5)/(prior_upper_clm5 - prior_lower_clm5)

para_list_mic =  c('diffus', 'cryo', 'q10', 'w_scaling', 
                   'tau4cwd', 'tau4l1', 'tau4l2', 'tau4s1', 'tau4s2_death', 'tau4s2_enz', 'tau4s3', 'tau4s4', 
                   'Km_assim', 'Km_decom', 
                   'fcwdl2',
                   'pl1s1', 'pl2s1', 'pl3s4', 'l1_cue', 'l2l3_cue', 
                   'mic_cue', 'pdeath2soc',
                   'b')

# results from mic model data assimilation
mic_mcmc_summary_para = readMat(paste(data_dir_output, 'mcmc_summary_', mic_model_name, '/', mic_model_name, '_para_mean', '.mat', sep = ''))
mic_mcmc_summary_para = mic_mcmc_summary_para$para.mean

mic_mcmc_summary_gr = readMat(paste(data_dir_output, 'mcmc_summary_', mic_model_name, '/', mic_model_name, '_para_gr', '.mat', sep = ''))
mic_mcmc_summary_gr = mic_mcmc_summary_gr$para.gr
mic_mcmc_summary_gr = apply(mic_mcmc_summary_gr, 1, mean, na.rm = TRUE)

mic_mcmc_summary_r2 = readMat(paste(data_dir_output, 'mcmc_summary_', mic_model_name, '/', mic_model_name, '_stat_r2', '.mat', sep = ''))
mic_mcmc_summary_r2 = mic_mcmc_summary_r2$stat.r2
mic_mcmc_summary_r2 = apply(mic_mcmc_summary_r2, 1, max, na.rm = TRUE)


# results from clm5 data assimilation
clm_mcmc_summary_para = readMat(paste(data_dir_output, 'mcmc_summary_', clm_model_name, '/', clm_model_name, '_para_mean', '.mat', sep = ''))
clm_mcmc_summary_para = clm_mcmc_summary_para$para.mean

clm_mcmc_summary_gr = readMat(paste(data_dir_output, 'mcmc_summary_', clm_model_name, '/', clm_model_name, '_para_gr', '.mat', sep = ''))
clm_mcmc_summary_gr = clm_mcmc_summary_gr$para.gr
clm_mcmc_summary_gr = apply(clm_mcmc_summary_gr, 1, mean, na.rm = TRUE)

clm_mcmc_summary_r2 = readMat(paste(data_dir_output, 'mcmc_summary_', clm_model_name, '/', clm_model_name, '_stat_r2', '.mat', sep = ''))
clm_mcmc_summary_r2 = clm_mcmc_summary_r2$stat.r2
clm_mcmc_summary_r2 = apply(clm_mcmc_summary_r2, 1, max, na.rm = TRUE)


eligible_loc = readMat(paste(data_dir_input, 'data4nn/eligible_profile_loc_0_cesm2_clm5_cen_vr_v2_whole_time.mat', sep = ''))
eligible_loc = eligible_loc$eligible.loc.0

clean_loc_mic = which(mic_mcmc_summary_r2 > 0 & mic_mcmc_summary_gr < 1.05)
clean_loc_clm = which(clm_mcmc_summary_r2 > 0 & clm_mcmc_summary_gr < 1.05)

valid_loc_mic = intersect(eligible_loc, clean_loc_mic)
valid_loc_clm = intersect(eligible_loc, clean_loc_clm)

scaled_default_mic = apply(mic_mcmc_summary_para[valid_loc_mic, ], 2, mean, na.rm = TRUE)
#################################################################################
# figure: distribution of mean values after data assimilation
#################################################################################
color_scheme = c('#e41a1c', '#984ea3', '#ff7f00', '#4daf4a', '#377eb8', '#a65628', 'grey3')

#------------------clm5
current_data_clm = c()
ipara = 1
for (ipara in 1:length(para_list_clm)) {
  middle_data = cbind(clm_mcmc_summary_para[valid_loc_clm, ipara], ipara)
  current_data_clm = rbind(current_data_clm, middle_data)
}
current_data_clm = data.frame(current_data_clm)
colnames(current_data_clm) = c('value', 'name')

p_ensemble_mu_clm =
  ggplot() + 
  geom_violin(data = current_data_clm, aes(y = value, x = as.factor(name)), size = 0.5, trim = FALSE, scale = 'width', fill = '#1f78b4') +
  geom_boxplot(data = current_data_clm, aes(y = value, x = as.factor(name)), outlier.colour = NA, size = 0.5, width = 0.2) +
  geom_point(aes(y = scaled_default_clm5, x = seq(1:length(para_list_clm))), shape = 8, size = 2, stroke = 1, color = 'red') + 
  scale_y_continuous(trans = 'identity', limits = c(0, 1)) +
  scale_x_discrete(labels = para_list_clm) +
  coord_flip() + 
  theme_classic() + 
  # add title
  labs(title = 'CLM5', x = '', y = expression(paste('Ensemble posterior ', mu, ' (scaled)', sep = ''))) + 
  # change the legend properties
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25))  +
  theme(legend.justification = c(0, 1), legend.position = c(0, 1), legend.background = element_rect(fill = NA)) + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, size = 35)) + 
  # modify the font size
  # modify the margin
  theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.5, 'inch')) +
  theme(plot.margin = unit(c(0., 0.2, 0.2, 0.2), 'inch')) +
  theme(axis.text=element_text(size = 20, color = 'black'), axis.text.x = element_text(angle = 0, vjust = 0.8), axis.title = element_text(size = 23), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) 

#------------------mic
current_data_mic = c()
ipara = 1
for (ipara in 1:length(para_list_mic)) {
  middle_data = cbind(mic_mcmc_summary_para[valid_loc_mic, ipara], ipara)
  current_data_mic = rbind(current_data_mic, middle_data)
}
current_data_mic = data.frame(current_data_mic)
colnames(current_data_mic) = c('value', 'name')

p_ensemble_mu_mic = 
  ggplot() + 
  geom_violin(data = current_data_mic, aes(y = value, x = as.factor(name)), size = 0.5, trim = FALSE, scale = 'width', fill = '#1f78b4') +
  geom_boxplot(data = current_data_mic, aes(y = value, x = as.factor(name)), outlier.colour = NA, size = 0.5, width = 0.2) +
  geom_point(aes(y = scaled_default_mic, x = seq(1:length(para_list_mic))), shape = 8, size = 2, stroke = 1, color = 'red') + 
  scale_y_continuous(trans = 'identity', limits = c(0, 1)) +
  scale_x_discrete(labels = para_list_mic) +
  coord_flip() + 
  theme_classic() + 
  # add title
  labs(title = 'COMPAS', x = '', y = expression(paste('Ensemble posterior ', mu, ' (scaled)', sep = ''))) + 
  # change the legend properties
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25))  +
  theme(legend.justification = c(0, 1), legend.position = c(0, 1), legend.background = element_rect(fill = NA)) + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, size = 35)) + 
  # modify the font size
  # modify the margin
  theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.5, 'inch')) +
  theme(plot.margin = unit(c(0., 0.2, 0.2, 0.2), 'inch')) +
  theme(axis.text=element_text(size = 20, color = 'black'), axis.text.x = element_text(angle = 0, vjust = 0.8), axis.title = element_text(size = 23), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) 

jpeg(paste('./Converged_SOC/ensemble_posterior_mu.jpeg', sep = ''), width = 15, height = 12, units = 'in', res = 300)
plot_grid(p_ensemble_mu_clm, p_ensemble_mu_mic, NULL,
          labels = c('a', 'b'),
          nrow = 1,
          rel_widths = c(1, 1, 0.1),
          label_size = 50,
          label_x = 0.02, label_y = 1.02,
          label_fontfamily = 'Arial',
          label_fontface = 'bold'
)
dev.off()



#################################################################################
# figure: posterior distribution at site level data assimilation
#################################################################################
profile_envinfo = readMat(paste(data_dir_input, '/wosis_2019_snap_shot/wosis_2019_snapshot_hugelius_mishra_env_info.mat', sep = ''))
profile_envinfo = profile_envinfo$EnvInfo
profile_full_id = profile_envinfo[ , 2]

wosis_soc_info = nc_open(paste(data_dir_input, 'wosis_2019_snap_shot/soc_data_integrate_wosis_2019_snapshot_hugelius_mishra.nc', sep = ''))
wosis_soc_info = ncvar_get(wosis_soc_info, 'data_soc_integrate')
colnames(wosis_soc_info) = c('profile_id', 'date', 'upper_depth', 'lower_depth', 'node_depth', 'soc_weight', 'soc_stock', 'bulk_denstiy', 'is_pedo')

wosis_profile_info = nc_open(paste(data_dir_input, 'wosis_2019_snap_shot/soc_profile_wosis_2019_snapshot_hugelius_mishra.nc', sep = ''))
wosis_profile_info = ncvar_get(wosis_profile_info, 'soc_profile_info')
colnames(wosis_profile_info) = c('profile_id', 'country_id', 'country_name', 'lon', 'lat', 'layer_num', 'date')


valid_profile = which(wosis_profile_info[ , 'layer_num'] > 6 & mic_mcmc_summary_r2 > 0.9 & clm_mcmc_summary_r2 > 0.9)
valid_profile[sample(1:5000, 8, replace = F)]

profile_list = c(32544, 50463, 32306, 49113, 45516, 38773, 51330, 23268) # c(2948, 9417, 9299, 1075, 101033, 471, 28585, 34372, 33977, 27950)

isample = 8
for (isample in 1:length(profile_list)) { 
  iprofile = profile_list[isample]
  iprofile_id = profile_full_id[iprofile]
  
  profile_lon_lat = profile_envinfo[profile_envinfo[, 2] == iprofile_id, c(4, 5)]
  profile_loc = which(wosis_soc_info[ , 'profile_id'] == iprofile_id)
  layer_depth = wosis_soc_info[profile_loc, 'node_depth']
  layer_obs = wosis_soc_info[profile_loc, 'soc_stock']
  
  #----------------clm5
  middle_data_path = paste(data_dir_output, clm_model_name, '_posterior/',
                           clm_model_name, '_result_profile_', (iprofile), '_', (iprofile_id), '.mat', sep = '')
  
  middle_data = readMat(middle_data_path)
  eval(parse(text = paste('middle_data = middle_data$', rownames(summary(middle_data)), sep = '')))
  post_dist = middle_data[[8]]
  valid_simu_num = middle_data[[10]]
  site_stat = middle_data[[9]]
  soc_layer_mod_clm5 = middle_data[[4]]
  soc_layer_std_clm5 = middle_data[[5]]
  
  best_chain = which(site_stat[ , 1] == max(site_stat[ , 1]))
  for (ichain in 1:3) { 
    post_dist[1, , 1:round(valid_simu_num[ichain]/2)] = NA
  }
  
  current_data = c()
  
  for (ipara in 1:21) {
    current_data = rbind(current_data,
                         cbind(post_dist[best_chain, ipara, ], ipara)
    )
  }
  
  current_data = data.frame(current_data)
  colnames(current_data) = c('value', 'name')
  
  p_clm = 
    ggplot() + 
    geom_violin(data = current_data, aes(y = value, x = as.factor(name)), size = 0.5, trim = FALSE, scale = 'width', fill = '#1f78b4') +
    geom_boxplot(data = current_data, aes(y = value, x = as.factor(name)), outlier.colour = NA, size = 0.5, width = 0.2) +
    geom_point(aes(y = scaled_default_clm5, x = seq(1:length(para_list_clm))), shape = 8, size = 1, stroke = 1, color = 'red') + 
    scale_y_continuous(trans = 'identity', limits = c(0, 1)) +
    scale_x_discrete(labels = para_list_clm) +
    scale_color_manual(name = '', values = c('#e41a1c', '#377eb8'), labels = c('MCMC posterior mean', 'SCE point estimate')) +
    theme_classic() + 
    coord_flip() +
    # add title
    labs(title = 'CLM5', x = '', y = expression(paste('Posterior ', theta, ' (scaled)', sep = ''))) + 
    # change the legend properties
    theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25))  +
    theme(legend.justification = c(0, 1), legend.position = c(0, 1), legend.background = element_rect(fill = NA)) + 
    # modify the position of title
    theme(plot.title = element_text(hjust = 0.5, size = 35)) + 
    # modify the font size
    # modify the margin
    theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.5, 'inch')) +
    theme(plot.margin = unit(c(0., 0.2, 0.2, 0.2), 'inch')) +
    theme(axis.text=element_text(size = 20, color = 'black'), axis.text.x = element_text(angle = 0, vjust = 0.8), axis.title = element_text(size = 23), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) 
  
  
  #----------------mic
  middle_data_path = paste(data_dir_output, mic_model_name, '_posterior/',
                           mic_model_name, '_result_profile_', (iprofile), '_', (iprofile_id), '.mat', sep = '')
  
  middle_data = readMat(middle_data_path)
  eval(parse(text = paste('middle_data = middle_data$', rownames(summary(middle_data)), sep = '')))
  post_dist = middle_data[[8]]
  valid_simu_num = middle_data[[10]]
  site_stat = middle_data[[9]]
  
  soc_layer_mod_mic = middle_data[[4]]
  soc_layer_std_mic = middle_data[[5]]
  
  best_chain = which(site_stat[ , 1] == max(site_stat[ , 1]))
  for (ichain in 1:3) { 
    post_dist[1, , 1:round(valid_simu_num[ichain]/2)] = NA
  }
  
  current_data = c()
  
  for (ipara in 1:23) {
    current_data = rbind(current_data,
                         cbind(post_dist[best_chain, ipara, ], ipara)
    )
  }
  
  current_data = data.frame(current_data)
  colnames(current_data) = c('value', 'name')
  
  p_mic = 
    ggplot() + 
    geom_violin(data = current_data, aes(y = value, x = as.factor(name)), size = 0.5, trim = FALSE, scale = 'width', fill = '#1f78b4') +
    geom_boxplot(data = current_data, aes(y = value, x = as.factor(name)), outlier.colour = NA, size = 0.5, width = 0.2) +
    geom_point(aes(y = scaled_default_mic, x = seq(1:length(para_list_mic))), shape = 8, size = 2, stroke = 1, color = 'red') + 
    scale_y_continuous(trans = 'identity', limits = c(0, 1)) +
    scale_x_discrete(labels = para_list_mic) +
    scale_color_manual(name = '', values = c('#e41a1c', '#377eb8'), labels = c('MCMC posterior mean', 'SCE point estimate')) +
    theme_classic() + 
    coord_flip() +
    # add title
    labs(title = 'COMPAS', x = '', y = expression(paste('Posterior ', theta, ' (scaled)', sep = ''))) + 
    # change the legend properties
    theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25))  +
    theme(legend.justification = c(0, 1), legend.position = c(0, 1), legend.background = element_rect(fill = NA)) + 
    # modify the position of title
    theme(plot.title = element_text(hjust = 0.5, size = 35)) + 
    # modify the font size
    # modify the margin
    theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.5, 'inch')) +
    theme(plot.margin = unit(c(0., 0.2, 0.2, 0.2), 'inch')) +
    theme(axis.text=element_text(size = 20, color = 'black'), axis.text.x = element_text(angle = 0, vjust = 0.8), axis.title = element_text(size = 23), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) 
  
  
  current_data = rbind(cbind(layer_depth, soc_layer_mod_mic/1000, soc_layer_std_mic/1000, 3), 
                       cbind(layer_depth, soc_layer_mod_clm5/1000, soc_layer_std_clm5/1000, 2), 
                       cbind(layer_depth, layer_obs/1000, NA, 1)
                       )
  current_data = data.frame(current_data)
  colnames(current_data) = c('depth', 'mod', 'std', 'model')
  
  p_da =
    ggplot() + 
    geom_point(data = current_data, aes(x = mod, y = -depth, color = as.factor(model)), shape = 2, stroke = 2, size = 7) +
    geom_errorbar(data = current_data, aes(xmin = mod - std, xmax = mod + std, y = -depth, color = as.factor(model)), size = 2, width = 2) +
    scale_color_manual(name = '', values = c('black', '#e41a1c', '#377eb8'), labels = c('Observations', 'CLM5', 'COMPAS')) +
    scale_x_continuous(position = 'top') +  
    theme_classic() + 
    # coord_flip() +
    # add title
    labs(title = '', y = 'Soil depth (cm)', x = expression(paste('SOC content (kg C m'^'-3', ')', sep = ''))) + 
    # change the legend properties
    theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25))  +
    theme(legend.justification = c(1, 0), legend.position = c(1, 0), legend.background = element_rect(fill = NA)) + 
    # modify the position of title
    theme(plot.title = element_text(hjust = 0.5, size = 35)) + 
    # modify the font size
    # modify the margin
    theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.5, 'inch')) +
    theme(plot.margin = unit(c(0., 0.2, 0.2, 0.2), 'inch')) +
    theme(axis.text=element_text(size = 20, color = 'black'), axis.text.x = element_text(angle = 0, vjust = 0.8), axis.title = element_text(size = 23), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) 
  
  
  jpeg(paste('./Converged_SOC/posterior_dist_', iprofile, '.jpeg', sep = ''), width = 22, height = 12, units = 'in', res = 300)
  print(plot_grid(p_clm, p_mic, p_da, NULL,
                  labels = c('a', 'b', 'c'),
                  nrow = 1,
                  rel_widths = c(1, 1.1, 1.2, 0.1),
                  label_size = 50,
                  label_x = 0.02, label_y = 1.02,
                  label_fontfamily = 'Arial',
                  label_fontface = 'bold'
  ))
  dev.off()
  
}

