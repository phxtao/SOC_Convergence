## Packages
library(R.matlab)
library(maps)
library(ggplot2)
library(R.matlab)
library(cowplot)
library(viridis)
library(cowplot)
library(FactoMineR)
library(factoextra)
library(ggrepel)
library(corrplot)
library(scales)

# dev.off()
##
rm(list = ls())

setwd('/Users/ft254/Google_Drive/R')

## Jet colorbar function
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
diff.colors <- colorRampPalette(c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "white", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"))

#############################################################################
# Data Path
#############################################################################
model_name = 'cesm2_clm5_mic_vr_v22'
nn_exp_name = 'exp_pc_cesm2_23'
time_domain = 'whole_time'

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
# COMPAS-sensitivity component control
#################################################################################
# color_scheme = c('#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7')
color_scheme = c('#e41a1c', '#ff7f00', '#4daf4a', '#377eb8', '#a65628', 'grey3')
# color_scheme = c('#88CCEE', '#CC6677', '#DDCC77', '#117733', '#332288', '#AA4499', '#44AA99', '#999933', '#882255', '#661100', '#6699CC', '#888888')
# color_scheme = c('#7F3C8D', '#11A579', '#3969AC', '#F2B701', '#E73F74', '#80BA5A', '#E68310', '#008695', '#CF1C90', '#f97b72', '#4b4b8f', '#A5AA99')


bootstrap_num = 200
stat_matric = c('coeff_efficiency', 'correlation', 'coeff_determination', 'mse')

component_name = c('H_NN', 'A_A', 'F_I', 'D_K', 'E_V', 'B_Xi', 'C_NPP', 'A_E', 'Ga_Homo')

depth_name = c('C_30cm', 'B_100cm', 'A_200cm', 'D_full_depth')


control_summary = readMat(paste(data_dir_output, 'world_simulation_analyses/compontent_control_summary_bootstrap_GCB_', model_name, '.mat', sep = ''))
control_summary = control_summary$control.summary
control_summary = data.frame(control_summary)
colnames(control_summary) = c('var_soc', 'upper_soc', 'lower_soc', 'var_stat', 'upper_stat', 'lower_stat', 'depth', 'component')

# control_summary_GCB = readMat(paste(data_dir_output, 'world_simulation_analyses/compontent_control_summary_bootstrap_GCB_', model_name, '.mat', sep = ''))
# control_summary_GCB = control_summary_GCB$control.summary
# control_summary_GCB = data.frame(control_summary_GCB)
# colnames(control_summary_GCB) = c('var_soc', 'upper_soc', 'lower_soc', 'var_stat', 'upper_stat', 'lower_stat', 'depth', 'component')
# 
# control_summary[control_summary$component == 2, ] = control_summary_GCB[control_summary_GCB$component == 2, ]
  
for (icomponent in 1:length(component_name)) {
  control_summary[control_summary[ , 'component'] == icomponent, 'component'] = component_name[icomponent]
}

for (idepth in 1:length(depth_name)) {
  control_summary[control_summary[ , 'depth'] == idepth, 'depth'] = depth_name[idepth]
}

### plot figure component only
current_data_range = control_summary[which(control_summary$depth == 'B_100cm' & control_summary$component == 'H_NN'), ]
current_data = control_summary[which(control_summary$depth == 'B_100cm' & control_summary$component != 'Ga_Homo' & control_summary$component != 'H_NN'), ]

abs(current_data_range$var_soc - current_data$var_soc)
mean(abs(current_data_range$var_soc - current_data$var_soc)[1]/abs(current_data_range$var_soc - current_data$var_soc)[2:7])

abs(current_data_range$var_stat - current_data$var_stat)
mean(abs(current_data_range$var_stat - current_data$var_stat)[1]/abs(current_data_range$var_stat - current_data$var_stat)[2:6])


line_label = c('Carbon transfer', 'Environmental modifer', 'Carbon input', 'Baseline decomposition', 'Vertical transport', 'Carbon input allocation')


S_sqrt = function(x){sign(x)*sqrt(abs(x))}
IS_sqrt = function(x){x^2*sign(x)}
S_sqrt_trans = function() trans_new("S_sqrt",S_sqrt,IS_sqrt)

current_data = current_data[current_data$component != 'A_E', ]

p_sensitivity_compas =
  ggplot(data = current_data) +
  geom_hline(yintercept = 0, size = 2, color = 'snow4', linetype = 'dashed', alpha = 1) +
  geom_vline(xintercept = 0, size = 2, color = 'snow4', linetype = 'dashed', alpha = 1) +
  geom_errorbar(aes(x = var_soc, ymin = lower_stat, ymax = upper_stat, color = component), size = 2, width = 0.7, stat = 'identity') +
  geom_errorbarh(aes(y = var_stat, xmin = lower_soc, xmax = upper_soc, color = component), size = 2, height = 0.02, stat = 'identity') +
  geom_point(aes(x = var_soc, y = var_stat, color = component, fill = component), shape = 15, size = 8) +
  scale_color_manual(name = '', labels = line_label, values = color_scheme) +
  scale_fill_manual(name = '', labels = line_label, values = color_scheme) +
  scale_x_continuous(limits = c(-20, NA), breaks = c(0, 90, 400, 1000, 2000), trans = 'S_sqrt') +
  scale_y_continuous(limits = c(NA, NA), breaks = -c(0, 0.01, 0.04, 0.09, 0.16, 0.25, 0.4), trans = 'S_sqrt') +
  theme_classic() +
  # change the legend properties
  # theme(panel.background = element_rect(fill = 'grey98'), plot.background = element_rect(fill = 'grey98'))+
  theme(legend.justification = c(0, 0), legend.position = c(0, 0), legend.background = element_rect(fill = NA), legend.text.align = 0) +
  theme(legend.text = element_text(size = 30), legend.title = element_text(size = 30))  +
  theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.8, 'inch')) +
  # add title
  labs(x = 'Total absolute error of SOC estimates (Pg C)', y = expression(paste('Deviation of modeling efficiency'), sep = '')) +
  # modify the position of title
  # modify the font size
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), 'inch'), plot.background = element_rect(fill = 'transparent', color = NULL)) +
  theme(axis.text = element_text(size = 30, color = 'black'), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1, color = 'black'), axis.ticks.length = unit(0.12, 'inch'))

#################################################################################
# COMPAS-sensitivity curve
#################################################################################
process_list =  c('A', 'K', 'Xi', 'V', 'I', 'E', 'NPP')
process_label = c('A_A', 'D_K', 'B_Xi', 'E_V', 'F_I', 'A_E', 'C_NPP')

process_description = c('Microbial CUE', 'Baseline Decomposition', 'Environmental Impacts', 'Vertical Transport', 'Input Allocation', 'Non-microbial carbon transfer', 'Carbon Input')
climate_list = c('A', 'B', 'C', 'D', 'E_all')

manage_proposal = seq(from = -0.2, to = 0.2, by = 0.04)*100
depth_name = c('C_30cm', 'B_100cm', 'A_200cm', 'D_full_depth')

marginal_change_total = readMat(paste(data_dir_output, 'world_simulation_analyses/marginal_change_summary_bootstrap_', model_name, '.mat', sep = ''))
marginal_change_total = marginal_change_total$marginal.change.total

marginal_change_total_GCB = readMat(paste(data_dir_output, 'world_simulation_analyses/marginal_change_summary_bootstrap_GCB_', model_name, '.mat', sep = ''))
marginal_change_total_GCB = marginal_change_total_GCB$marginal.change.total
####
marginal_change_total_mean = apply(marginal_change_total, c(1, 2, 3), mean, na.rm = TRUE)
marginal_change_total_upper = apply(marginal_change_total, c(1, 2, 3), quantile, probs = 0.16, na.rm = TRUE)
marginal_change_total_lower = apply(marginal_change_total, c(1, 2, 3), quantile, probs = 0.84, na.rm = TRUE)

marginal_change_total_mean_GCB = apply(marginal_change_total_GCB, c(1, 2, 3), mean, na.rm = TRUE)
marginal_change_total_upper_GCB = apply(marginal_change_total_GCB, c(1, 2, 3), quantile, probs = 0.16, na.rm = TRUE)
marginal_change_total_lower_GCB = apply(marginal_change_total_GCB, c(1, 2, 3), quantile, probs = 0.84, na.rm = TRUE)


col_list = c('stock', 'stock_upper', 'stock_lower', 'process', 'process_upper', 'process_lower', 'process_name')
marginal_change_summary = data.frame(array(NA, dim = c(length(process_list)*length(manage_proposal), length(col_list))))
colnames(marginal_change_summary) = col_list

counter = 1

for (iprocess in 1) {

  for (imanage in 1:length(manage_proposal)) {
    marginal_change_summary[counter, 'process_name'] = process_label[iprocess]
    marginal_change_summary[counter, 'stock'] = marginal_change_total_mean_GCB[1, iprocess, imanage]
    marginal_change_summary[counter, 'stock_upper'] = marginal_change_total_upper_GCB[1, iprocess, imanage]
    marginal_change_summary[counter, 'stock_lower'] = marginal_change_total_lower_GCB[1, iprocess, imanage]

    marginal_change_summary[counter, 'process'] = marginal_change_total_mean_GCB[3, iprocess, imanage]
    marginal_change_summary[counter, 'process_upper'] = marginal_change_total_upper_GCB[3, iprocess, imanage]
    marginal_change_summary[counter, 'process_lower'] = marginal_change_total_lower_GCB[3, iprocess, imanage]

    counter = counter + 1
  }
}


for (iprocess in 2:length(process_list)) {
  
  for (imanage in 1:length(manage_proposal)) {
    marginal_change_summary[counter, 'process_name'] = process_label[iprocess]
    marginal_change_summary[counter, 'stock'] = marginal_change_total_mean[1, iprocess, imanage]
    marginal_change_summary[counter, 'stock_upper'] = marginal_change_total_upper[1, iprocess, imanage]
    marginal_change_summary[counter, 'stock_lower'] = marginal_change_total_lower[1, iprocess, imanage]
    
    marginal_change_summary[counter, 'process'] = marginal_change_total_mean[3, iprocess, imanage]
    marginal_change_summary[counter, 'process_upper'] = marginal_change_total_upper[3, iprocess, imanage]
    marginal_change_summary[counter, 'process_lower'] = marginal_change_total_lower[3, iprocess, imanage]
    
    counter = counter + 1
  }
}


#-----------------------------------------------------------------

stock_reference = marginal_change_summary$stock[6]
marginal_change_summary$stock = (marginal_change_summary$stock/stock_reference-1)*100
marginal_change_summary$stock_upper = (marginal_change_summary$stock_upper/stock_reference-1)*100
marginal_change_summary$stock_lower = (marginal_change_summary$stock_lower/stock_reference-1)*100

marginal_change_summary = marginal_change_summary[marginal_change_summary$process_name != 'A_E', ]

p_curve_compas =
  ggplot(data = marginal_change_summary) +
  geom_ribbon(aes(x = process*100, ymin = stock_lower, ymax = stock_upper, fill = process_name), alpha = 0.15) +
  geom_line(aes(x = process*100, y = stock, color = process_name), alpha = 1, size = 2) +
  geom_hline(yintercept = 0, color = 'grey4', size = 2, linetype = 'dashed', alpha = 0.7) +
  # geom_hline(yintercept = marginal_change_summary$stock[6]*(1+0.004)**30, color = 'red', size = 2, linetype = 'dashed', alpha = 1) +
  geom_vline(xintercept = 0, color = 'grey4', size = 2, linetype = 'dashed', alpha = 0.7) +
  scale_x_continuous(trans = 'identity', n.breaks = 7) +
  scale_y_continuous(trans = 'identity', n.breaks = 7) +
  coord_cartesian(xlim = c(-10, 10), ylim = c(-30, 40)) +
  scale_color_manual(name = '', labels = line_label, values = color_scheme) +
  scale_fill_manual(name = '', labels = line_label, values = color_scheme) +
  # change the background to black and white
  theme_classic() +
  # theme(legend.position = 'None') +
  # theme(panel.background = element_rect(fill = 'grey98'), plot.background = element_rect(fill = 'grey98'))+
  theme(legend.justification = 'right', legend.position = 'none', legend.background = element_rect(fill = NA), legend.text.align = 0) +
  theme(legend.text = element_text(size = 30), legend.title = element_text(size = 30))  +
  theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.8, 'inch')) +
  # add title
  labs(y = 'Changes of global SOC stock (%)', x = paste('Proportional changes of components (%)')) +
  # modify the position of title
  # modify the font sizea
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), 'inch'), plot.background = element_rect(fill = 'transparent', color = NULL)) +
  theme(axis.text=element_text(size = 30, color = 'black'), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1, color = 'black'), axis.ticks.length = unit(0.12, 'inch'))


#################################################################################
# CLM5 - importance of mechanisms to soc storage and spatial variation
#################################################################################
model_name = 'cesm2_clm5_cen_vr_v2'
nn_exp_name = 'exp_pc_cesm2_23'
time_domain = 'whole_time'
bootstrap_num = 200

data_dir_output = '/Users/ft254/DATAHUB/ENSEMBLE/OUTPUT_DATA/'
data_dir_input = '/Users/ft254/DATAHUB/ENSEMBLE/INPUT_DATA/'

component_name = c('H_NN', 'A_A', 'F_I', 'D_K', 'E_V', 'B_Xi', 'C_NPP', 'Ga_Homo', 'Gb_Default')

depth_name = c('C_30cm', 'B_100cm', 'A_200cm', 'D_full_depth')
control_summary = readMat(paste(data_dir_output, 'world_simulation_analyses/compontent_control_summary_bootstrap_GCB_', model_name, '.mat', sep = ''))
control_summary = control_summary$control.summary
control_summary = data.frame(control_summary)
colnames(control_summary) = c('var_soc', 'upper_soc', 'lower_soc', 'var_stat', 'upper_stat', 'lower_stat', 'depth', 'component')


for (icomponent in 1:length(component_name)) {
  control_summary[control_summary[ , 'component'] == icomponent, 'component'] = component_name[icomponent]
}

for (idepth in 1:length(depth_name)) {
  control_summary[control_summary[ , 'depth'] == idepth, 'depth'] = depth_name[idepth]
}

### plot figure component only
current_data_range = control_summary[which(control_summary$depth == 'B_100cm' & control_summary$component == 'H_NN'), ]
current_data = control_summary[which(control_summary$depth == 'B_100cm' & control_summary$component != 'Ga_Homo' & control_summary$component != 'H_NN'), ]

abs(current_data_range$var_soc - current_data$var_soc)
mean(abs(current_data_range$var_soc - current_data$var_soc)[1]/abs(current_data_range$var_soc - current_data$var_soc)[2:7])

abs(current_data_range$var_stat - current_data$var_stat)
mean(abs(current_data_range$var_stat - current_data$var_stat)[1]/abs(current_data_range$var_stat - current_data$var_stat)[2:6])


line_label = c('Carbon transfer', 'Environmental modifer', 'Carbon input', 'Baseline decomposition', 'Vertical transport', 'Carbon input allocation')


S_sqrt = function(x){sign(x)*sqrt(abs(x))}
IS_sqrt = function(x){x^2*sign(x)}
S_sqrt_trans = function() trans_new("S_sqrt",S_sqrt,IS_sqrt)

current_data = current_data[current_data$component != 'Gb_Default', ]

p_sensitivity_clm5 =
  ggplot(data = current_data) +
  geom_hline(yintercept = 0, size = 2, color = 'snow4', linetype = 'dashed', alpha = 1) +
  geom_vline(xintercept = 0, size = 2, color = 'snow4', linetype = 'dashed', alpha = 1) +
  geom_errorbar(aes(x = var_soc, ymin = lower_stat, ymax = upper_stat, color = component), size = 2, width = 0., stat = 'identity') +
  geom_errorbarh(aes(y = var_stat, xmin = lower_soc, xmax = upper_soc, color = component), size = 2, height = 0.0, stat = 'identity') +
  geom_point(aes(x = var_soc, y = var_stat, color = component, fill = component), shape = 15, size = 8) +
  scale_color_manual(name = '', labels = line_label, values = color_scheme) +
  scale_fill_manual(name = '', labels = line_label, values = color_scheme) +
  scale_x_continuous(trans = 'identity') +
  scale_y_continuous(trans = 'identity') +
  coord_cartesian(ylim = c(NA, 0.01), xlim = c(NA, NA)) + 
  theme_classic() +
  # change the legend properties
  # theme(panel.background = element_rect(fill = 'grey98'), plot.background = element_rect(fill = 'grey98'))+
  theme(legend.justification = 'right', legend.position = 'none', legend.background = element_rect(fill = NA), legend.text.align = 0) +
  theme(legend.text = element_text(size = 30), legend.title = element_text(size = 30))  +
  theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.8, 'inch')) +
  # add title
  labs(x = 'Total absolute error of SOC estimates (Pg C)', y = expression(paste('Deviation of modeling efficiency'), sep = '')) +
  # modify the position of title
  # modify the font size
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), 'inch'), plot.background = element_rect(fill = 'transparent', color = NULL)) +
  theme(axis.text = element_text(size = 30, color = 'black'), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1, color = 'black'), axis.ticks.length = unit(0.12, 'inch'))


#################################################################################
# CLM5 - sensitivity curve
#################################################################################
process_list =  c('A', 'K', 'Xi', 'V', 'I', 'NPP')
process_label = c('A_A', 'D_K', 'B_Xi', 'E_V', 'F_I', 'C_NPP')

process_description = c('CUE', 'Baseline Decomposition', 'Environmental Impacts', 'Vertical Transport', 'Input Allocation', 'Carbon Input')
climate_list = c('A', 'B', 'C', 'D', 'E_all')

manage_proposal = seq(from = -0.2, to = 0.2, by = 0.04)*100
depth_name = c('C_30cm', 'B_100cm', 'A_200cm', 'D_full_depth')

grid_lon_lat = readMat(paste(data_dir_output, 'world_simulation_analyses/marginal_sensitivity_proda_', process_list[1], '_', model_name, '_', nn_exp_name, '_bootstrap_', as.character(1), '.mat', sep = ''))

grid_lon_lat = grid_lon_lat$var.data.middle
grid_lon_lat = grid_lon_lat[[1]][[1]][ , c(1:2)]


############################# global soc stock
resolution = 0.5
lat_seq = grid_lon_lat[ , 2]
# area of grid 
radius = 6371008.8
length_top = (2*pi*radius*cos(abs(lat_seq+resolution/2)/180*pi)/360)*resolution
length_down = (2*pi*radius*cos(abs(lat_seq-resolution/2)/180*pi)/360)*resolution
height = (pi*radius/180)*resolution
lat_grid_area = (length_top + length_down)*height/2

marginal_change_summary = read.table(paste(data_dir_output, 'world_simulation_analyses/marginal_change_summary_bootstrap.csv', sep = ''), sep = ',', header = TRUE)

line_label = c('Carbon transfer', 'Environmental modifer', 'Carbon input', 'Baseline decomposition', 'Vertical transport', 'Carbon input allocation')


stock_reference = marginal_change_summary$stock[6]
marginal_change_summary$stock = (marginal_change_summary$stock/stock_reference-1)*100
marginal_change_summary$stock_upper = (marginal_change_summary$stock_upper/stock_reference-1)*100
marginal_change_summary$stock_lower = (marginal_change_summary$stock_lower/stock_reference-1)*100

p_curve_clm5 =
  ggplot(data = marginal_change_summary) +
  geom_ribbon(aes(x = process*100, ymin = stock_lower, ymax = stock_upper, fill = process_name), alpha = 0.15) +
  geom_line(aes(x = process*100, y = stock, color = process_name), alpha = 1, size = 2) +
  geom_hline(yintercept = 0, color = 'grey4', size = 2, linetype = 'dashed', alpha = 0.7) +
  # geom_hline(yintercept = marginal_change_summary$stock[6]*(1+0.004)**30, color = 'red', size = 2, linetype = 'dashed', alpha = 1) +
  geom_vline(xintercept = 0, color = 'grey4', size = 2, linetype = 'dashed', alpha = 0.7) +
  scale_x_continuous(trans = 'identity', n.breaks = 7) +
  scale_y_continuous(trans = 'identity', n.breaks = 7) +
  coord_cartesian(xlim = c(-10, 10), ylim = c(-30, 40)) +
  # annotate('text', x = -10, y = marginal_change_summary$stock[6]*(1+0.004)**30, label = '"4per1000" by 2050 target', hjust = 0, vjust = -0.5, size = 10, color = 'red') +
  scale_color_manual(name = '', labels = line_label, values = color_scheme) +
  scale_fill_manual(name = '', labels = line_label, values = color_scheme) +
  # change the background to black and white
  theme_classic() +
  # theme(legend.position = 'None') +
  # theme(panel.background = element_rect(fill = 'grey98'), plot.background = element_rect(fill = 'grey98'))+
  theme(legend.justification = 'right', legend.position = 'none', legend.background = element_rect(fill = NA), legend.text.align = 0) +
  theme(legend.text = element_text(size = 30), legend.title = element_text(size = 30))  +
  theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.8, 'inch')) +
  # add title
  labs(y = 'Changes of global SOC stock (%)', x = paste('Proportional changes of components (%)')) +
  # modify the position of title
  # modify the font sizea
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), 'inch'), plot.background = element_rect(fill = 'transparent', color = NULL)) +
  theme(axis.text=element_text(size = 30, color = 'black'), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1, color = 'black'), axis.ticks.length = unit(0.12, 'inch'))


jpeg(paste('./Converged_SOC/process_relative_importance.jpeg', sep = ''), width = 22, height = 22, units = 'in', res = 300)
plot_grid(p_sensitivity_clm5, p_curve_clm5,
          p_sensitivity_compas, p_curve_compas,
          nrow = 2, ncol = 2,
          rel_widths = c(1, 1),
          labels = c('a', 'b', 'c', 'd'),
          label_size = 50,
          label_x = 0.0, label_y = 1.02,
          label_fontfamily = 'Arial',
          label_fontface = 'bold'
)

dev.off()
