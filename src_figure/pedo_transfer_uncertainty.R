library(ncdf4)
library(R.matlab)
library(ggplot2)
# library(qlcMatrix)
# library(jcolors)
library(cowplot)
library(scales)

# dev.off()
rm(list = ls())
setwd('/Users/ft254/Google_Drive/R')

#######################################
# wosis database
#######################################
profile_info_colname = c('profile_id', 'country_id', 'country_name', 'lat', 'lon', 'layer_num', 'date')
soc_data_colname = c('profile_id', 'date', 'upper_depth', 'lower_depth', 'node_depth', 'soc_weight', 'soc_stock', 'bulk_denstiy', 'is_pedo')

ncfname = '/Users/ft254/DATAHUB/ENSEMBLE/INPUT_DATA/wosis_2019_snap_shot/soc_profile_wosis_2019_snapshot_hugelius_mishra.nc'
profile_info = nc_open(ncfname)
profile_info = ncvar_get(profile_info)
colnames(profile_info) = profile_info_colname

ncfname = '/Users/ft254/DATAHUB/ENSEMBLE/INPUT_DATA/wosis_2019_snap_shot/soc_data_integrate_wosis_2019_snapshot_hugelius_mishra.nc'
profile_data = nc_open(ncfname)
profile_data = ncvar_get(profile_data)
colnames(profile_data) = soc_data_colname

fun_bulk_density = function(soc_weight) {
  fun_bulk_density = 0.32 + 1.30*exp(-0.0089*soc_weight*1.724)
  return(fun_bulk_density)
}

pedo_bd = fun_bulk_density(profile_data[ , 'soc_weight'])

validation_data = cbind(pedo_bd[profile_data[ , 'is_pedo'] == 0], 
                        profile_data[profile_data[ , 'is_pedo'] == 0, 'bulk_denstiy'],
                        pedo_bd[profile_data[ , 'is_pedo'] == 0]*profile_data[profile_data[ , 'is_pedo'] == 0, 'soc_weight'],
                        profile_data[profile_data[ , 'is_pedo'] == 0, 'soc_stock']/1000
                        )
colnames(validation_data) = c('bd_pedo', 'bd_obs', 'soc_pedo', 'soc_obs')


pred_bd_diff = validation_data[ , 'bd_pedo'] - validation_data[ , 'bd_obs']
t.test(pred_bd_diff)

pred_soc_diff = validation_data[ , 'soc_pedo'] - validation_data[ , 'soc_obs'] # unit kg C m-3
t.test(pred_soc_diff)



summary(lm(validation_data[ , 'bd_pedo'] ~ validation_data[ , 'bd_obs']))
summary(lm(log10(validation_data[validation_data[ , 'soc_pedo'] > 0, 'soc_pedo']) ~ log10(validation_data[validation_data[ , 'soc_pedo'] > 0, 'soc_obs'])))
summary(lm((validation_data[validation_data[ , 'soc_pedo'] > 0, 'soc_pedo']) ~ (validation_data[validation_data[ , 'soc_pedo'] > 0, 'soc_obs'])))


jpeg(paste('./Converged_SOC/pedo_uncertainty.jpeg', sep = ''), width = 10, height = 10, units = 'in', res = 300)

ggplot() +  
  geom_point(aes(x = pedo_bd[profile_data[ , 'is_pedo'] == 0]*profile_data[profile_data[ , 'is_pedo'] == 0, 'soc_weight'], y = profile_data[profile_data[ , 'is_pedo'] == 0, 'soc_stock']/1000), shape = 16, size = 0.5, alpha = 0.4, color = '#737373') +
  geom_abline(slope = 1, intercept = 0, color = '#99000d', size = 1) + 
  scale_x_continuous(trans = 'log10', limits = c(0.01, 1000), n.breaks = 6, labels = trans_format('log10', math_format(10^.x))) +
  scale_y_continuous(trans = 'log10', limits = c(0.01, 1000), n.breaks = 6, labels = trans_format('log10', math_format(10^.x))) +
  # change the background to black and white
  theme_classic() + 
  # add title
  labs(title = '', x = expression(paste('Observations (kg C m'^'-3', ')', sep = '')), y = expression(paste('Pedo-predictions (kg C m'^'-3', ')', sep = ''))) + 
  # change the legend properties
  theme(legend.text = element_text(size = 30), legend.title = element_text(size = 30))  +
  theme(legend.justification = c(1, 1), legend.position = c(1, 1.05), legend.background = element_rect(fill = NA)) + 
  theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.5, 'inch')) +
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, size = 50)) + 
  # modify the font size
  # modify the margin
  # theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), 'inch')) +
  theme(axis.text=element_text(size = 30, color = 'black'), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1, color = 'black'), axis.ticks.length = unit(0.12, 'inch'))

dev.off()
