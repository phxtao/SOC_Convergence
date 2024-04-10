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

para_name = c('bio', 'cryo', 'q10', 'w_scaling', 
              'tau4cwd', 'tau4l1', 'tau4l2', 'tau4s1', 'tau4s2_death', 'tau4s2_enz', 'tau4s3', 'tau4s4', 
              'mm_const_assim', 'mm_const_decom', 
              'fcwdl2', 
              'pl1s1', 'pl2s1', 'pl3s4', 'l1_cue', 'l2l3_cue', 
              'mic_cue', 'pdeath2soc', 
              'beta')

n_para = length(para_name)
n_site = 5000

current_data = c()

isite = 1
for (isite in 1:n_site) {
  site_results_dir = paste(data_dir_output, 'cesm2_clm5_mic_vr_v22_recover/cesm2_clm5_mic_vr_v22_result_profile_', isite, '.mat', sep = '')
  
  if (file.exists(site_results_dir) == TRUE) {
    mcmc_site = readMat(site_results_dir)
    mcmc_site = mcmc_site[[1]]
    
    current_data = rbind(current_data, cbind(t(mcmc_site[[1]]), 
                                             mcmc_site[[3]], 
                                             mcmc_site[[4]], 
                                             c(1:23), 
                                             isite))
    
  }
}

colnames(current_data) = c('origin', 'recover', 'std', 'para', 'site')
current_data = data.frame(current_data)


ipara = 1
for (ipara in 1:n_para) {
  p =
    ggplot() + 
    geom_abline(slope = 1, intercept = 0, size = 1) +
    geom_point(data = current_data[current_data$para == ipara, ], aes(x = origin, y = recover), color = 'black', alpha = 0.7, shape = 16, size = 4) + 
    scale_x_continuous(limits = c(0.0, 1), trans = 'identity') +
    scale_y_continuous(limits = c(0.0, 1), trans = 'identity') +
    theme_classic() + 
    # add title
    labs(title = para_name[ipara], x = expression(paste('Origin', sep = '')), y = expression(paste('Recover', sep = ''))) + 
    # change the legend properties
    guides(fill = guide_colorbar(direction = 'horizontal', barwidth = 15, barheight = 2.5, title.position = 'left', title.hjust = 0, title.vjust = 0.3, label.position = 'top', label.hjust = 0.5, frame.linewidth = 0), reverse = FALSE) +
    theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25))  +
    theme(legend.justification = c(1, 0), legend.position = 'none', legend.background = element_rect(fill = NA)) + 
    # modify the position of title
    theme(plot.title = element_text(hjust = 0.5, size = 30)) + 
    # modify the font size
    # modify the margin
    # theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
    theme(plot.margin = unit(c(0., 0.2, 0.2, 0.2), 'inch')) +
    theme(axis.text=element_text(size = 20, color = 'black'), axis.title = element_text(size = 20), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1, color = 'black'), axis.ticks.length = unit(0.12, 'inch')) 
  
  eval(parse(text = paste('p', ipara, ' = p', sep = '')))
  
}



jpeg(paste('./Converged_SOC/para_recovery.jpeg', sep = ''), width = 30, height = 20, units = 'in', res = 300)
plot_grid(p1, p2, p3, p4, p5, p6, 
          p7, p8, p9, p10, p11, p12, 
          p13, p14, p15, p16, p17, p18, 
          p19, p20, p21, p22, p23, NULL, 
          nrow = 4, ncol = 6 ,
          rel_widths = c(1, 1, 1, 1, 1, 1),
          label_size = 70,
          label_x = 0.05, label_y = 1.05,
          label_fontfamily = 'Arial',
          label_fontface = 'bold'
)
dev.off()


