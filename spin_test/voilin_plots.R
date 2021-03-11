library(tidyr)
library(ggplot2)
library(reshape)
#################
### set home directory
homedir <- "/Users/eballer/BBL/imco/pmacs/PMACS_remote"
#homedir <- "/project/imco"
model = "gam_sex"
#model= "lm_exec_accuracy"
network_names <- c("VIS", "MOT", "DA", "VA", "LIM", "FP", "DM")

#set yeo colors manually - this can be obtained through https://surfer.nmr.mgh.harvard.edu/fswiki/CorticalParcellation_Yeo2011
   #I typed the RGB values into a rgb -> hex converter and stored the values here. Works!

yeo_colors <- c(
    `VIS` = "#781286",
    `MOT` = "#4682b4",
    `DA` = "#00760e",
    `VA` = "#c43afa",
    `LIM` = "#dcf8a4",
    `FP` = "#e69422",
    `DM` = "#cd3e56")

lh_spin_df <- data.frame(read.table(paste0(homedir, "/baller/results/coupling_accuracy//spin_test_results/lh_spin_test_", model, "_proportions.csv"), sep = ",")
)

actual_results <- lh_spin_df[,1]
lh_spin_without_target_col <- lh_spin_df[,2:1001]
#melt dataframe so it is in a good format for violin plotting.
#melt df so it is in a good position to be plotted
melted_df_network_name <- rep(x = network_names, each = 1000)
melted_df_network_num <- rep(x = seq(1:7), each = 1000)
melted_df_spin_results <- rbind(t(lh_spin_without_target_col[1,]), 
                                t(lh_spin_without_target_col[2,]),
                                t(lh_spin_without_target_col[3,]),
                                t(lh_spin_without_target_col[4,]),
                                t(lh_spin_without_target_col[5,]),
                                t(lh_spin_without_target_col[6,]),
                                t(lh_spin_without_target_col[7,]))

melted_df <- as.data.frame(cbind(melted_df_network_name, melted_df_network_num,melted_df_spin_results))
names(melted_df) <- c("network_name", "network_num", "spin")
melted_df$spin <- as.numeric(as.character(melted_df$spin))

#with mean lines, different fonts
#save images
jpeg(file=paste0(homedir, "/baller/results/images/spin_", model, "_t_fdr05.jpeg"))
ggplot(melted_df, aes(x = factor(network_name, level = network_names), y = spin, fill = network_name)) +  
  scale_fill_manual(values=yeo_colors) + 
  geom_violin(trim = TRUE) + 
  xlab("Yeo 7 Network") + ylab(paste0("Proportion")) +
  scale_y_continuous(limits=c(-0.1, .75)) + 
  #scale_y_continuous(limits=c(-0.1, .4)) + 
  geom_violin(trim=FALSE) + 
  theme(legend.position = "none",
        legend.title = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 10)) +
  stat_summary(fun.y = mean, geom = "errorbar", 
              aes(ymax = ..y.., ymin = ..y.., group = factor(network_name)),
               width = 0.5, linetype = "dashed", position = position_dodge(0.9)) + 
  #geom_boxplot(width = 0.15, position = position_dodge(0.9)) + 
  geom_segment(aes(x = 0.5, y = actual_results[1], xend = 1.5, yend = actual_results[1])) + 
  geom_segment(aes(x = 1.5, y = actual_results[2], xend = 2.5, yend = actual_results[2])) +
  geom_segment(aes(x = 2.5, y = actual_results[3], xend = 3.5, yend = actual_results[3])) +
  geom_segment(aes(x = 3.5, y = actual_results[4], xend = 4.5, yend = actual_results[4])) +
  geom_segment(aes(x = 4.5, y = actual_results[5], xend = 5.5, yend = actual_results[5])) +
  geom_segment(aes(x = 5.5, y = actual_results[6], xend = 6.5, yend = actual_results[6])) +
  geom_segment(aes(x = 6.5, y = actual_results[7], xend = 7.5, yend = actual_results[7])) +
  ggtitle(paste0("Spin Test Permutation Results ", model))
dev.off()

#try bars
jpeg(file=paste0(homedir,"/baller/results/images/spin_", model, "_t_fdr05_bar_graph_means.jpeg"))
means <- as.vector(rowMeans(lh_spin_without_target_col))
actual_results <- as.vector(actual_results)
spin_results <- rbind(actual_results,means)
melted_for_plot <- melt(spin_results)
melted_for_plot$network_name <- rep(network_names, each = 2)
names(melted_for_plot) <- c("group_name", "group_num", "value", "net")
melted_for_plot$value <- as.numeric(as.character(melted_for_plot$value))


#ggplot(melted_for_plot, aes(fill=group_name, y=value, x=net)) + 
ggplot(data=melted_for_plot, aes(fill=group_name, y=value, x=factor(net, level = network_names))) + 
#ggplot(data=melted_for_plot, aes(fill=yeo_colors, y=value, x=factor(net, level = network_names))) + 

  geom_bar(position="dodge", stat="identity") + 
ggtitle(paste0("Spin Test Permutation Results ", model))
dev.off()
