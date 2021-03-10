library(tidyr)
library(ggplot2)
library(reshape)
#################
### set home directory
homedir <- "/Users/eballer/BBL/imco/pmacs/PMACS_remote"
#homedir <- "/project/imco"
model = "gam_sex"
network_names <- c("VIS", "MOT", "DA", "VA", "LIM", "FP", "DM")

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
jpeg(file=paste0(homedir, "/baller/results/images/spin_gam_sex_t_fdr05.jpeg"))
ggplot(melted_df, aes(x = network_name, y = spin, fill = network_name)) +  
  geom_violin(trim = TRUE) + 
  xlab("Yeo 7 Network") + ylab(paste0("Proportion")) +
  scale_y_continuous(limits=c(-0.1, 1)) + 
  geom_violin(trim=FALSE) + 
  theme(legend.position = "right",
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
  ggtitle(paste0("Spin Test Permutation Results"))
dev.off()

#try bars
jpeg(file=paste0(homedir, "/baller/results/images/spin_gam_sex_t_fdr05_bar_graph_means.jpeg"))
means <- as.vector(rowMeans(lh_spin_without_target_col))
actual_results <- as.vector(actual_results)
spin_results <- rbind(actual_results,means)
melted_for_plot <- melt(spin_results)
melted_for_plot$network_name <- rep(network_names, each = 2)
names(melted_for_plot) <- c("group_name", "group_num", "value", "net")
melted_for_plot$value <- as.numeric(as.character(melted_for_plot$value))


ggplot(melted_for_plot, aes(fill=group_name, y=value, x=net)) + 
  geom_bar(position="dodge", stat="identity")
dev.off()
