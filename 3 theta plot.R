library(lmerTest)
library(rstatix)
library(dplyr)
library(bruceR)
library(gridExtra)


set.wd()
source('summarySE.R')

stat_data<-read.csv("3 theta_stat/mean_theta_params.csv",header = TRUE)
stat_data <- stat_data %>% filter(!is.na(CF))
# View(stat_data)
stat_data$age<-as.factor(stat_data$age)


############
# Age group
############
CF_sub_1 <- stat_data %>% summarySE(measurevar = "CF", groupvars=c('PID', "age")) 
CF_sub_1 <- CF_sub_1 %>% 
  rename(y_mean = CF_mean, y_median = CF_median) %>% 
  mutate(P_para = "CF")
PW_sub_1 <- stat_data %>% summarySE(measurevar = "PW", groupvars=c('PID', "age")) 
PW_sub_1 <- PW_sub_1 %>% 
  rename(y_mean = PW_mean, y_median = PW_median) %>% 
  mutate(P_para = "PW")
sub_p <- bind_rows(CF_sub_1, PW_sub_1)

CF_1 <- summarySE(CF_sub_1, measurevar = "y_mean", groupvars=c("age"))%>% 
  mutate(P_para = "CF")
PW_1 <- summarySE(PW_sub_1, measurevar = "y_mean", groupvars=c("age"))%>% 
  mutate(P_para = "PW")
mean_p <- bind_rows(CF_1, PW_1)

CF_data <- mean_p %>% filter(P_para == "CF")
CF_data_sub <- sub_p %>% filter(P_para == "CF")

CF_data <- CF_data %>% mutate(P_para = paste(P_para, age, sep = "_"))
CF_data_sub <- CF_data_sub %>% mutate(P_para = paste(P_para, age, sep = "_"))

p1 <- ggplot(CF_data, aes(x = P_para, y = y_mean_mean)) +
  geom_errorbar(aes(x = P_para, y = y_mean_mean,
                    ymin = y_mean_mean-se, ymax = y_mean_mean+se), width = .2, position = position_nudge(x = -0.15)) + 
  geom_point(aes(x = P_para, y = y_mean_mean), color="#551C6E", fill="#551C6E",
             position = position_nudge(x = -0.15), shape=21, size=3) +
  geom_jitter(data = CF_data_sub, aes(x = P_para, y = y_mean), color="#551C6E", fill="#551C6E",
              position = position_jitter(width = 0), shape=21, size = 2, alpha=0.4, show.legend=F) +
  # labs(title = 'CF by age') +
  ylab('Central frequency (Hz)')+
  scale_y_continuous(breaks = seq(0, max(CF_data_sub$y_mean))) +
  xlab('Age')+
  scale_x_discrete(labels = c("CF_0" = "Children", "CF_1" = "Adults")) + # Add x-axis labels
  theme_bw() + # Change theme to black and white
  theme(
  panel.grid.major.x = element_blank(),
  panel.grid.minor = element_blank(),
  legend.title = element_blank(),
  legend.position = "none",
  strip.text = element_text(size = 20), # Increase size for strip text
  legend.text = element_text(size = 20), # Increase size for legend text
  axis.text.x = element_text(size = 20), # Increase size for axis text
  plot.title = element_text(size = 24, hjust = 0.5), # Increase size for plot title
  axis.title = element_text(size = 20) # Increase size for axis title
)
p1

# For the power plot
CF_data <- mean_p %>% filter(P_para == "PW")
CF_data_sub <- sub_p %>% filter(P_para == "PW")

# Modify the data to include 'age' in 'P_para'
CF_data <- CF_data %>% mutate(P_para = paste(P_para, age, sep = "_"))
CF_data_sub <- CF_data_sub %>% mutate(P_para = paste(P_para, age, sep = "_"))

p2 <- ggplot(CF_data, aes(x = P_para, y = y_mean_mean)) +
  geom_errorbar(aes(x = P_para, y = y_mean_mean, ymin = y_mean_mean-se,
                    ymax = y_mean_mean+se), width = .2, position = position_nudge(x = -0.15)) + 
  geom_point(aes(x = P_para, y = y_mean_mean), color="#551C6E", fill="#551C6E",
             position = position_nudge(x = -0.15), shape=21, size=3) +
  geom_jitter(data = CF_data_sub, aes(x = P_para, y = y_mean),color="#551C6E", fill="#551C6E",
              position = position_jitter(width = 0), shape=21, size = 2, alpha=0.4, show.legend=F) +
  #   labs(title = 'PW by age') +
  ylab(expression('Peak power ' * '(' * mu * 'V'^'2' * ')'))+
  scale_y_continuous(breaks = seq(0, max(CF_data_sub$y_mean))) +
  xlab('Age')+
  # Note p significance manually
  geom_text(data = CF_data, x = 1.5, y = max(CF_data_sub$y_mean),
           label = "**", size = 8) + # p = 0.003
  scale_x_discrete(labels = c("PW_0" = "Children", "PW_1" = "Adults")) + # Add x-axis labels
  theme_bw() + # Change theme to black and white
  theme(
  panel.grid.major.x = element_blank(),
  panel.grid.minor = element_blank(),
  legend.title = element_blank(),
  legend.position = "none",
  strip.text = element_text(size = 20), # Increase size for strip text
  legend.text = element_text(size = 20), # Increase size for legend text
  axis.text.x = element_text(size = 20), # Increase size for axis text
  plot.title = element_text(size = 24, hjust = 0.5), # Increase size for plot title
  axis.title = element_text(size = 20) # Increase size for axis title
)
p2



#############
## load group
#############
CF_sub_1 <- stat_data %>% summarySE(measurevar = "CF", groupvars=c('PID', "load")) 
CF_sub_1 <- CF_sub_1 %>% 
  rename(y_mean = CF_mean, y_median = CF_median) %>% 
  mutate(P_para = "CF")
PW_sub_1 <- stat_data %>% summarySE(measurevar = "PW", groupvars=c('PID', "load")) 
PW_sub_1 <- PW_sub_1 %>% 
  rename(y_mean = PW_mean, y_median = PW_median) %>% 
  mutate(P_para = "PW")
sub_p <- bind_rows(CF_sub_1, PW_sub_1)

CF_1 <- summarySE(CF_sub_1, measurevar = "y_mean", groupvars=c("load"))%>% 
  mutate(P_para = "CF")
PW_1 <- summarySE(PW_sub_1, measurevar = "y_mean", groupvars=c("load"))%>% 
  mutate(P_para = "PW")
mean_p <- bind_rows(CF_1, PW_1)

# CF by load
para_data <- mean_p %>% filter(P_para == "CF")
para_data_sub <- sub_p %>% filter(P_para == "CF")

# Modify the data to include 'age' in 'P_para'
para_data <- para_data %>% mutate(P_para = paste(P_para, load, sep = "_"))
para_data_sub <- para_data_sub %>% mutate(P_para = paste(P_para, load, sep = "_"))

p3 <- ggplot(para_data, aes(x = P_para, y = y_mean_mean)) +
  geom_errorbar(aes(x = P_para, y = y_mean_mean, ymin = y_mean_mean-se,
                    ymax = y_mean_mean+se), width = .2, position = position_nudge(x = -0.15)) + 
  geom_point(aes(x = P_para, y = y_mean_mean, group = load, fill=load),
             position = position_nudge(x = -0.15), shape=21, size=3) +
  geom_jitter(data = para_data_sub, aes(x = P_para, y = y_mean, fill = load),
              position = position_jitter(width = 0), shape=21, size = 2, alpha=0.4, show.legend=F) +
  scale_fill_manual(values = c("#93B237", "#551C6E"),
                    labels = c("1b" = "1-back", "2b" = "2-back"))+
  labs(title = 'CF by load') +
  ylab('Central frequency (Hz)')+
  scale_y_continuous(breaks = seq(0, max(para_data_sub$y_mean))) +
  xlab('Load')+
  scale_x_discrete(labels = c("CF_1b" = "1-back", "CF_2b" = "2-back")) + # Add x-axis labels
  theme_bw() + # Change theme to black and white
  theme(
  panel.grid.major.x = element_blank(),
  panel.grid.minor = element_blank(),
  legend.title = element_blank(),
  legend.position = "none",
  strip.text = element_text(size = 20), # Increase size for strip text
  legend.text = element_text(size = 20), # Increase size for legend text
  axis.text.x = element_text(size = 20), # Increase size for axis text
  plot.title = element_text(size = 24, hjust = 0.5), # Increase size for plot title
  axis.title = element_text(size = 20) # Increase size for axis title
)
p3

# power plot
para_data <- mean_p %>% filter(P_para == "PW")
para_data_sub <- sub_p %>% filter(P_para == "PW")

# Modify the data to include 'age' in 'P_para'
para_data <- para_data %>% mutate(P_para = paste(P_para, load, sep = "_"))
para_data_sub <- para_data_sub %>% mutate(P_para = paste(P_para, load, sep = "_"))

p4 <- ggplot(para_data, aes(x = P_para, y = y_mean_mean)) +
  geom_errorbar(aes(x = P_para, y = y_mean_mean, ymin = y_mean_mean-se,
                    ymax = y_mean_mean+se), width = .2, position = position_nudge(x = -0.15)) + 
  geom_point(aes(x = P_para, y = y_mean_mean), color="#551C6E", fill="#551C6E",
             position = position_nudge(x = -0.15), shape=21, size=3) +
  geom_jitter(data = para_data_sub, aes(x = P_para, y = y_mean), color="#551C6E", fill="#551C6E",
              position = position_jitter(width = 0), shape=21, size = 2, alpha=0.4, show.legend=F) +
  # labs(title = 'PW by load') +
  ylab(expression('Peak power ' * '(' * mu * 'V'^'2' * ')'))+
  scale_y_continuous(breaks = seq(0, max(para_data_sub$y_mean))) +
  xlab('Load')+
  # Note p significance manually
  geom_text(data = para_data, x = 1.5, y = max(para_data_sub$y_mean),
           label = "*", size = 8) + # p = 0.016
  scale_x_discrete(labels = c("PW_1b" = "1-back", "PW_2b" = "2-back")) + # Add x-axis labels
  theme_bw() + # Change theme to black and white
  theme(
  panel.grid.major.x = element_blank(),
  panel.grid.minor = element_blank(),
  legend.title = element_blank(),
  legend.position = "none",
  strip.text = element_text(size = 20), # Increase size for strip text
  legend.text = element_text(size = 20), # Increase size for legend text
  axis.text.x = element_text(size = 20), # Increase size for axis text
  plot.title = element_text(size = 24, hjust = 0.5), # Increase size for plot title
  axis.title = element_text(size = 20) # Increase size for axis title
)
p4


png("4 theta pic/Fig2 b power.png", width = 8, height = 8, units = "in", res = 350)
grid.arrange(p2, p4, ncol = 2)
dev.off()


##############
# exp & power
##############
CF_sub_1 <- stat_data %>% summarySE(measurevar = "CF", groupvars=c('PID', "exp"))
CF_sub_1 <- CF_sub_1 %>%
  rename(y_mean = CF_mean, y_median = CF_median) %>%
  mutate(P_para = "CF")
PW_sub_1 <- stat_data %>% summarySE(measurevar = "PW", groupvars=c('PID', "exp"))
PW_sub_1 <- PW_sub_1 %>%
  rename(y_mean = PW_mean, y_median = PW_median) %>%
  mutate(P_para = "PW")
sub_p <- bind_rows(CF_sub_1, PW_sub_1)

CF_1 <- summarySE(CF_sub_1, measurevar = "y_mean", groupvars=c("exp"))%>%
  mutate(P_para = "CF")
PW_1 <- summarySE(PW_sub_1, measurevar = "y_mean", groupvars=c("exp"))%>%
  mutate(P_para = "PW")
mean_p <- bind_rows(CF_1, PW_1)


para_data <- mean_p %>% filter(P_para == "PW")
para_data_sub <- sub_p %>% filter(P_para == "PW")

# Modify the data to include 'age' in 'P_para'
para_data <- para_data %>% mutate(P_para = paste(P_para, exp, sep = "_"))
para_data_sub <- para_data_sub %>% mutate(P_para = paste(P_para, exp, sep = "_"))

p6 <- ggplot(para_data, aes(x = P_para, y = y_mean_mean)) +
  geom_errorbar(aes(x = P_para, y = y_mean_mean, ymin = y_mean_mean-se,
                    ymax = y_mean_mean+se), width = .2, position = position_nudge(x = -0.15)) +
  geom_point(aes(x = P_para, y = y_mean_mean), color="#551C6E", fill="#551C6E", 
             position = position_nudge(x = -0.15), shape=21, size=3) +
  geom_jitter(data = para_data_sub, aes(x = P_para, y = y_mean), color="#551C6E", fill="#551C6E",
              position = position_jitter(width = 0), shape=21, size = 2, alpha=0.4, show.legend=F) +
  # scale_fill_manual(values = c("#93B237", "#551C6E"),
  #                   labels = c("wmc" = "Chinese character", "wmv" = "visual pattern"))+
  # labs(title = 'PW by exp') +
  ylab(expression('Peak power ' * '(' * mu * 'V'^'2' * ')'))+
  scale_y_continuous(breaks = seq(0, max(para_data_sub$y_mean))) +
  xlab('Stimulus type')+
  # 分条件给图注释p值
  geom_text(data = para_data, x = 1.5, y = max(para_data_sub$y_mean),
            label = "*", size = 8) + # p = 0.043
  scale_x_discrete(labels = c("PW_wmc" = "Chinese character", "PW_wmv" = "Visual pattern")) + # Add x-axis labels
  theme_bw() + # Change theme to black and white
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.position = "none",
    strip.text = element_text(size = 20), # Increase size for strip text
    legend.text = element_text(size = 20), # Increase size for legend text
    axis.text.x = element_text(size = 20), # Increase size for axis text
    plot.title = element_text(size = 24, hjust = 0.5), # Increase size for plot title
    axis.title = element_text(size = 20) # Increase size for axis title
  )
p6
ggsave('4 theta pic/Fig2 d PW exp.png', width = 6, height = 8, plot = p6, dpi=350)


##############
# exp & offset
##############
CF_sub_1 <- stat_data %>% summarySE(measurevar = "exponent", groupvars=c('PID', "age")) 
CF_sub_1 <- CF_sub_1 %>% 
  rename(y_mean = exponent_mean, y_median = exponent_median) %>% 
  mutate(P_para = "exponent")
PW_sub_1 <- stat_data %>% summarySE(measurevar = "offset", groupvars=c('PID', "age")) 
PW_sub_1 <- PW_sub_1 %>% 
  rename(y_mean = offset_mean, y_median = offset_median) %>% 
  mutate(P_para = "offset")
sub_p <- bind_rows(CF_sub_1, PW_sub_1)

CF_1 <- summarySE(CF_sub_1, measurevar = "y_mean", groupvars=c("age"))%>% 
  mutate(P_para = "exponent")
PW_1 <- summarySE(PW_sub_1, measurevar = "y_mean", groupvars=c("age"))%>% 
  mutate(P_para = "offset")
mean_p <- bind_rows(CF_1, PW_1)

para_data <- mean_p %>% filter(P_para == "exponent")
para_data_sub <- sub_p %>% filter(P_para == "exponent")

# Modify the data to include 'age' in 'P_para'
para_data <- para_data %>% mutate(P_para = paste(P_para, age, sep = "_"))
para_data_sub <- para_data_sub %>% mutate(P_para = paste(P_para, age, sep = "_"))

p7 <- ggplot(para_data, aes(x = P_para, y = y_mean_mean)) +
  geom_errorbar(aes(x = P_para, y = y_mean_mean, ymin = y_mean_mean-se,
                    ymax = y_mean_mean+se), width = .2, position = position_nudge(x = -0.15)) + 
  geom_point(aes(x = P_para, y = y_mean_mean), color="#551C6E", fill="#551C6E",
             position = position_nudge(x = -0.15), shape=21, size=3) +
  geom_jitter(data = para_data_sub, aes(x = P_para, y = y_mean), color="#551C6E", fill="#551C6E",
              position = position_jitter(width = 0), shape=21, size = 2, alpha=0.4, show.legend=F) +
  # labs(title = 'exponent by age') +
  ylab(expression('Exponent '*'('*mu*'V'^'2'*'Hz'['-1']*')'))+
  scale_y_continuous(breaks = seq(0, max(para_data_sub$y_mean))) +
  xlab('Age')+
  # 分条件给图注释p值
  geom_text(data = para_data, x = 1.5, y = max(para_data_sub$y_mean),
            label = "***", size = 8) + # p < 0.001
  scale_x_discrete(labels = c("exponent_0" = "Children", "exponent_1" = "Adults")) + # Add x-axis labels
  theme_bw() + # Change theme to black and white
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.position = "none",
    strip.text = element_text(size = 20), # Increase size for strip text
    legend.text = element_text(size = 20), # Increase size for legend text
    axis.text.x = element_text(size = 20), # Increase size for axis text
    plot.title = element_text(size = 24, hjust = 0.5), # Increase size for plot title
    axis.title = element_text(size = 20) # Increase size for axis title
  )
p7

para_data <- mean_p %>% filter(P_para == "offset")
para_data_sub <- sub_p %>% filter(P_para == "offset")

# Modify the data to include 'age' in 'P_para'
para_data <- para_data %>% mutate(P_para = paste(P_para, age, sep = "_"))
para_data_sub <- para_data_sub %>% mutate(P_para = paste(P_para, age, sep = "_"))

p8 <- ggplot(para_data, aes(x = P_para, y = y_mean_mean)) +
  geom_errorbar(aes(x = P_para, y = y_mean_mean, ymin = y_mean_mean-se,
                    ymax = y_mean_mean+se), width = .2, position = position_nudge(x = -0.15)) + 
  geom_point(aes(x = P_para, y = y_mean_mean), color="#551C6E", fill="#551C6E",
             position = position_nudge(x = -0.15), shape=21, size=3) +
  geom_jitter(data = para_data_sub, aes(x = P_para, y = y_mean), color="#551C6E", fill="#551C6E",
              position = position_jitter(width = 0), shape=21, size = 2, alpha=0.4, show.legend=F) +
  # scale_fill_manual(values = c("#93B237", "#551C6E"),
  #                   labels = c("0" = "Child", "1" = "Adult"))+
  # labs(title = 'offset by age') +
  ylab(expression('Offset ' * '(' * mu * 'V'^'2' * ')'))+
  scale_y_continuous(breaks = seq(0, max(para_data_sub$y_mean))) +
  xlab('Age')+
  # 分条件给图注释p值
  geom_text(data = para_data, x = 1.5, y = max(para_data_sub$y_mean),
            label = "***", size = 8) + # p = 0.001
  scale_x_discrete(labels = c("offset_0" = "Children", "offset_1" = "Adults")) + # Add x-axis labels
  theme_bw() + # Change theme to black and white
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.position = "none",
    strip.text = element_text(size = 20), # Increase size for strip text
    legend.text = element_text(size = 20), # Increase size for legend text
    axis.text.x = element_text(size = 20), # Increase size for axis text
    plot.title = element_text(size = 24, hjust = 0.5), # Increase size for plot title
    axis.title = element_text(size = 20) # Increase size for axis title
  )
p8

png("4 theta pic/Fig2 f ap.png", width = 8, height = 8, units = "in", res = 350)
grid.arrange(p7, p8, ncol = 2)
dev.off()

##############
# exp by load
##############

# exponent by load
CF_sub_1 <- stat_data %>% summarySE(measurevar = "PW", groupvars=c('PID', "load", "age"))
CF_sub_1 <- CF_sub_1 %>%
  rename(y_mean = PW_mean, y_median = PW_median) %>%
  mutate(P_para = "PW")

CF_1 <- summarySE(CF_sub_1, measurevar = "y_mean", groupvars=c("load", 'age'))%>%
  mutate(P_para = "PW")

para_data <- CF_1 #%>% filter(P_para == "PW")
para_data_sub <- CF_sub_1 #%>% filter(P_para == "PW")


p9 <- ggplot(para_data, aes(x = age, y = y_mean_mean, fill = load)) +
  geom_errorbar(aes(x = age, y = y_mean_mean,
                    group = load, ymin = y_mean_mean-se,
                    ymax = y_mean_mean+se), width = .2, position = position_dodge(width = 1.2)) +
  geom_point(aes(x = age, y = y_mean_mean, group = load, fill=load),
             position = position_dodge(width = 1.2), shape=21, size=3) +
  geom_jitter(data = para_data_sub, aes(x = age, y = y_mean, fill = load),
              position = position_jitterdodge(jitter.width = 0, dodge.width = 0.7), shape=21, size = 2, alpha=0.4, show.legend=F) +
  scale_fill_manual(values = c("#93B237", "#551C6E"),
                    labels = c("1b" = "1-back", "2b" = "2-back"))+
  # labs(title = 'PW by load and age') +
  ylab(expression('Peak power ' * '(' * mu * 'V'^'2' * ')'))+
  scale_x_discrete(labels = c("0" = "Child", "1" = "Adult")) +
  scale_y_continuous(breaks = seq(0, max(para_data_sub$y_mean))) +
  xlab('')+
  geom_text(data = para_data, x = 2, y = max(para_data_sub$y_mean),
           label = "***", size = 8) + # p < 0.001
  theme_bw() + # Change theme to black and white
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.position = "bottom",
    strip.text = element_text(size = 20), # Increase size for strip text
    legend.text = element_text(size = 20), # Increase size for legend text
    axis.text.x = element_text(size = 20), # Increase size for axis text
    plot.title = element_text(size = 24, hjust = 0.5), # Increase size for plot title
    axis.title = element_text(size = 20) # Increase size for axis title
  )
p9
ggsave('4 theta pic/Fig2 e exponent load age.png',width = 5, height = 8, plot = p9, dpi=350)





# a<-stat_data%>%filter(age==1&exp=="wmv"&load=="1b")%>%select(CF)
# b<-stat_data%>%filter(age==1&exp=="wmv"&load=="2b")%>%select(CF)
# c<-as.data.frame(cbind(a,b)[complete.cases(cbind(a,b)),])
# s[2]<-t.test(as.vector(c[,1]),as.vector(c[,2]),paired = TRUE)$p.value

# a<-stat_data%>%filter(age==0&exp=="wmc"&load=="1b")%>%select(CF)
# b<-stat_data%>%filter(age==0&exp=="wmc"&load=="2b")%>%select(CF)
# c<-as.data.frame(cbind(a,b)[complete.cases(cbind(a,b)),])
# s[3]<-t.test(as.vector(c[,1]),as.vector(c[,2]),paired = TRUE)$p.value

# a<-stat_data%>%filter(age==1&exp=="wmc"&load=="1b")%>%select(CF)
# b<-stat_data%>%filter(age==1&exp=="wmc"&load=="2b")%>%select(CF)
# c<-as.data.frame(cbind(a,b)[complete.cases(cbind(a,b)),])
# s[4]<-t.test(as.vector(c[,1]),as.vector(c[,2]),paired = TRUE)$p.value
# p.adjust(s, "fdr")
# ####follow up
# res.aov <- anova_test(
#   data=stat_data, dv = offset, wid = PID,between = age,
#   within = c(exp, load)
# )
# get_anova_table(res.aov)
# a<-stat_data%>%filter(age==1&exp=="wmc"&load=="1b")%>%select(offset)
# b<-stat_data%>%filter(age==1&exp=="wmc"&load=="2b")%>%select(offset)
# t.test(a$offset,b$offset,paired = TRUE)
# "ANOVA Table (type III tests)

#         Effect DFn DFd       F        p p<.05      ges
# 1          age   1  54 353.482 2.32e-25     * 8.32e-01
# 2          exp   1  54   0.061 8.06e-01       1.61e-04
# 3         load   1  54   0.647 4.25e-01       8.47e-04
# 4      age:exp   1  54   0.018 8.95e-01       4.69e-05
# 5     age:load   1  54   0.432 5.14e-01       5.65e-04
# 6     exp:load   1  54   0.802 3.75e-01       4.25e-04
# 7 age:exp:load   1  54   1.160 2.86e-01       6.15e-04"

# s<-rep(1:4)
# a<-stat_data%>%filter(age==0&exp=="wmv"&load=="1b")%>%select(offset)
# b<-stat_data%>%filter(age==0&exp=="wmv"&load=="2b")%>%select(offset)
# c<-as.data.frame(cbind(a,b)[complete.cases(cbind(a,b)),])
# s[1]<-t.test(as.vector(c[,1]),as.vector(c[,2]),paired = TRUE)$p.value


# a<-stat_data%>%filter(age==1&exp=="wmv"&load=="1b")%>%select(offset)
# b<-stat_data%>%filter(age==1&exp=="wmv"&load=="2b")%>%select(offset)
# c<-as.data.frame(cbind(a,b)[complete.cases(cbind(a,b)),])
# s[2]<-t.test(as.vector(c[,1]),as.vector(c[,2]),paired = TRUE)$p.value

# a<-stat_data%>%filter(age==0&exp=="wmc"&load=="1b")%>%select(offset)
# b<-stat_data%>%filter(age==0&exp=="wmc"&load=="2b")%>%select(offset)
# c<-as.data.frame(cbind(a,b)[complete.cases(cbind(a,b)),])
# s[3]<-t.test(as.vector(c[,1]),as.vector(c[,2]),paired = TRUE)$p.value

# a<-stat_data%>%filter(age==1&exp=="wmc"&load=="1b")%>%select(offset)
# b<-stat_data%>%filter(age==1&exp=="wmc"&load=="2b")%>%select(offset)
# c<-as.data.frame(cbind(a,b)[complete.cases(cbind(a,b)),])
# s[4]<-t.test(as.vector(c[,1]),as.vector(c[,2]),paired = TRUE)$p.value
# p.adjust(s, "fdr")









# res.aov <- anova_test(
#   data=stat_data, dv = BW, wid = PID,between = age,
#   within = c(exp, load)
# )
# get_anova_table(res.aov)
# "NOVA Table (type III tests)

#         Effect DFn DFd     F     p p<.05      ges
# 1          age   1  54 1.846 0.180       0.016000
# 2          exp   1  54 1.762 0.190       0.005000
# 3         load   1  54 0.097 0.756       0.000258
# 4      age:exp   1  54 0.319 0.575       0.000892
# 5     age:load   1  54 0.401 0.529       0.001000
# 6     exp:load   1  54 0.281 0.599       0.001000
# 7 age:exp:load   1  54 0.050 0.824       0.000227"



# res.aov <- anova_test(
#   data=stat_data, dv = exponent, wid = PID,between = age,
#   within = c(exp, load)
# )
# get_anova_table(res.aov)
# "ANOVA Table (type III tests)

#         Effect DFn DFd       F        p p<.05      ges
# 1          age   1  54 101.746 5.03e-14     * 6.09e-01
# 2          exp   1  54   0.005 9.44e-01       7.01e-06
# 3         load   1  54   0.092 7.63e-01       9.29e-05
# 4      age:exp   1  54   0.001 9.74e-01       1.50e-06
# 5     age:load   1  54   3.855 5.50e-02       4.00e-03
# 6     exp:load   1  54   0.705 4.05e-01       5.39e-04
# 7 age:exp:load   1  54   0.032 8.59e-01       2.43e-05"
# s<-rep(1:2)
# a<-stat_data%>%filter(age==1&load=="1b")%>%select(exponent)
# b<-stat_data%>%filter(age==1&load=="2b")%>%select(exponent)
# s[1]<-t.test(a$exponent,b$exponent,paired = TRUE)$p.value
# a<-stat_data%>%filter(age==0&load=="1b")%>%select(exponent)
# b<-stat_data%>%filter(age==0&load=="2b")%>%select(exponent)
# s[2]<-t.test(a$exponent,b$exponent,paired = TRUE)$p.value
# p.adjust(s,"fdr")