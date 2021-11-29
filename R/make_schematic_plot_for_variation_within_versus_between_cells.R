library(dplyr)
library(ggplot2)
library(patchwork)
library(mvtnorm)
#https://stackoverflow.com/questions/50303683/overlay-two-contours-of-bivariate-gaussian-distribution-using-ggplot2
set.seed(123)
m1 <- c(8,12)
sigma1 <- matrix(c(14,-6,-6,7), nrow=2)
m2 <- c(12,8)
sigma2 <- matrix(c(14,6,6,7), nrow=2)
data.grid <- expand.grid(s.1 = seq(-25, 25, length.out=200), s.2 = seq(-25, 
                                                                       25, length.out=200))
q2.samp = cbind(data.grid, prob = mvtnorm::dmvnorm(data.grid, mean = m2, sigma=sigma2))
q1.samp = cbind(data.grid, prob = mvtnorm::dmvnorm(data.grid, mean = m1, sigma=sigma1))

theta1 <- rmvnorm(10, mean = m1, sigma = sigma1)
theta2 <- rmvnorm(10, mean = m2, sigma = sigma2)
df <- as.data.frame(rbind(theta1,theta2)) %>% mutate(cell=c(rep(1,10),rep(2,10)))

plt1 <- ggplot(df) + 
  geom_contour(data=q1.samp,aes(x=s.1,y=s.2,z=prob),col=RColorBrewer::brewer.pal(n=4,"PiYG")[2]) +    
  geom_contour(data=q2.samp,aes(x=s.1,y=s.2,z=prob),col=RColorBrewer::brewer.pal(n=4,"PiYG")[3]) + 
  geom_point(aes(V1,V2,color=factor(cell),
                 shape=factor(cell)),size=5) + 
  theme_classic() + 
  theme(axis.title = element_text(size=20)) +
  scale_x_discrete(labels = NULL) +
  scale_y_discrete(labels = NULL) +
  scale_color_manual(values=RColorBrewer::brewer.pal(n=4,"PiYG")[c(1,4)]) + 
  labs(x=expression(theta[1]),y=expression(theta[2]),color="cell",shape="cell")

plt_cell1 <- ggplot(df %>% filter(cell==1)) + 
  geom_contour(data=q1.samp,aes(x=s.1,y=s.2,z=prob),col=RColorBrewer::brewer.pal(n=4,"PiYG")[2]) +    
  geom_point(aes(V1,V2,color=factor(cell),
                 shape=factor(cell)),size=5) + 
  theme_classic() + 
  theme(axis.title = element_text(size=20)) +
  scale_x_discrete(labels = NULL) +
  scale_y_discrete(labels = NULL) +
  scale_color_manual(values=RColorBrewer::brewer.pal(n=4,"PiYG")[c(1,4)]) + 
  labs(x=expression(theta[1]),y=expression(theta[2]),color="cell",shape="cell")

plt_cell2 <- ggplot(df %>% filter(cell==2)) + 
  geom_contour(data=q2.samp,aes(x=s.1,y=s.2,z=prob),col=RColorBrewer::brewer.pal(n=4,"PiYG")[3]) + 
  geom_point(aes(V1,V2,color=factor(cell),
                 shape=factor(cell)),size=5) + 
  theme_classic() + 
  theme(axis.title = element_text(size=20)) +
  scale_x_discrete(labels = NULL) +
  scale_y_discrete(labels = NULL) +
  scale_color_manual(values=RColorBrewer::brewer.pal(n=4,"PiYG")[c(4)]) + 
  labs(x=expression(theta[1]),y=expression(theta[2]),color="cell",shape="cell")

plt_between_cells <- ggplot(df) + 
  geom_contour(data=q1.samp,aes(x=s.1,y=s.2,z=prob),
               col=RColorBrewer::brewer.pal(n=4,"PiYG")[2],
               alpha=0.2) +    
  geom_contour(data=q2.samp,aes(x=s.1,y=s.2,z=prob),
               col=RColorBrewer::brewer.pal(n=4,"PiYG")[3],
               alpha=0.2) + 
  geom_point(aes(V1,V2,color=factor(cell),
                 shape=factor(cell)),size=5,alpha=0.3) + 
  geom_point(data=df %>% group_by(cell) %>%
               summarise(V1=mean(V1),
                         V2=mean(V2)),
             aes(V1,V2,color=factor(cell),
                 shape=factor(cell)),size=7,alpha=1) +
  theme_classic() + 
  theme(axis.title = element_text(size=20)) +
  scale_x_discrete(labels = NULL) +
  scale_y_discrete(labels = NULL) +
  scale_color_manual(values=RColorBrewer::brewer.pal(n=4,"PiYG")[c(1,4)]) + 
  labs(x=expression(theta[1]),y=expression(theta[2]),color="cell",shape="cell")

bottom_row <- (plt_cell1 | plt_cell2 | plt_between_cells) +  plot_layout(widths = c(1,1,2))
{plot_spacer() | plt1}/bottom_row
ggsave(here::here("plots/within_vs_between_cells_simplified.eps"),device=cairo_ps,width=210,height=150,units="mm")
