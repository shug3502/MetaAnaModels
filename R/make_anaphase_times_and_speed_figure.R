make_anaphase_times_and_speed_figure <- function(){
  #get all the interesting stats, and simple positions by frame df
green_val <- RColorBrewer::brewer.pal(4,"PiYG")[4]
plt1 <- ggplot(all_the_interesting_stats_df,
               aes(av_radius,relative_t_ana)) + 
  geom_point(alpha=0.3) + 
  geom_smooth(method="lm",se=FALSE,color=green_val) + 
  theme_bw() + 
  scale_x_continuous(limits=c(0,6)) +
  scale_y_continuous(limits=c(-100,100)) + 
  labs(x="Average radius (um)\nduring metaphase",
       y="Anaphase onset time (s)\nof each pair relative\nto median for the cell",
       title="KT pairs at the edge\n of the metaphase plate\nseparate earlier"
  )

plt2 <- ggplot(all_the_interesting_stats_df, 
               aes(av_radius,v_ana)) + 
  geom_point(alpha=0.3) + 
  geom_smooth(method="lm",se=FALSE,color=green_val) + 
  scale_x_continuous(limits=c(0,6)) +
  theme_bw() + 
  #  scale_y_continuous(limits=c(-100,100)) + 
  labs(x="Average radius (um)\n during metaphase",
       y="Anaphase speed (um/s)",
       title="KT pairs at the edge of\nthe metaphase plate\nhave lower speed v_A"
  )

plt3 <- simple_positions_by_frame_treatment_df %>% 
  filter(Frames_since_start > 0) %>% 
  group_by(kittracking_file_str,SisterPairID,SisterID,treatment) %>%
  mutate(framewise_dist_travelled=c(0,sqrt(diff(Position_1)^2+diff(Position_2)^2+diff(Position_3)^2))) %>%
  ungroup() %>% drop_na()%>%
  group_by(kittracking_file_str,SisterPairID,SisterID,treatment) %>%
  summarise(speed_3d=median(framewise_dist_travelled/dt,na.rm=T)) %>%
  inner_join(all_the_interesting_stats_df,by=c("kittracking_file_str", "SisterPairID", "SisterID", "treatment")) %>%
  filter(av_radius<=8) %>%
  ggplot(aes(av_radius,speed_3d)) + 
  geom_point(alpha=0.3) + 
  geom_smooth(method="lm",se=FALSE,color=green_val) + 
  theme_bw() + 
  labs(x="Average radius (um)\nduring metaphase",
       y="Average framewise speed\nin anaphase (um/s)",
       title="Framewise speed in\nanaphase does not show\nsame spatial effect")

plt1 | plt2 | plt3 
ggsave(here::here("plots/anaphase_time_and_speed_versus_radius.eps"),device=cairo_ps,width=210,height=150,units="mm")
}
