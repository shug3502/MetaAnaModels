make_tracking_figure <- function(){
jobset_str_list <- list.files(path="data/",pattern="*.csv",
                              recursive = TRUE,full.names = TRUE)
Data_2s <- purrr::map(jobset_str_list,
                      function(x) process_jobset(x,max_missing=0.95,K=Inf,
                                                 plot_opt=FALSE)) %>%
  bind_rows(.id="cell") %>%
  mutate(filename=jobset_str_list[as.integer(cell)])

Data_raw <- purrr::map(jobset_str_list,
                       function(x) read.csv(x,header=TRUE)) %>%
  bind_rows(.id="cell") %>%
  mutate(filename=jobset_str_list[as.integer(cell)])

tracked_summary_df <- Data_raw %>%
  group_by(filename,SisterPairID,SisterID) %>%
  summarise(missing_time_pts=sum(is.na(Position_1)),
            frames_for_cell=max(Frame),
            proportion_tracked=1-missing_time_pts/frames_for_cell)

df <- tracked_summary_df %>% group_by(filename) %>%
  arrange(1-proportion_tracked) %>%
  mutate(rn = row_number())
ddx <- 0.1
num_kts_plt <- ggplot(df,aes(x=proportion_tracked,y=rn,group=filename)) +
  geom_step(alpha=0.5,color=RColorBrewer::brewer.pal(n=4,"PiYG")[c(2)]) +
  geom_hline(aes(yintercept=92),color="grey",linetype="dashed") +
  geom_step(data=df %>% mutate(proportion_binned=cut(proportion_tracked,
                                                     right=FALSE,
                                                     breaks=seq(from=0,to=1+ddx,by=ddx))) %>%
              group_by(proportion_binned) %>%
              summarise(rn=median(rn)),
            aes(x=c(seq(from=ddx/2,to=1-ddx/2,by=ddx),1),y=rn,group="1")) +
  theme_bw() +
  theme(legend.position="none") +
  scale_y_continuous(limits=c(0,130)) +
  labs(x="Proportion of movie",
       y="Number of kinetochores")

tracked_summary_df <- Data_raw %>%
  group_by(filename,SisterPairID,Frame) %>%
  summarise(is_missing=any(is.na(Position_1))) %>%
  group_by(filename,SisterPairID) %>%
  summarise(missing_time_pts=sum(is_missing),
            frames_for_cell=max(Frame),
            proportion_tracked=1-missing_time_pts/frames_for_cell)

df <- tracked_summary_df %>% group_by(filename) %>%
  arrange(1-proportion_tracked) %>%
  mutate(rn = row_number())

num_pairs_plt <- ggplot(df,aes(x=proportion_tracked,y=rn,group=filename)) +
  geom_step(alpha=0.5,color=RColorBrewer::brewer.pal(n=4,"PiYG")[c(2)]) +
  geom_hline(aes(yintercept=46),color="grey",linetype="dashed") +
  geom_step(data=df %>% mutate(proportion_binned=cut(proportion_tracked,
                                                     right=FALSE,
                                                     breaks=seq(from=0,to=1+ddx,by=ddx))) %>%
              group_by(proportion_binned) %>%
              summarise(rn=median(rn)),
            aes(x=c(seq(from=ddx/2,to=1-ddx/2,by=ddx),1),y=rn,group="1")) +
  theme_bw() +
  theme(legend.position="none") +
  scale_y_continuous(limits=c(0,65)) +
  labs(x="Proportion of movie",
       y="Number of pairs")

# #Look at this to select nicely tracked cell
# Data_2s %>% group_by(filename,SisterPairID) %>%
#   summarise(proportionNaN=first(proportionNaN),nFrames=max(Frame)) %>%
#   ungroup() %>% filter(proportionNaN<0.5) %>%
#   group_by(filename) %>%
#   summarise(nPairs=n(),nFrames=first(nFrames)) %>% View()

 jobset_str <- Data_2s$filename %>% unique() %>% nth(13)
position_time_plt <- Data_2s %>%
  filter(filename==jobset_str) %>%
  ggplot(aes((Frame-1)*dt,Position_1,color=factor(SisterID),
             group=interaction(SisterID,SisterPairID))) +
  geom_line(alpha=0.5) +
  scale_color_manual(values=RColorBrewer::brewer.pal(n=4,"PiYG")[c(2,4)]) +
  theme_bw() +
  labs(x="Time (s)",y="x Position (um)",color="SisterID") +
  theme(legend.position="none")

plate_position_plt <- Data_2s %>%
  filter(filename==jobset_str) %>%
  filter(Frame<=281) %>%
  filter(Frame>=182) %>%
  ggplot(aes(Position_2,Position_3,color=factor(SisterID))) +
  geom_line(aes(group=interaction(SisterPairID,SisterID)),alpha=0.5) +
  geom_point(data = Data_2s %>%
               filter(filename==jobset_str) %>%
               filter(Frame==281), alpha=1,size=2) +
  scale_color_manual(values=RColorBrewer::brewer.pal(n=4,"PiYG")[c(2,4)]) +
  theme_bw() +
  labs(x="y Position (um)",y="z Position (um)",color="SisterID") +
  theme(legend.position="none")

av_radius_df <- Data_2s %>%
  filter(filename==jobset_str) %>%
  group_by(SisterPairID) %>%
  summarise(av_position_2=median(Position_2,na.rm=T),
            av_position_3=median(Position_3,na.rm=T),
            `y (um)` = cut(av_position_2,seq(from=-6,to=6,by=1)),
            `z (um)` = cut(av_position_3,seq(from=-6,to=6,by=1)))

tracks_within_metaphase_plate_plt <- Data_2s %>%
  filter(filename==jobset_str) %>%
  inner_join(av_radius_df,by="SisterPairID") %>%
  filter(proportionNaN<0.5) %>%
  ggplot(aes((Frame-1)*dt,Position_1,color=factor(SisterID),
             group=interaction(SisterPairID,SisterID))) +
  geom_line() +
  facet_grid(`z (um)`~`y (um)`,
             labeller=function(x) label_both(x,sep="\n")) +
  # theme_bw() +
  labs(x="Time (s)",y="x Position (um)",color="SisterID") +
  theme(legend.position="bottom") +
  theme_classic() +
  theme(strip.background=element_blank(),
        text = element_text(size = 8),
        axis.text = element_text(size = 8),
        legend.position = "none",
        strip.text = element_text(size = 8)) +
  scale_color_manual(values=RColorBrewer::brewer.pal(n=4,"PiYG")[c(2,4)])

{num_kts_plt | num_pairs_plt | plate_position_plt | position_time_plt} /
  tracks_within_metaphase_plate_plt + plot_layout(heights=c(1,2))
ggsave(here::here("plots/tracking_fig_panel.eps"),device=cairo_ps,
       width=210,height=180,units="mm")
}
################
