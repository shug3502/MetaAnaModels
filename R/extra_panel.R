
#generate ana times df and complex positions etc
ana_times_df <- draws %>%
  mutate(t_ana = if_else(param=="t_ana",theta,NA_real_)) %>%
  group_by(cell,filename,SisterPairID) %>%
  summarise(t_ana = median(t_ana,na.rm=TRUE)) %>%
  ungroup() %>%
  mutate(cell=as.integer(cell),
         kittracking_file_str = filename %>% stringr::str_split('/') %>% purrr::map(last)) %>%
  tidyr::unnest()

#align to anaphase for each trajectory rather than the median for each cell
complex_positions_by_frame_df <- Data_2s %>%
  full_join(ana_times_df,by=c("kittracking_file_str","SisterPairID")) %>%
  mutate(Frames_since_start = Frame - t_ana%/%dt) %>%
  filter(Frames_since_start < nframes_early) %>%
  filter(Frames_since_start >= -nframes_early) %>%
  group_by(kittracking_file_str,Frames_since_start,SisterPairID,SisterID) %>%
  summarise(Position_1 = first(Position_1),
            Position_2 = first(Position_2),
            Position_3 = first(Position_3)) %>%
  group_by(kittracking_file_str) %>%
  filter(length(unique(SisterPairID))>min_num_sisters) %>%
  mutate(filename=kittracking_file_str)

dx=0.5
breaks=seq(from=-8,to=8,by=dx)
midpts <- seq(from=(breaks[1]+dx/2),to=(breaks[length(breaks)]-dx/2),by=dx)
complex_positions_by_frame_df %>%
filter(Frames_since_start>0) %>%
mutate(rad=sqrt(Position_2^2+Position_3^2)) %>% 
#group_by(SisterID,SisterPairID) %>%
#mutate(rad=rad/first(rad)) %>%
ungroup() %>%
mutate(x = cut(Position_1,breaks=breaks,labels=FALSE),
       pts = midpts[x]) %>%
group_by(pts) %>%
summarise(lb=quantile(rad,0.25,na.rm=T),ub=quantile(rad,0.75,na.rm=T),
	  rad=median(rad,na.rm=T)) %>%
ggplot(aes(x=pts,y=rad)) + geom_ribbon(aes(ymin=lb,ymax=ub),fill="grey85",color="black") + geom_line(size=2) + theme_bw() +
labs(x="x Position (um)",y="Radius (um)")
ggsave(here::here("plots/rad_fig_example_unnormalized.eps"),width=3,height=3)

