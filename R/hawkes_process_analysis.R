hawkes_process_analysis <- function(sigma_sim,jobset_str,dt){
#fit hawkes process model to directional switching events 

Data <- process_jobset(jobset_str,K=Inf,max_missing=0.25)
K = max(Data$Frame)
pairIDs <- Data$SisterPairID %>% unique()
switch_probs_df <- get_switch_probs_df(sigma_sim,pairIDs,K)
xyt_df <- get_directional_switch_events(switch_probs_df, Data, dt=dt, min_prob=0.5) #only event with a posterior probability of more than 0.5

pairIDs <- xyt_df[["SisterPairID"]] %>% unique() %>% sort()
data <- list(Space=xyt_df[,c("Position_2","Position_3")],
     Time=xyt_df[["Time"]]-min(xyt_df[["Time"]]), 
     n=nrow(xyt_df),
     time_window=K*dt-min(xyt_df[["Time"]]), 
     space_window=pi*5^2,
     SisterPair=match(xyt_df[["SisterPairID"]],pairIDs),
     nPairs=length(pairIDs),
     prior_space=c(0,10),
     prior_time=c(0,1)
)
#used this in seth flaxman's code for hawkes process for gunshot contagion
m = stan_model(here::here("src/stan_files/hawkes_switching_analysis.stan"))
fit = sampling(m,data,warmup=500,iter=1000,chains=4)
print(fit,c("lengthscaleS","lengthscaleT","lengthscale","a","mu0"))

g <- plot_hawkes_process_analysis(xyt_df,fit,jobset_str)
h <- plot_kernel_space_time(fit, jobset_str)

#posterior=as.array(fit)
#np <- nuts_params(fit)
#library(bayesplot)
#scatter_theta_cp <- mcmc_scatter(
#  posterior, 
#  pars = c("lengthscaleS", "lengthscale"), 
#  np = np, 
#  size = 1
#)
#scatter_theta_cp
#mcmc_trace(posterior, pars = "mu0", np = np)

return(fit)
}

plot_hawkes_process_analysis <- function(xyt_df,fit,jobset_str){

job_id = stringr::str_split(jobset_str,"kittracking")[[1]][2]
edited_job_id = job_id %>% stringr::str_replace_all("\\.","")
out <- rstan::extract(fit)
llmask = colMeans(out$ll) > quantile(colMeans(out$ll),1-mean(out$a))
llshapes = rep(22,length(llmask))
llshapes[llmask == T] = 20
llalphas = rep(.25,length(llmask))
llalphas[llmask == T] = .6
df = data.frame(x=xyt_df$Position_2,y=xyt_df$Position_3,
                ll=factor(llmask,c("Background","Excitatory")),
		time=cut(xyt_df[["Time"]],seq(from=0,to=4*150,by=150)))
g = ggplot(data=df,aes(x,y))
g = g + geom_point(aes(shape=llshapes,alpha=llalphas),size=2)
g = g + facet_wrap(.~time,labeller="label_both") 
g = g + scale_shape_identity()
g = g + scale_alpha_identity()
g = g + theme_bw() + labs(x="y Position (um)",y="z Position (um)")
g
ggsave(here::here(paste0("plots/switching_background_vs_excitory_from_hawkes_process_",edited_job_id,".eps")),device=cairo_ps,width=4,height=4)
return(g)
}

plot_kernel_space_time <- function(fit, jobset_str){
job_id = stringr::str_split(jobset_str,"kittracking")[[1]][2]
edited_job_id = job_id %>% stringr::str_replace_all("\\.","")
out <- rstan::extract(fit)

#(a * lengthscaleT * exp(timeD[i,j] * lengthscaleT) * gaussian(spaceD[i,j],lengthscaleS))

pdf(here::here(paste0("plots/hawkes_kernels_",edited_job_id,".pdf")))
par(mfrow=c(1,2))
#NB scale from seconds to minutes
x=seq(from=0,to=5,by=0.01)
plot(x,out$a[1]*dexp(0,rate=out$lengthscaleT[1]*60)*dnorm(x, mean = 0, sd = out$lengthscaleS[1]),type='l',lwd=0.25,ylim=c(0,1.5),
     xlab="Distance (um)",ylab="Intensity (#/um min)")
for (iDraw in 2:100){
  lines(x,out$a[iDraw]*dexp(0,rate=out$lengthscaleT[iDraw]*60)*dnorm(x, mean = 0, sd = out$lengthscaleS[iDraw]),lwd=0.25)
}
lines(x,mean(out$a)*dexp(0,rate=mean(out$lengthscaleT*60))*dnorm(x,sd=mean(out$lengthscaleS)),lwd=4,col="magenta")

tt=seq(from=0,to=1,by=0.02)
plot(tt,out$a[1]*dexp(tt, rate=out$lengthscaleT[1]*60)*dnorm(0, mean = 0, sd = out$lengthscaleS[1]),type='l',lwd=0.25,ylim=c(0,1.5), 
     xlab="Time (min)",ylab="Intensity (#/um min)")
for (iDraw in 2:100){
  lines(tt,out$a[iDraw]*dexp(tt, rate = out$lengthscaleT[iDraw]*60)*dnorm(0, mean = 0, sd = out$lengthscaleS[iDraw]),lwd=0.25)
}
lines(tt,mean(out$a)*dexp(tt,rate=mean(out$lengthscaleT*60))*dnorm(0, mean = 0, sd = mean(out$lengthscaleS)),lwd=4,col="magenta")
dev.off()

my_kernel <- function(x,tt,a,lengthscaleT,lengthscaleS) a*dexp(tt, rate=lengthscaleT*60)*dnorm(x, mean = 0, sd = lengthscaleS)
grid_df <- expand.grid(x=x,time=tt) %>%
		mutate(intensity = my_kernel(x,time,mean(out$a),mean(out$lengthscaleT),mean(out$lengthscaleS)))
h <- ggplot(grid_df,aes(x=time,y=x,fill=intensity)) + 
geom_tile() + theme_bw() +
scale_fill_distiller(palette="Greys",direction=1) +
labs(x="Time (min)",y="Distance (um)",fill="Intensity\n(#/um min)")
h
ggsave(here::here(paste0("plots/hawkes_kernel_in_2d_",edited_job_id,".eps")),device=cairo_ps,width=4,height=4)
return(h)
}
