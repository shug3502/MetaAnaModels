
interpolate_missing_data <- function(Position,Time){
  Y <- zoo::zoo(cbind(Time,Position)) #see https://stackoverflow.com/questions/7188807/interpolate-na-values
  index(Y) <- Y[,1]
  Y_approx <- na.approx(Y)
  return(as.numeric(Y_approx[,2]))
}

extract_long_tracks <- function(Data,K=Inf,T0=0,max_missing=0){
  #take str giving reference to a jobset converted to csv and return all tracks beyond a certain length
  #defaults to using all the available data
  how_much_missing <- Data %>%
    dplyr::filter(Frame<=(T0+K), Frame>T0) %>% #assess only on first K frames, starting at frame T0 + 1
    dplyr::group_by(SisterPairID,SisterID) %>% #assess each track individually
    dplyr::summarise(proportionNaN = sum(is.na(Position_1))/length(Position_1)) %>%
    dplyr::group_by(SisterPairID) %>% #combine to consider pairs together
    dplyr::summarise(proportionNaN = max(proportionNaN))
  Data <- dplyr::left_join(Data,how_much_missing) %>%
    dplyr::group_by(SisterPairID,SisterID) %>%
    dplyr::filter(proportionNaN <= max_missing,
           Frame<=(T0+K), Frame>T0) %>%
    dplyr::ungroup()

  return(Data)
}

reorder_sisters <- function(Data){
  #ensures that sister 1 has on average the larger x (Position_1) coordinate
  Data %>%
    dplyr::group_by(SisterPairID,SisterID) %>%
    dplyr::summarise(x = mean(Position_1,na.rm=T)) %>%
    dplyr::group_by(SisterPairID) %>%
    dplyr::summarise(need_to_switch = (first(x) - last(x))<0) %>%
    dplyr::right_join(Data) %>%
    dplyr::mutate(SisterID = as.integer(dplyr::case_when(
      need_to_switch & (SisterID==1) ~ 2,
      need_to_switch & (SisterID==2) ~ 1,
      !need_to_switch ~ as.double(SisterID),
      TRUE ~ as.double(NA)))) %>%
    dplyr::select(-need_to_switch)
}

process_jobset <- function(jobset_str,K=Inf,max_missing=0,start_from=0,plot_opt=0){
  #this function makes it easier to read tracking output of tracked kinetochores
Data <- read.csv(jobset_str,header=TRUE) %>%
  dplyr::group_by(SisterPairID,SisterID) %>%
  dplyr::mutate(Position_1=interpolate_missing_data(Position_1,Time)) %>%
  dplyr::mutate(Position_2=interpolate_missing_data(Position_2,Time)) %>%
  dplyr::mutate(Position_3=interpolate_missing_data(Position_3,Time)) %>%
    dplyr::ungroup() %>%
    extract_long_tracks(K,start_from,max_missing) %>%
    reorder_sisters()
  if (plot_opt){
    g <- ggplot(Data, aes(x=Time, y=Position_1,color=factor(SisterID))) +
      geom_line() +
      facet_wrap(.~SisterPairID) +
      theme_bw() +
      theme(legend.position = "None",
      strip.text.x = element_text(size = 8),
      axis.text.x = element_text(angle=90)) +
      labs(x="Time (s)",y="Position (um)")
    print(g)
    ggsave(stringr::str_replace(jobset_str,".csv",
           paste("_processed_tracks_length",K,".eps",sep="")))
  }
return(Data)
}

ind_to_binary <- function(ind,nStates=5){
  aux = rep(0,nStates)
  aux[ind] = 1.0
  return(aux)
}

odeUpdateMatrix <- function(theta,nStates=5){
  M = matrix(rep(0,2*(nStates+2)),nrow=2)
  M[1,1] = -theta$kappa - theta$alpha;
  M[1,2] = theta$kappa;
  M[1,3] = -theta$v_plus;
  M[1,4] = -theta$v_plus;
  M[1,5] = -theta$v_minus;
  M[1,6] = -theta$v_minus;
  M[2,1] = theta$kappa;
  M[2,2] = -theta$kappa - theta$alpha;
  M[2,3] = theta$v_plus;
  M[2,4] = theta$v_minus;
  M[2,5] = theta$v_plus;
  M[2,6] = theta$v_minus;
  if (nStates==5){
    M[1,7] = 0;
    M[2,7] = 0;
  }
  return(M)
}

odeUpdateVector <- function(theta, angleTheta=0.0){
  mu = rep(0,2)
  mu[1] = theta$kappa*theta$L*cos(angleTheta);
  mu[2] = -theta$kappa*theta$L*cos(angleTheta);
  return(mu)
}

generate_synthetic_anaphase_data <- function(theta, sigma0, x0,
                                             K=170, dt=2.0, nStates=5){
  stopifnot(length(x0)==2)
  stopifnot(length(sigma0)==nStates)
  stopifnot(is.list(theta))
  tau = theta$tau
  p_icoh = theta$p_icoh
  p_coh = theta$p_coh
  v_ana = theta$v_ana
  t_ana = theta$t_ana

  P = matrix(rep(0,nStates^2),nrow=nStates)
  q_coh = 1 - p_coh
  q_icoh = 1 - p_icoh
  P[1,1] = p_icoh*p_icoh;
  P[2,1] = p_coh*q_coh;
  P[3,1] = p_coh*q_coh;
  P[4,1] = q_icoh*q_icoh;
  P[5,1] = 0;
  P[1,2] = p_icoh*q_icoh;
  P[2,2] = p_coh*p_coh;
  P[3,2] = q_coh*q_coh;
  P[4,2] = p_icoh*q_icoh;
  P[5,2] = 0;
  P[1,3] = p_icoh*q_icoh;
  P[2,3] = q_coh*q_coh;
  P[3,3] = p_coh*p_coh;
  P[4,3] = p_icoh*q_icoh;
  P[5,3] = 0;
  P[1,4] = q_icoh*q_icoh;
  P[2,4] = p_coh*q_coh;
  P[3,4] = p_coh*q_coh;
  P[4,4] = p_icoh*p_icoh;
  P[5,4] = 0;
  P[1,5] = 0;
  P[2,5] = 0;
  P[3,5] = 0;
  P[4,5] = 0;
  P[5,5] = 1;

  x = matrix(rep(0,2*K),nrow=2)
  sigma = matrix(rep(0,nStates*K),nrow=nStates)
  x[1:2,1] = x0
  sigma[,1] = sigma0
  for (t in 2:K){
    if (t*dt < t_ana){
      #generate hidden states
      sigma[,t] = ind_to_binary(sample.int(nStates,size=1,prob=P[which(sigma[,t-1]>0),]),nStates=nStates)
      #fwd euler step
      x[1:2,t] = x[1:2,t-1] + dt*odeUpdateMatrix(theta,nStates=5)%*%c(x[1:2,t-1],sigma[,t]) +
        dt*odeUpdateVector(theta) + rnorm(2)*sqrt(dt/tau)
    } else if (t*dt >= t_ana) {
      sigma[,t] = ind_to_binary(nStates,nStates=nStates)
      x[1:2,t] = x[1:2,t-1] + dt*c(v_ana,-v_ana) + rnorm(2)*sqrt(dt/tau)
    }
  }
  return(list(x=x,sigma=sigma))
}

get_single_pair <- function(Data,id) {
  Y1 <- Data %>% dplyr::filter(SisterPairID==id,
                        SisterID==1) %>%
    dplyr::pull(Position_1)
  Y2 <- Data %>% dplyr::filter(SisterPairID==id,
                        SisterID==2) %>%
    dplyr::pull(Position_1)
  Y = cbind(Y1,Y2)
  return(Y)
}

prepare_for_stan_format <- function(Data,IDs=NA){
  if (is.na(IDs)){
    #use default of all available sisters
    IDs <- unique(Data$SisterPairID)
  }
  Y_list <- purrr::map(IDs,function(id) get_single_pair(Data,id))
  return(Y_list)
}

get_start_end_of_nonmissing_data <- function(y_missing){
#convert a binary vector of whether data is missing into two integers T0 and T1
#which indicate the final missing point at the start, and the final existing point at the end
#T0>=0, T1<=K
######################
y_not_missing_ind <- which(y_missing<=0)
T0 = y_not_missing_ind[1] - 1
T1 = y_not_missing_ind[length(y_not_missing_ind)]
return(list(T0=T0,T1=T1))
}

compute_angle_phi <- function(Position_1,Position_2,Position_3){
#computes at a single time
#assumes Position_1 is a length 2 vector with positions for each sister
#returns cos of the angle from kk-axis to normal of metaphase plate
##########
stopifnot(length(Position_1)==2)
if (any(is.na(c(Position_1,Position_2,Position_3)))) {return(NA)}
inter_kt_vec = c(diff(Position_1),diff(Position_2),diff(Position_3))
cos_phi = inter_kt_vec[1]/norm(inter_kt_vec,type="2")
return(cos_phi)
}

get_cos_phi <- function(Data,id){
  cos_phi <- Data %>%
    dplyr::filter(SisterPairID==id) %>% group_by(Time) %>%
    summarise(cos_phi = compute_angle_phi(Position_1,Position_2,Position_3)) %>%
    pull(cos_phi)
}

