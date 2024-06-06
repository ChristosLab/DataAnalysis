library(R.matlab)
library(mgcv)
library(stats)
library(matlib)

# fn <- c("MSNG_avg_wfs.mat", "ODRdist_avg_wfs.mat", "ODRdistVar_avg_wfs.mat")
fn <- c("ODRdistVar_avg_wfs.mat")
# fn <- c("ODRdist_sig_missing_avg_wfs.mat")
spikewidths_fit <- numeric()
spikewidths_spline <- numeric()
ns <- 0
bs <- 0
all_fits <- matrix(nrow = 0, ncol = 41)
all_splinefits <- matrix(nrow = 0, ncol = 641)
# wfs <- numeric()
# fits <- numeric()
for (name in fn){
  mat_data <- readMat(name)
  
  # str(mat_data)
  time <- 0:69
  # results <- data.frame()
  for (n in 1:nrow(mat_data$all.avg.wf)){
  # for (n in c(257)){
    # for (n in 504){
    
    wf<-mat_data$all.avg.wf[n, 15:55]
    
    # plot(time[15:55],wf,type = 'n', xlab = 'Timepoints',
    #      ylab = 'Normalized Amplitude', ylim = c(-1, 1))
    # points(time[15:55], wf, pch = 'o')
    # print(wf)
    time_t = time[15:55]
    # mult <- 10
    df <- data.frame(time = time[15:55], wf = wf)
    gam_model <- gam(wf ~ s(time, bs = 'cr'), data = df)
    fit <- gam_model$fitted.values
    
    minidx <- 12
    maxidx <- which.max(fit==max(fit[minidx+1:length(wf)], na.rm=TRUE))
    width <- (maxidx-minidx)*25
    
    time_t = time[15:55]
    spline_result <- spline(time_t, fit, n = 641, 
                            method = "natural")
    spline_fit <- spline_result$y
    spline_time <- spline_result$x
    
    minidx <- 177 #12*mult #which.min(fit)
    maxidx <- which.max(spline_fit==max(spline_fit[minidx+1:length(spline_fit)],
                                        na.rm=TRUE))
    
    prev_window <- (length(fit) - 1 )*25
    curr_window <- (prev_window)/(length(spline_fit) - 1)
    
    width_c <- (maxidx-minidx)*curr_window
    if (abs(width-width_c) > 50) {
      print("Wrong")
    }
    spikewidths_fit <- c(spikewidths_fit, width)
    spikewidths_spline <- c(spikewidths_spline, width_c)
    # wfs <- c(wfs, wf)
    # fits <-c(fits, fit)
    # print((maxidx-minidx)*25)
    # plot(time_t,wf,type = 'n', xlab = 'Timepoints',
    #      ylab = 'Normalized Amplitude', ylim = c(-1, 1)) 
    # title(main=as.character(n))
    # points(time_t, wf, pch = 'o')
    # lines(time_t, fit, col='red', ylim = c(-1, 1), lwd=2)
    # lines(spline_time, spline_fit, col='blue', ylim = c(-1, 1), lwd=2)
    
    
    if (n==1){
      if (width<=300){
        # plot(spline_time, spline_fit, col='red', ylim = c(-1, 1),
        #      xlab = 'Timepoints', ylab = 'Normalized Amplitude')
        # title(main = substr(name, 1, nchar(name) - 12), col.main = "red", 
              # font.main = 4)
        # ns_wf <- wf
        # ns_fit <- spline_fit
        ns <- ns + 1
      }
      else{
        # plot(spline_time, spline_fit, col='lightblue', ylim = c(-1, 1),
        #      xlab = 'Timepoints', ylab = 'Normalized Amplitude')
        # title(main = substr(name, 1, nchar(name) - 12), col.main = "red", 
        #       font.main = 4)
        # bs_wf <- wf
        # bs_fit <- spline_fit
        bs <- bs + 1
      }

    }

    else{
    if (width<=300){
      # lines(spline_time, spline_fit, col='red', ylim = c(-1, 1), lwd=2)
      # ns_wf <- wf
      # ns_fit <- spline_fit
      ns <- ns + 1
    }
    else{
      # lines(spline_time, spline_fit, col='lightblue', ylim = c(-1, 1), lwd=2)
      # bs_wf <- wf
      # bs_fit <- spline_fit
      bs <- bs + 1
    }
    }
    print(n)
    # Sys.sleep(5)
  }
  # legend("topright",
  #        legend = c("BS", "NS"),
  #        col = c("lightblue", "red"),
  #        lwd = 2)
}

# -------------------------------------------------------------------------
# NS_BS_wf <-data.frame(ns_wf=ns_wf, ns_fit=ns_fit, bs_wf=bs_wf, bs_fit=bs_fit)
# writeMat("NS_BS_wf.mat",  x =NS_BS_wf)

# plot(time,wf,type = 'l', xlab = 'Timepoints', ylab = 'Normalized Amplitude')
# MSNG_widths <- spikewidths[1: 1436]
# ODRdist_widths <- spikewidths[(1436+1): (1436+578)]
# ODRdistVar_widths <- spikewidths[(1436+578+1): length(spikewidths)]
hist(spikewidths_spline,breaks = 10, main = "Histogram of Spikewidths",
     xlab = "Spikewidth", ylab = "Number of Neurons", col = "lightblue",
     border = "black")
# 
MSNGWidths <-data.frame(Spikewidths_fit=spikewidths_fit, 
                        spikewidths_spline=spikewidths_spline)
writeMat("ODRdistVar_widths_v1.mat",  x =MSNGWidths)

# AllWidths <-data.frame(Spikewidth = spikewidths)
# writeMat("AllWidths_26toMax.mat",  x =AllWidths)

