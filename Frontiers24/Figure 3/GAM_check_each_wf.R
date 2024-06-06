library(matlib)
library(R.matlab)
library(mgcv)
library(stats)

# fn <- c("MSNG_avg_wfs.mat", "ODRdist_avg_wfs.mat", "ODRdistVar_avg_wfs.mat")
fn <- c("Average_waveforms/GRU178_1_2019_avg_wf.mat")
spikewidths <- numeric()
all_min_idx <- numeric()
all_max_idx <- numeric()
# wfs <- numeric()
# fits <- numeric()
for (name in fn){
  mat_data <- readMat(name)
  
  # str(mat_data)
  time <- 0:69
  # results <- data.frame()
  # for (n in 101:150){#nrow(mat_data$all.avg.wf)){
  # for (n in c(1220)){
    # for (n in 504){
    n <- 1
    wf<-mat_data$mean.wf[n, 15:55]
    
    # plot(time[15:55],wf,type = 'n', xlab = 'Timepoints',
    # ylab = 'Normalized Amplitude', ylim = c(-1, 1))
    # points(time[15:55], wf, pch = 'o')
    # print(wf)
    
    df <- data.frame(time = time[15:55], wf = wf)
    
    gam_model <- gam(wf ~ s(time, bs = 'cr'), data = df)
    fit <- gam_model$fitted.values
    
    minidx <- which.min(fit)
    maxidx <- which.max(fit==max(fit[minidx+1:length(wf)], na.rm=TRUE))
    width <- (maxidx-minidx)*25
    spikewidths <- c(spikewidths, width)
    all_min_idx <- c(all_min_idx, minidx)
    all_max_idx <- c(all_max_idx, maxidx)
    
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
    # wfs <- c(wfs, wf)
    # fits <-c(fits, fit)
    # print((maxidx-minidx)*25)
    
    # if (n==1){
    #   if (width<=300){
    #     plot(time[15:55], fit, col='red', ylim = c(-1, 1),
    #          xlab = 'Timepoints', ylab = 'Normalized Amplitude')
    #     title(main = substr(name, 1, nchar(name) - 12), col.main = "red", font.main = 4)
    #   }
    #   else{
    #     plot(time[15:55], fit, col='lightblue', ylim = c(-1, 1),
    #          xlab = 'Timepoints', ylab = 'Normalized Amplitude')
    #     title(main = substr(name, 1, nchar(name) - 12), col.main = "red", font.main = 4)
    #   }
    # 
    # }
    # 
    # else{
    # if (width<=300){
      plot(time[15:55],wf,type = 'n', xlab = 'Timepoints',
           ylab = 'Normalized Amplitude', ylim = c(-1, 1)) 
      title(main=as.character(n))
      points(time[15:55], wf, pch = 'o')
      lines(time[15:55], fit, col='red', ylim = c(-1, 1), lwd=2)
      lines(spline_time, spline_fit, col='blue', ylim = c(-1, 1), lwd=2)
      
      ns_wf <- wf
      ns_fit <- fit
      print(n)
    # }
    # else{
    #   lines(time[15:55], fit, col='blue', ylim = c(-1, 1), lwd=2)
    #   bs_wf <- wf
    #   bs_fit <- fit
    # }
    # }
    # print(n)
    # Sys.sleep(5)
  # }
  # legend("topright",
  #        legend = c("BS", "NS"),
  #        col = c("lightblue", "red"),
  #        lwd = 2)
}

# -------------------------------------------------------------------------
NS_wf <-data.frame(ns_wf=ns_wf, ns_fit=ns_fit)
writeMat("NS_wf.mat",  x =NS_wf)

# plot(time,wf,type = 'l', xlab = 'Timepoints', ylab = 'Normalized Amplitude')
# MSNG_widths <- spikewidths[1: 1436]
# ODRdist_widths <- spikewidths[(1436+1): (1436+578)]
# ODRdistVar_widths <- spikewidths[(1436+578+1): length(spikewidths)]
# hist(spikewidths,breaks = 10, main = "Histogram of Spikewidths", 
#      xlab = "Spikewidth", ylab = "Number of Neurons", col = "lightblue", border = "black")
# 
# ODRdistWidths <-data.frame(MinVector=all_min_idx, MaxVector=all_max_idx, Spikewidth = spikewidths)
# writeMat("ODRdistWidthInfo.mat",  x =ODRdistWidths)

# AllWidths <-data.frame(Spikewidth = spikewidths)
# writeMat("AllWidths.mat",  x =AllWidths)

