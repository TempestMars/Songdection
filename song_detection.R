## song process

library(seewave)
library(tuneR)
library(spectral)


## cutting function, beware to change the the ffmpeg directioon
songsplit <- function(inputfile, start1, end1, type){
  count1 <- length(start1)
  char_count <- nchar(inputfile)
  filename_1 <- strsplit(inputfile, "/")
  filename_1 <- filename_1[[1]]
  filename_1 <- filename_1[length(filename_1)]
  filename_1 <- strsplit(filename_1, ".wav")[[1]][1]
  
  if(type == "motif"){
    for (i in 1 : (count1)){
      name1 <- paste("motif",as.character(i), "_", filename_1, ".wav", sep = "")
      system(
        paste(
          "C:/Users/bibio/Desktop/ffmpeg-20190426-0fc4646-win64-static/bin/ffmpeg.exe", 
          " -i ",
          inputfile,
          "-ss ",
          as.character(start1[i]/44100),
          " -to ",
          as.character(end1[i]/44100),
          " -c copy",
          name1
        ), 
        show.output.on.console = FALSE
      )
    }
  }
  if(type == "song"){
    for (i in 1 : (count1)){
      name1 <- paste("song",as.character(i), "_", filename_1, ".wav", sep = "")
      system(
        paste(
          "C:/Users/bibio/Desktop/ffmpeg-20190426-0fc4646-win64-static/bin/ffmpeg.exe", 
          " -i ",
          inputfile,
          "-ss ",
          as.character(start1[i]/44100),
          " -to ",
          as.character(end1[i]/44100),
          " -c copy",
          name1
        ), 
        show.output.on.console = FALSE)
    }
  }
  if(type == "syllable"){
    for (i in 1 : (count1)){
      name1 <- paste("syllable",as.character(i), "_", filename_1, ".wav", sep = "")
      system(
        paste(
          "C:/Users/bibio/Desktop/ffmpeg-20190426-0fc4646-win64-static/bin/ffmpeg.exe", 
          " -i ",
          inputfile,
          "-ss ",
          as.character(start1[i]/44100),
          " -to ",
          as.character(end1[i]/44100),
          " -c copy",
          name1
        ), 
        show.output.on.console = FALSE
      )
    }   
  }
}

## work flow

## function, segmentate songs according amplitude threshold and minimal interval. Can add a margin to the segment. 
song_segmentation <- function(path, threshold_song = 0.45, threshold_amp = 1.2, margin = 0.015, 
                              plot = 0.01, segt = TRUE, type = "song"){
  ## read the wav file
  path_1 <- path
  wav1 <- readWave(path_1)
  ## extract the acoustic signals from 250 to 12000 Hz
  pas1 <- ffilter(wav1,from = 250, to = 12000, bandpass = TRUE,
                  custom = NULL, wl = 1024, ovlp = 99, wn = "hanning", fftw = FALSE,
                  rescale=FALSE, listen=FALSE, output="Wave")
  
  ## extract amplitudes
  env1 <- env(pas1, envt = "abs", 
              msmooth = NULL, ksmooth = NULL, ssmooth = FALSE,
              fftw = FALSE, norm = FALSE,
              plot = FALSE, k = 1, j = 1)
  env1 <- (env1[-length(env1)])
  
  duration_song <- duration(pas1)
  env1 <- data.frame(env1)
  env1$counts <- c(1:nrow(env1))
  ## find all points that are above the amplitude threshold 
  thres1 <- threshold_song*44100
  env2 <- env1[env1$env1 >= threshold_amp,]
  counts1 <- c(0, env2$counts, duration_song*44100)
  n1 <- length(counts1)
  ## procede to the next loop if no sound detected
  if(n1 == 0){
    print("no sound detected")
    next
  }
  counts2 <- c(counts1[c(2:n1,n1)])
  
  ## based on the set interval, find suitable starts for cutting 
  cut_start_motif <- counts2[which(counts2 - counts1 >= thres1)] 
  ## get the whold recording if the detected sounds are continuous 
  if(length(cut_start_motif) == 0){
    cut_start_motif <- 0
    cut_end_motif <- duration_song*44100
    
  }else{
    if(length(cut_start_motif) == 1){
      cut_end_motif <- duration_song*44100
    }else{
      n2 <- length(cut_start_motif)
      cuts_temp <- counts1[which(counts2 - counts1 >= thres1)]
      cut_end_motif <- c(cuts_temp[1:n2], duration_song*44100)
      cut_start_motif <- c(0, cut_start_motif)  
    }
  }
  
  
  ## check if the cuts are too short
  index_cut <- which(cut_end_motif - cut_start_motif >= threshold_song*44100)
  cut_start_motif <- cut_start_motif[index_cut]
  cut_end_motif <- cut_end_motif[index_cut]
  n2 <- length(cut_start_motif)
  ## add margin
  if(n2 == 1){
    if(cut_start_motif[1] - margin*44100 < 0){
      cut_start_motif[1] <- 0
    }else{
      cut_start_motif <- (cut_start_motif - margin*44100)
    }
    if(cut_end_motif[1] + margin*44100 > duration_song*44100){
      cut_end_motif[1] <- duration_song*44100
    }else{
      cut_end_motif <- (cut_end_motif + margin*44100)
    }
    
  }
  if(n2 > 1){
    if(cut_start_motif[1] - margin*44100 < 0){
      cut_start_motif[1] <- 0
      cut_start_motif[2:n2] <- (cut_start_motif[2:n2] - margin*44100)
    }else{
      cut_start_motif <- (cut_start_motif - margin*44100)
    }
    if(cut_end_motif[n2] + margin*44100 > duration_song*44100){
      cut_end_motif[n2] <- duration_song*44100
      cut_end_motif[1:n2-1] <- (cut_end_motif[1:n2-1] + margin*44100)
    }else{
      cut_end_motif <- (cut_end_motif + margin*44100)
    }
  }
  
  ## depending on the input, plot the results
  plot_index <- runif(1)
  if(plot_index <= plot){
    spectro(wav1, fastdisp = TRUE,
            axisX = TRUE, axisY = TRUE, scale = FALSE, grid = FALSE,
            flim = c(0, 10))
    points(counts1/44100, y = rep(4, length(counts1)))
    abline(v = cut_start_motif/44100, col = "red")
    abline(v = cut_end_motif/44100, col = "blue")
  }
  ## depending on the input, cut the segements
  if(segt == TRUE){
    if(length(cut_start_motif) == 0){
      print(paste(path, ":", "no cutting point", sep = ""))
    }else{
      songsplit(path_1, cut_start_motif, cut_end_motif, type)
    }
  }
  gc()
}




## cut recordings to fine songs
## set output direction
setwd("C:/Users/bibio/Desktop/MS_4/506_cut")
## set inpput files
fn <- list.files("C:/Users/bibio/Desktop/MS_4/506", full.names = TRUE, pattern = ".wav")

lapply(fn, song_segmentation, 
       threshold_song = 0.45, threshold_amp = 2, margin = 0.015, 
       plot = 0, segt = TRUE, type = "song")


## test on one file
song_segmentation("H:/Song_recordings/father_songs/524/524Feb_10-6-129.wav", 
                  threshold_song = 0.15, threshold_amp = 1.2, margin = 0.015, 
                  plot = 1, segt = FALSE, type = "song")

library(purrr)

batch_song <- possibly(song_segmentation, otherwise = "something wrong here")
system.time(
  map(fn, batch_song)
)
song_segmentation(fn[1], plot = 1)

path <- "C:/Users/bibio/Desktop/New folder/motif1_M111.wav"
threshold <- 0.035

motifTsyllable <- function(path, threshold_syllable = 0.020, threshold_int = 0.005, threshold_amp = 1, margin = 0.008){
  not_process <- c()
  path_1 <- path
  wav1 <- readWave(path_1)
  pas1 <- ffilter(wav1,from = 250, to = 12000, bandpass = TRUE,
                  custom = NULL, wl = 1024, ovlp = 99, wn = "hanning", fftw = FALSE,
                  rescale=FALSE, listen=FALSE, output="Wave")
  
  
  env1 <- env(pas1, envt = "abs", 
              msmooth = NULL, ksmooth = NULL, ssmooth = FALSE,
              fftw = FALSE, norm = FALSE,
              plot = FALSE, k = 1, j = 1)
  duration_motif <- duration(pas1)*44100
  env1 <- (env1[-length(env1)])
  env1 <- data.frame(env1)
  env1$counts <- c(1:nrow(env1))
  
  thres1 <- threshold_int*44100
  env2 <- env1[env1$env1 >= 1, ] 
  counts1 <- env2$counts
  n1 <- length(counts1)
  counts2 <- c(counts1[c(2:n1,n1)])
  
  cut_start_syllable <- c(counts1[1], counts2[which(counts2 - counts1 >= thres1)]) 
  cut_end_syllable <- c(counts1[which(counts2 - counts1 >= thres1)], counts1[n1])
  n2 <- length(cut_start_syllable)
  if(is.na(counts1)){
    not_process <- c(not_process, path)
    print(not_process)
    break
  }
  
  if(counts1[1] - margin*44100 < 0){
    cut_start_syllable[1] <- 0
    cut_start_syllable[2:n2] <- (cut_start_syllable[2:n2] - margin*44100)
  }else{
    cut_start_syllable <- (cut_start_syllable - margin*44100)
  }
  if(counts1[n1] + margin*44100 > duration_motif){
    cut_end_syllable[n2] <- duration_motif
    cut_end_syllable[1:n2-1] <- (cut_end_syllable[1:n2-1] + margin*44100)
  }else{
    cut_end_syllable <- (cut_end_syllable + margin*44100)
  }
  
  
  spectro(wav1, fastdisp = FALSE,
          axisX = TRUE, axisY = TRUE, scale = FALSE, grid = FALSE,
          flim = c(0, 10))
  points(counts1/44100, rep(6, length(counts1)))
  abline(v = cut_start_syllable/44100, col = "blue")
  abline(v = cut_end_syllable/44100, col = "red")
  
  
  songsplit(path_1, cut_start_syllable, cut_end_syllable, "syllable")
}


par(mfrow=c(5,5))
par(mar=c(1,1,1,1))
spectro(wav1, fastdisp = TRUE,
        axisX = TRUE, axisY = TRUE, scale = FALSE, grid = FALSE,
        flim = c(0, 10))

## process songs

fn_song <- list.files("C:/Users/bibio/Desktop/MS4/SongM", pattern = ".wav",
                      full.names = TRUE)
setwd("C:/Users/bibio/Desktop/MS4/motifs_M")
system.time(
  lapply(fn_song[1:25], songTmotif, threshold_song = 0.5, threshold_int = 0.005, threshold_amp = 1.2, margin = 0.015,  
         plot = TRUE, segt = FALSE)
)


library(purrr)
function_test <- function(t,b = 2){
  t <- t + b
}
funtion_print <- function(x){
  print(x)
}
test_1 <- possibly(function_test, otherwise = "Wrong")
map("a", test_1, b = 3)

batch_song <- possibly(songTmotif, otherwise = "something wrong here")
system.time(
  map(fn_song, batch_song)
)
path <- fn_song[208]
t1 <- readWave(fn_song[208])
spectro(t1, flim = c(0, 8))

fn_song[c(89,122,125,129,131,142,143,200,202,204,207,208)]
fn_song[c(24,29,30,46,47,65,81,93,95,97,99,100, 154, 165)]
fn_song_notP <- fn_song[c(24,29,30,46,47,65,81,93,95,97,99,100, 154, 165)]


map(fn_song_notP, batch_song)


a1 <- system.time(fn_song <- list.files("C:/Users/bibio/Desktop/MS4/SongH", pattern = ".wav",
                                        full.names = TRUE))


fn_motifs <- list.files("C:/Users/bibio/Desktop/New folder", pattern = ".wav",
                        full.names = TRUE)
setwd("C:/Users/bibio/Desktop/MS4/syllables")
lapply(fn_motifs , motifTsyllable, threshold_syllable = 0.020, threshold_int = 0.005, threshold_amp = 1.2, margin = 0.008)




t1 <- readWave("C:/Users/bibio/Desktop/MS4/SongH/M122-1.wav")
spectro(t1)
wav1 <- t1


## get abs value

fn_motifs <- list.files("C:/Users/bibio/Desktop/New folder", pattern = ".wav",
                        full.names = TRUE)

range_record<-c()
for (i in 1 : length(fn_motifs)){
  
  path_1 <- fn_motifs[i]
  
  wav1 <- readWave(path_1)
  ## filter 
  pas1 <- ffilter(wav1,from = 250, to = 12000, bandpass = TRUE,
                  custom = NULL, wl = 1024, ovlp = 99, wn = "hanning", fftw = FALSE,
                  rescale=FALSE, listen=FALSE, output="Wave")
  
  ## get amplitude data
  env1 <- env(pas1, envt = "abs", 
              msmooth = NULL, ksmooth = NULL, ssmooth = FALSE,
              fftw = FALSE, norm = FALSE,
              plot = FALSE, k = 1, j = 1)
  
  env1 <- (env1[-length(env1)])
  print(range(env1))
  range_record <- c(range_record, range(env1)[2])
  
}



## timing 
wav1 <- readWave("//vuw/Personal$/Homes/L/liuq/Desktop/MS4/Acoustic_ML/song7.wav")
env1 <- env(wav1, envt = "hil", 
            msmooth = NULL, ksmooth = NULL, ssmooth = FALSE,
            fftw = FALSE, norm = FALSE,
            plot = FALSE, k = 1, j = 1)
env1[which(is.na(env1))] <- 0
envN_1 <-((env1 - min(env1))/(max(env1)-min(env1)))
plot(density(envN_1))
sum(envN_1)
envN_1_1 <- envN_1/sum(envN_1)
plot(envN_1_1)
sum(envN_1_1)
cum_env1 <- cumsum(envN_1_1)
which(cum_env1 >= 0.5)[1]/44100
length(env1)/44100
duration(wav1)
mean(envN_1_1)

str(envN_1_1)
framet <- data.frame(x = c(1:length(envN_1_1)), y = envN_1_1)

df <- approxfun(framet)
plot(envN_1_1)
xnew <- c(1:3000)
points(xnew,df(xnew),col = "red")

length(env1$env1)
hist(env1$env1)
hist((env1$env1 - mean(env1$env1))/(sum((env1$env1)^2)/length(env1$env1)))




add_silence <- function(file){
  file_1 <- file
  name_1 <- paste("S", basename(file), sep = "")   
  system(
    paste(  
      "C:/Users/bibio/Desktop/ffmpeg-20190426-0fc4646-win64-static/bin/ffmpeg.exe", 
      "-i C:/Users/bibio/Desktop/MS4/silence/silence_100ms.wav",
      "-i",
      file_1,
      "-i C:/Users/bibio/Desktop/MS4/silence/silence_100ms.wav -filter_complex \"[0:0][1:0][2:0]concat=n=3:v=0:a=1[out]\" -map \"[out]\"" ,
      name_1
    ),
    show.output.on.console = FALSE
  ) 
}

setwd("C:/Users/bibio/Desktop/MS4/syllable_M/Test_M/STest_M")
fn_1 <- list.files("C:/Users/bibio/Desktop/MS4/syllable_M/Test_M", full.names = TRUE)
lapply(fn_1, add_silence)

setwd("C:/Users/bibio/Desktop/MS4/syllable_H/Test_H/STest_H")
fn_1 <- list.files("C:/Users/bibio/Desktop/MS4/syllable_H/Test_H", full.names = TRUE)
lapply(fn_1, add_silence)


## fill the data

Fill_frame <- read_csv("C:/Users/bibio/Desktop/MS4/M_test.csv")
summary(Fill_frame)
mean_F2 <- mean(Fill_frame$F2[!is.na(Fill_frame$F2)])
Fill_frame$F2[is.na(Fill_frame$F2)] <- mean_F2

mean_cvfund <- mean(Fill_frame$cvfund[!is.na(Fill_frame$cvfund)])
Fill_frame$cvfund[is.na(Fill_frame$cvfund)] <- mean_cvfund

mean_fund <- mean(Fill_frame$fund[!is.na(Fill_frame$fund)])
Fill_frame$fund[is.na(Fill_frame$fund)] <- mean_fund

mean_maxfund <- mean(Fill_frame$maxfund[!is.na(Fill_frame$maxfund)])
Fill_frame$maxfund[is.na(Fill_frame$maxfund)] <- mean_maxfund

mean_minfund <- mean(Fill_frame$minfund[!is.na(Fill_frame$minfund)])
Fill_frame$minfund[is.na(Fill_frame$minfund)] <- mean_minfund
summary(Fill_frame)

write.csv(Fill_frame, "M_test.csv")
getwd()


fn <- list.files("//vuw/Personal$/Homes/L/liuq/Desktop/Projects/MS4/Songs_cut_H",
                 pattern  = ".wav", full.names = TRUE)


## check song or calls

index <- c()
for(i in 1 : length(fn)){
  wav1 <- readWave(fn[i])
  pas1 <- ffilter(wav1,from = 250, to = 12000, bandpass = TRUE,
                  custom = NULL, wl = 1024, ovlp = 99, wn = "hanning", fftw = FALSE,
                  rescale=FALSE, listen=FALSE, output="Wave")
  
  ## extract amplitudes
  env1 <- env(pas1, envt = "abs", 
              msmooth = NULL, ksmooth = NULL, ssmooth = FALSE,
              fftw = FALSE, norm = FALSE,
              plot = FALSE, k = 1, j = 1)
  env1 <- (env1[-length(env1)])
  
  duration_song <- duration(pas1)
  env1 <- data.frame(env1)
  env1$counts <- c(1:nrow(env1))
  ## find all points that are above the amplitude threshold 
  env2 <- env1[env1$env1 >= 1.2,]
  counts1 <- c(0, env2$counts, duration_song*44100)
  n1 <- length(counts1)
  
  counts2 <- c(counts1[c(2:n1,n1)])
  max1 <- max(counts2 - counts1)/44100
  if(max1 >= 0.15){
    index <- c(index, i)
  }
  
}

## visual check

pas1 <- ffilter(wav1,from = 250, to = 12000, bandpass = TRUE,
                custom = NULL, wl = 1024, ovlp = 99, wn = "hanning", fftw = FALSE,
                rescale=FALSE, listen=FALSE, output="Wave")

wav1 <- readWave("//vuw/Personal$/Homes/L/liuq/Desktop/call_song.wav")
spectro(wav1)

