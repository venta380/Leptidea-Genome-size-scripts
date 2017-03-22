kmer<-read.table(file='/proj/b2014034/nobackup/jellyfish_genome_size_predict/genome_size_newdata_trimmed_K17_subsample_2/test.hist', header=FALSE)
sample_details<-read.table(file='/proj/b2014034/nobackup/jellyfish_genome_size_predict/genome_size_newdata_trimmed_K17_subsample_2/reads.txt', header=FALSE)
K=17
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


library(ggplot2)
Kmers <- kmer$V1
abundence <- kmer$V2
sample <- kmer$V3
group <- kmer$V4
data = data.frame(Kmers, abundence, sample, group)




sample_id <- sample_details$V1
read_pairs <- sample_details$V2
avg_read_length <- sample_details$V3
calculation_table = data.frame(sample_id, read_pairs, avg_read_length)


names(data)
Irish_juvernica <- subset(data, group=="Irish_Juvernica")
Swedish_sinapis <- subset(data, group=="Swedish_sinapis")
Sinapis_106 <- subset(data, group=="Sinapis_106")
Reali <- subset(data, group=="Reali")
Juvernica <- subset(data, group=="Juvernica")
Sinapis_56 <- subset(data, group=="Sinapis_56")


Sinapis_56_1 <- subset(Sinapis_56, sample=="1")
Sinapis_56_2 <- subset(Sinapis_56, sample=="2")
Sinapis_56_3 <- subset(Sinapis_56, sample=="3")
Sinapis_56_4 <- subset(Sinapis_56, sample=="6")
Sinapis_56_5 <- subset(Sinapis_56, sample=="7")
Sinapis_56_6 <- subset(Sinapis_56, sample=="8")
Sinapis_56_7 <- subset(Sinapis_56, sample=="9")
Sinapis_56_8 <- subset(Sinapis_56, sample=="10")
Sinapis_56_9 <- subset(Sinapis_56, sample=="11")
Sinapis_56_10 <- subset(Sinapis_56, sample=="12")

Sinapis_56_peakes <- c((subset(Sinapis_56_1, abundence == max(Sinapis_56_1$abundence[3:50])))$Kmers, (subset(Sinapis_56_2, abundence == max(Sinapis_56_2$abundence[3:50])))$Kmers, (subset(Sinapis_56_3, abundence == max(Sinapis_56_3$abundence[3:50])))$Kmers, (subset(Sinapis_56_4, abundence == max(Sinapis_56_4$abundence[3:50])))$Kmers, (subset(Sinapis_56_5, abundence == max(Sinapis_56_5$abundence[3:50])))$Kmers, (subset(Sinapis_56_6, abundence == max(Sinapis_56_6$abundence[3:50])))$Kmers, (subset(Sinapis_56_7, abundence == max(Sinapis_56_7$abundence[3:50])))$Kmers, (subset(Sinapis_56_8, abundence == max(Sinapis_56_8$abundence[3:50])))$Kmers, (subset(Sinapis_56_9, abundence == max(Sinapis_56_9$abundence[3:50])))$Kmers, (subset(Sinapis_56_10, abundence == max(Sinapis_56_10$abundence[3:50])))$Kmers )
Sinapis_56_samples <- c("Sample_1", "Sample_2", "Sample_3", "Sample_6", "Sample_7", "Sample_8", "Sample_9", "Sample_10", "Sample_11", "Sample_12")
Sinapis_56_read_pairs <- c((subset(calculation_table, sample_id == "1"))$read_pairs, (subset(calculation_table, sample_id == "2"))$read_pairs, (subset(calculation_table, sample_id == "3"))$read_pairs, (subset(calculation_table, sample_id == "6"))$read_pairs, (subset(calculation_table, sample_id == "7"))$read_pairs, (subset(calculation_table, sample_id == "8"))$read_pairs, (subset(calculation_table, sample_id == "9"))$read_pairs, (subset(calculation_table, sample_id == "10"))$read_pairs, (subset(calculation_table, sample_id == "11"))$read_pairs, (subset(calculation_table, sample_id == "12"))$read_pairs) 
Sinapis_56_avg_read_length <- c((subset(calculation_table, sample_id == "1"))$avg_read_length, (subset(calculation_table, sample_id == "2"))$avg_read_length, (subset(calculation_table, sample_id == "3"))$avg_read_length, (subset(calculation_table, sample_id == "6"))$avg_read_length, (subset(calculation_table, sample_id == "7"))$avg_read_length, (subset(calculation_table, sample_id == "8"))$avg_read_length, (subset(calculation_table, sample_id == "9"))$avg_read_length, (subset(calculation_table, sample_id == "10"))$avg_read_length, (subset(calculation_table, sample_id == "11"))$avg_read_length, (subset(calculation_table, sample_id == "12"))$avg_read_length) 
Sinapis_56_total_bases <- c(Sinapis_56_avg_read_length*(Sinapis_56_read_pairs*2))
Sinapis_56_depth <- c((Sinapis_56_peakes*Sinapis_56_avg_read_length)/(Sinapis_56_avg_read_length-(K+1)))
Sinapis_56_G_size <- c(Sinapis_56_total_bases/Sinapis_56_depth)

Sinapis_56_read_pairs

Irish_juvernica1 <- subset(Irish_juvernica, sample=="Ire-juv-1C")
Irish_juvernica2 <- subset(Irish_juvernica, sample=="Ire-juv-21C")
Irish_juvernica3 <- subset(Irish_juvernica, sample=="Ire-juv-22C")
Irish_juvernica4 <- subset(Irish_juvernica, sample=="Ire-juv-2C")
Irish_juvernica5 <- subset(Irish_juvernica, sample=="Ire-juv-41C")
Irish_juvernica6 <- subset(Irish_juvernica, sample=="Ire-juv-42C")
Irish_juvernica7 <- subset(Irish_juvernica, sample=="Ire-juv-61C")
Irish_juvernica8 <- subset(Irish_juvernica, sample=="Ire-juv-62C")
Irish_juvernica9 <- subset(Irish_juvernica, sample=="Ire-juv-81C")
Irish_juvernica10 <- subset(Irish_juvernica, sample=="Ire-juv-82C")

Irish_juvernica_peakes <- c((subset(Irish_juvernica1, abundence == max(Irish_juvernica1$abundence[3:50])))$Kmers, (subset(Irish_juvernica2, abundence == max(Irish_juvernica2$abundence[3:50])))$Kmers, (subset(Irish_juvernica3, abundence == max(Irish_juvernica3$abundence[3:50])))$Kmers, (subset(Irish_juvernica4, abundence == max(Irish_juvernica4$abundence[3:50])))$Kmers, (subset(Irish_juvernica5, abundence == max(Irish_juvernica5$abundence[3:50])))$Kmers, (subset(Irish_juvernica6, abundence == max(Irish_juvernica6$abundence[3:50])))$Kmers, (subset(Irish_juvernica7, abundence == max(Irish_juvernica7$abundence[3:50])))$Kmers, (subset(Irish_juvernica8, abundence == max(Irish_juvernica8$abundence[3:50])))$Kmers, (subset(Irish_juvernica9, abundence == max(Irish_juvernica9$abundence[3:50])))$Kmers, (subset(Irish_juvernica10, abundence == max(Irish_juvernica10$abundence[3:50])))$Kmers )
Irish_juvernica_samples <- c("Ire-juv-1C", "Ire-juv-21C", "Ire-juv-22C", "Ire-juv-2C", "Ire-juv-41C", "Ire-juv-42C", "Ire-juv-61C", "Ire-juv-62C", "Ire-juv-81C", "Ire-juv-82C")
Irish_juvernica_read_pairs <- c((subset(calculation_table, sample_id == "Ire-juv-1C"))$read_pairs, (subset(calculation_table, sample_id == "Ire-juv-21C"))$read_pairs, (subset(calculation_table, sample_id == "Ire-juv-22C"))$read_pairs, (subset(calculation_table, sample_id == "Ire-juv-2C"))$read_pairs, (subset(calculation_table, sample_id == "Ire-juv-41C"))$read_pairs, (subset(calculation_table, sample_id == "Ire-juv-42C"))$read_pairs, (subset(calculation_table, sample_id == "Ire-juv-61C"))$read_pairs, (subset(calculation_table, sample_id == "Ire-juv-62C"))$read_pairs, (subset(calculation_table, sample_id == "Ire-juv-81C"))$read_pairs, (subset(calculation_table, sample_id == "Ire-juv-82C"))$read_pairs) 
Irish_juvernica_avg_read_length <- c((subset(calculation_table, sample_id == "Ire-juv-1C"))$avg_read_length, (subset(calculation_table, sample_id == "Ire-juv-21C"))$avg_read_length, (subset(calculation_table, sample_id == "Ire-juv-22C"))$avg_read_length, (subset(calculation_table, sample_id == "Ire-juv-2C"))$avg_read_length, (subset(calculation_table, sample_id == "Ire-juv-41C"))$avg_read_length, (subset(calculation_table, sample_id == "Ire-juv-42C"))$avg_read_length, (subset(calculation_table, sample_id == "Ire-juv-61C"))$avg_read_length, (subset(calculation_table, sample_id == "Ire-juv-62C"))$avg_read_length, (subset(calculation_table, sample_id == "Ire-juv-81C"))$avg_read_length, (subset(calculation_table, sample_id == "Ire-juv-82C"))$avg_read_length) 
Irish_juvernica_total_bases <- c(Irish_juvernica_avg_read_length*(Irish_juvernica_read_pairs*2))
Irish_juvernica_depth <- c((Irish_juvernica_peakes*Irish_juvernica_avg_read_length)/(Irish_juvernica_avg_read_length-(K+1)))
Irish_juvernica_G_size <- c(Irish_juvernica_total_bases/Irish_juvernica_depth)
Irish_juvernica_G_size

Sinapis_106_1 <- subset(Sinapis_106, sample=="37")
Sinapis_106_2 <- subset(Sinapis_106, sample=="38")
Sinapis_106_3 <- subset(Sinapis_106, sample=="39")
Sinapis_106_4 <- subset(Sinapis_106, sample=="40")
Sinapis_106_5 <- subset(Sinapis_106, sample=="41")
Sinapis_106_6 <- subset(Sinapis_106, sample=="42")
Sinapis_106_7 <- subset(Sinapis_106, sample=="43")
Sinapis_106_8 <- subset(Sinapis_106, sample=="44")
Sinapis_106_9 <- subset(Sinapis_106, sample=="46")
Sinapis_106_10 <- subset(Sinapis_106, sample=="47")

Sinapis_106_peakes <- c((subset(Sinapis_106_1, abundence == max(Sinapis_106_1$abundence[3:50])))$Kmers, (subset(Sinapis_106_2, abundence == max(Sinapis_106_2$abundence[3:50])))$Kmers, (subset(Sinapis_106_3, abundence == max(Sinapis_106_3$abundence[3:50])))$Kmers, (subset(Sinapis_106_4, abundence == max(Sinapis_106_4$abundence[3:50])))$Kmers, (subset(Sinapis_106_5, abundence == max(Sinapis_106_5$abundence[3:50])))$Kmers, (subset(Sinapis_106_6, abundence == max(Sinapis_106_6$abundence[3:50])))$Kmers, (subset(Sinapis_106_7, abundence == max(Sinapis_106_7$abundence[3:50])))$Kmers, (subset(Sinapis_106_8, abundence == max(Sinapis_106_8$abundence[3:50])))$Kmers, (subset(Sinapis_106_9, abundence == max(Sinapis_106_9$abundence[3:50])))$Kmers, (subset(Sinapis_106_10, abundence == max(Sinapis_106_10$abundence[3:50])))$Kmers )
Sinapis_106_samples <- c("Sample_37", "Sample_38", "Sample_39", "Sample_40", "Sample_41", "Sample_42", "Sample_43", "Sample_44", "Sample_46", "Sample_47")
Sinapis_106_read_pairs <- c((subset(calculation_table, sample_id == "37"))$read_pairs, (subset(calculation_table, sample_id == "38"))$read_pairs, (subset(calculation_table, sample_id == "39"))$read_pairs, (subset(calculation_table, sample_id == "40"))$read_pairs, (subset(calculation_table, sample_id == "41"))$read_pairs, (subset(calculation_table, sample_id == "42"))$read_pairs, (subset(calculation_table, sample_id == "43"))$read_pairs, (subset(calculation_table, sample_id == "44"))$read_pairs, (subset(calculation_table, sample_id == "46"))$read_pairs, (subset(calculation_table, sample_id == "47"))$read_pairs) 
Sinapis_106_avg_read_length <- c((subset(calculation_table, sample_id == "37"))$avg_read_length, (subset(calculation_table, sample_id == "38"))$avg_read_length, (subset(calculation_table, sample_id == "39"))$avg_read_length, (subset(calculation_table, sample_id == "40"))$avg_read_length, (subset(calculation_table, sample_id == "41"))$avg_read_length, (subset(calculation_table, sample_id == "42"))$avg_read_length, (subset(calculation_table, sample_id == "43"))$avg_read_length, (subset(calculation_table, sample_id == "44"))$avg_read_length, (subset(calculation_table, sample_id == "46"))$avg_read_length, (subset(calculation_table, sample_id == "47"))$avg_read_length) 
Sinapis_106_total_bases <- c(Sinapis_106_avg_read_length*(Sinapis_106_read_pairs*2))
Sinapis_106_depth <- c((Sinapis_106_peakes*Sinapis_106_avg_read_length)/(Sinapis_106_avg_read_length-(K+1)))
Sinapis_106_G_size <- c(Sinapis_106_total_bases/Sinapis_106_depth)
Sinapis_106_G_size


Swedish_sinapis1 <- subset(Swedish_sinapis, sample=="Swe-sin-101C")
Swedish_sinapis2 <- subset(Swedish_sinapis, sample=="Swe-sin-102C")
Swedish_sinapis3 <- subset(Swedish_sinapis, sample=="Swe-sin-1C")
Swedish_sinapis4 <- subset(Swedish_sinapis, sample=="Swe-sin-2C")
Swedish_sinapis5 <- subset(Swedish_sinapis, sample=="Swe-sin-31C")
Swedish_sinapis6 <- subset(Swedish_sinapis, sample=="Swe-sin-32C")
Swedish_sinapis7 <- subset(Swedish_sinapis, sample=="Swe-sin-61C")
Swedish_sinapis8 <- subset(Swedish_sinapis, sample=="Swe-sin-62C")
Swedish_sinapis9 <- subset(Swedish_sinapis, sample=="Swe-sin-91C")
Swedish_sinapis10 <- subset(Swedish_sinapis, sample=="Swe-sin-92C")


Swedish_sinapis_peakes <- c((subset(Swedish_sinapis1, abundence == max(Swedish_sinapis1$abundence[3:50])))$Kmers, (subset(Swedish_sinapis2, abundence == max(Swedish_sinapis2$abundence[3:50])))$Kmers, (subset(Swedish_sinapis3, abundence == max(Swedish_sinapis3$abundence[3:50])))$Kmers, (subset(Swedish_sinapis4, abundence == max(Swedish_sinapis4$abundence[3:50])))$Kmers, (subset(Swedish_sinapis5, abundence == max(Swedish_sinapis5$abundence[3:50])))$Kmers, (subset(Swedish_sinapis6, abundence == max(Swedish_sinapis6$abundence[3:50])))$Kmers, (subset(Swedish_sinapis7, abundence == max(Swedish_sinapis7$abundence[3:50])))$Kmers, (subset(Swedish_sinapis8, abundence == max(Swedish_sinapis8$abundence[3:50])))$Kmers, (subset(Swedish_sinapis9, abundence == max(Swedish_sinapis9$abundence[3:50])))$Kmers, (subset(Swedish_sinapis10, abundence == max(Swedish_sinapis10$abundence[3:50])))$Kmers )
Swedish_sinapis_samples <- c("Swe-sin-101C", "Swe-sin-102C", "Swe-sin-1C", "Swe-sin-2C", "Swe-sin-31C", "Swe-sin-32C", "Swe-sin-61C", "Swe-sin-62C", "Swe-sin-91C", "Swe-sin-92C")
Swedish_sinapis_read_pairs <- c((subset(calculation_table, sample_id == "Swe-sin-101C"))$read_pairs, (subset(calculation_table, sample_id == "Swe-sin-102C"))$read_pairs, (subset(calculation_table, sample_id == "Swe-sin-1C"))$read_pairs, (subset(calculation_table, sample_id == "Swe-sin-2C"))$read_pairs, (subset(calculation_table, sample_id == "Swe-sin-31C"))$read_pairs, (subset(calculation_table, sample_id == "Swe-sin-32C"))$read_pairs, (subset(calculation_table, sample_id == "Swe-sin-61C"))$read_pairs, (subset(calculation_table, sample_id == "Swe-sin-62C"))$read_pairs, (subset(calculation_table, sample_id == "Swe-sin-91C"))$read_pairs, (subset(calculation_table, sample_id == "Swe-sin-92C"))$read_pairs) 
Swedish_sinapis_avg_read_length <- c((subset(calculation_table, sample_id == "Swe-sin-101C"))$avg_read_length, (subset(calculation_table, sample_id == "Swe-sin-102C"))$avg_read_length, (subset(calculation_table, sample_id == "Swe-sin-1C"))$avg_read_length, (subset(calculation_table, sample_id == "Swe-sin-2C"))$avg_read_length, (subset(calculation_table, sample_id == "Swe-sin-31C"))$avg_read_length, (subset(calculation_table, sample_id == "Swe-sin-32C"))$avg_read_length, (subset(calculation_table, sample_id == "Swe-sin-61C"))$avg_read_length, (subset(calculation_table, sample_id == "Swe-sin-62C"))$avg_read_length, (subset(calculation_table, sample_id == "Swe-sin-91C"))$avg_read_length, (subset(calculation_table, sample_id == "Swe-sin-92C"))$avg_read_length) 
Swedish_sinapis_total_bases <- c(Swedish_sinapis_avg_read_length*(Swedish_sinapis_read_pairs*2))
Swedish_sinapis_depth <- c((Swedish_sinapis_peakes*Swedish_sinapis_avg_read_length)/(Swedish_sinapis_avg_read_length-(K+1)))
Swedish_sinapis_G_size <- c(Swedish_sinapis_total_bases/Swedish_sinapis_depth)
Swedish_sinapis_G_size



Reali_1 <- subset(Reali, sample=="26")
Reali_2 <- subset(Reali, sample=="27")
Reali_3 <- subset(Reali, sample=="29")
Reali_4 <- subset(Reali, sample=="30")
Reali_5 <- subset(Reali, sample=="31")
Reali_6 <- subset(Reali, sample=="32")
Reali_7 <- subset(Reali, sample=="33")
Reali_8 <- subset(Reali, sample=="34")
Reali_9 <- subset(Reali, sample=="35")
Reali_10 <- subset(Reali, sample=="36")


Reali_peakes <- c((subset(Reali_1, abundence == max(Reali_1$abundence[3:50])))$Kmers, (subset(Reali_2, abundence == max(Reali_2$abundence[3:50])))$Kmers, (subset(Reali_3, abundence == max(Reali_3$abundence[3:50])))$Kmers, (subset(Reali_4, abundence == max(Reali_4$abundence[3:50])))$Kmers, (subset(Reali_5, abundence == max(Reali_5$abundence[3:50])))$Kmers, (subset(Reali_6, abundence == max(Reali_6$abundence[3:50])))$Kmers, (subset(Reali_7, abundence == max(Reali_7$abundence[3:50])))$Kmers, (subset(Reali_8, abundence == max(Reali_8$abundence[3:50])))$Kmers, (subset(Reali_9, abundence == max(Reali_9$abundence[3:50])))$Kmers, (subset(Reali_10, abundence == max(Reali_10$abundence[3:50])))$Kmers )
Reali_samples <- c("Sample_26", "Sample_27", "Sample_29", "Sample_30", "Sample_31", "Sample_32", "Sample_33", "Sample_34", "Sample_35", "Sample_36")
Reali_read_pairs <- c((subset(calculation_table, sample_id == "26"))$read_pairs, (subset(calculation_table, sample_id == "27"))$read_pairs, (subset(calculation_table, sample_id == "29"))$read_pairs, (subset(calculation_table, sample_id == "30"))$read_pairs, (subset(calculation_table, sample_id == "31"))$read_pairs, (subset(calculation_table, sample_id == "32"))$read_pairs, (subset(calculation_table, sample_id == "33"))$read_pairs, (subset(calculation_table, sample_id == "34"))$read_pairs, (subset(calculation_table, sample_id == "35"))$read_pairs, (subset(calculation_table, sample_id == "36"))$read_pairs) 
Reali_avg_read_length <- c((subset(calculation_table, sample_id == "26"))$avg_read_length, (subset(calculation_table, sample_id == "27"))$avg_read_length, (subset(calculation_table, sample_id == "29"))$avg_read_length, (subset(calculation_table, sample_id == "30"))$avg_read_length, (subset(calculation_table, sample_id == "31"))$avg_read_length, (subset(calculation_table, sample_id == "32"))$avg_read_length, (subset(calculation_table, sample_id == "33"))$avg_read_length, (subset(calculation_table, sample_id == "34"))$avg_read_length, (subset(calculation_table, sample_id == "35"))$avg_read_length, (subset(calculation_table, sample_id == "36"))$avg_read_length) 
Reali_total_bases <- c(Reali_avg_read_length*(Reali_read_pairs*2))
Reali_depth <- c((Reali_peakes*Reali_avg_read_length)/(Reali_avg_read_length-(K+1)))
Reali_G_size <- c(Reali_total_bases/Reali_depth)
Reali_G_size




Juvernica_1 <- subset(Juvernica, sample=="13")
Juvernica_2 <- subset(Juvernica, sample=="14")
Juvernica_3 <- subset(Juvernica, sample=="15")
Juvernica_4 <- subset(Juvernica, sample=="17")
Juvernica_5 <- subset(Juvernica, sample=="18")
Juvernica_6 <- subset(Juvernica, sample=="20")
Juvernica_7 <- subset(Juvernica, sample=="21")
Juvernica_8 <- subset(Juvernica, sample=="22")
Juvernica_9 <- subset(Juvernica, sample=="23")
Juvernica_10 <- subset(Juvernica, sample=="24")




Juvernica_peakes <- c((subset(Juvernica_1, abundence == max(Juvernica_1$abundence[3:50])))$Kmers, (subset(Juvernica_2, abundence == max(Juvernica_2$abundence[3:50])))$Kmers, (subset(Juvernica_3, abundence == max(Juvernica_3$abundence[3:50])))$Kmers, (subset(Juvernica_4, abundence == max(Juvernica_4$abundence[3:50])))$Kmers, (subset(Juvernica_5, abundence == max(Juvernica_5$abundence[3:50])))$Kmers, (subset(Juvernica_6, abundence == max(Juvernica_6$abundence[3:50])))$Kmers, (subset(Juvernica_7, abundence == max(Juvernica_7$abundence[3:50])))$Kmers, (subset(Juvernica_8, abundence == max(Juvernica_8$abundence[3:50])))$Kmers, (subset(Juvernica_9, abundence == max(Juvernica_9$abundence[3:50])))$Kmers, (subset(Juvernica_10, abundence == max(Juvernica_10$abundence[3:50])))$Kmers )
Juvernica_samples <- c("Sample_13", "Sample_14", "Sample_15", "Sample_17", "Sample_18", "Sample_20", "Sample_21", "Sample_22", "Sample_23", "Sample_24")
Juvernica_read_pairs <- c((subset(calculation_table, sample_id == "13"))$read_pairs, (subset(calculation_table, sample_id == "14"))$read_pairs, (subset(calculation_table, sample_id == "15"))$read_pairs, (subset(calculation_table, sample_id == "17"))$read_pairs, (subset(calculation_table, sample_id == "18"))$read_pairs, (subset(calculation_table, sample_id == "20"))$read_pairs, (subset(calculation_table, sample_id == "21"))$read_pairs, (subset(calculation_table, sample_id == "22"))$read_pairs, (subset(calculation_table, sample_id == "23"))$read_pairs, (subset(calculation_table, sample_id == "24"))$read_pairs) 
Juvernica_avg_read_length <- c((subset(calculation_table, sample_id == "13"))$avg_read_length, (subset(calculation_table, sample_id == "14"))$avg_read_length, (subset(calculation_table, sample_id == "15"))$avg_read_length, (subset(calculation_table, sample_id == "17"))$avg_read_length, (subset(calculation_table, sample_id == "18"))$avg_read_length, (subset(calculation_table, sample_id == "20"))$avg_read_length, (subset(calculation_table, sample_id == "21"))$avg_read_length, (subset(calculation_table, sample_id == "22"))$avg_read_length, (subset(calculation_table, sample_id == "23"))$avg_read_length, (subset(calculation_table, sample_id == "24"))$avg_read_length) 
Juvernica_total_bases <- c(Juvernica_avg_read_length*(Juvernica_read_pairs*2))
Juvernica_depth <- c((Juvernica_peakes*Juvernica_avg_read_length)/(Juvernica_avg_read_length-(K+1)))
Juvernica_G_size <- c(Juvernica_total_bases/Juvernica_depth)
Juvernica_G_size




genome_sizes <- c(Irish_juvernica_G_size/1000000, Swedish_sinapis_G_size/1000000, Sinapis_106_G_size/1000000, Reali_G_size/1000000, Juvernica_G_size/1000000, Sinapis_56_G_size/1000000)
depth <- c(Irish_juvernica_depth, Swedish_sinapis_depth, Sinapis_106_depth, Reali_depth, Juvernica_depth, Sinapis_56_depth)
populations <- c(rep("LjIre", 10), rep("LsSwe", 10) , rep("LsSpa", 10), rep("LrSpa", 10), rep("LjKaz", 10), rep("LsKaz", 10))
samples <- c(Irish_juvernica_samples, Swedish_sinapis_samples, Sinapis_106_samples, Reali_samples, Juvernica_samples, Sinapis_56_samples )
colours <- c(rep("#C8C800",10),rep("#006400",10), rep("#0000FF",10), rep("#FF8C00",10), rep("#C04000",10), rep("#FF0000",10))
final_data = data.frame(genome_sizes, depth, populations, samples, colours)

final_data$colours= factor(final_data$colours,levels = c("#C04000", "#FF8C00", "#FF0000", "#0000FF", "#C8C800", "#006400"),ordered = TRUE)
final_data$populations= factor(final_data$populations,levels = c("LsSwe", "LsKaz", "LsSpa", "LrSpa", "LjIre", "LjKaz"),ordered = TRUE)


final_table_scaled<-subset(final_data ,final_data$population=='LsSwe')


scaling_val1=mean(final_table_scaled$genome_sizes)
scaling_val2=642750272/1000000

genome_sizes_scaled=c(((final_data$genome_sizes)/scaling_val1)*scaling_val2)

scaled_new=data.frame(final_data,genome_sizes_scaled)

final_data=scaled_new

irish_plot <- ggplot(Irish_juvernica, aes(x=Kmers, y=abundence, colour=group, group=sample)) + 
  geom_line(color='firebrick1') + geom_vline(xintercept = mean(Irish_juvernica_peakes), colour="firebrick1") +
  ggtitle("K-mer plot for Irish juvernica")
irish_plot

swe_sin_plot <- ggplot(Swedish_sinapis, aes(x=Kmers, y=abundence, colour=group, group=sample)) + 
  geom_line(color='hotpink') + geom_vline(xintercept = (mean(Swedish_sinapis_peakes)), colour="hotpink") +
  ggtitle("K-mer plot for Swedish Sinapis")
swe_sin_plot

Sinapis_106_plot <- ggplot(Sinapis_106, aes(x=Kmers, y=abundence, colour=group, group=sample)) + 
  geom_line(color='cyan') + geom_vline(xintercept = (mean(Sinapis_106_peakes)), colour="cyan") +
  ggtitle("K-mer plot for Sinapis_106")
Sinapis_106_plot

Reali_plot <- ggplot(Reali, aes(x=Kmers, y=abundence, colour=group, group=sample)) + 
  geom_line(color='forestgreen') + geom_vline(xintercept = (mean(Reali_peakes)), colour="forestgreen") +
  ggtitle("K-mer plot for Reali")
Reali_plot

Juvernica_plot <- ggplot(Juvernica, aes(x=Kmers, y=abundence, colour=group, group=sample)) + 
  geom_line(color='gold4') + geom_vline(xintercept = (mean(Juvernica_peakes)), colour="gold4") +
  ggtitle("K-mer plot for Juvernica")
Juvernica_plot

Sinapis_56_plot <- ggplot(Sinapis_56, aes(x=Kmers, y=abundence, colour=group, group=sample)) + 
  geom_line(color='cornflowerblue') + geom_vline(xintercept = (mean(Sinapis_56_peakes)), colour="cornflowerblue") +
  ggtitle("K-mer plot for Sinapis_56")
Sinapis_56_plot

#combigned_plot <- ggplot(data, aes(x=Kmers, y=abundence, colour=group, group=sample)) + geom_vline(x = (mean(Irish_juvernica_peakes)), colour="firebrick1") + geom_vline(x = (mean(Swedish_sinapis_peakes)), colour="hotpink") + geom_vline(x = (mean(Sinapis_106_peakes)), colour="cyan") + geom_vline(x = (mean(Sinapis_56_peakes)), colour="cornflowerblue") + geom_vline(x = (mean(Reali_peakes)), colour="forestgreen") + geom_vline(x = (mean(Juvernica_peakes)), colour="gold4") + 
#  geom_line() +
#  ggtitle("K-mer plot (combigned)")
#combigned_plot

boxplot1 <- ggplot(final_data, aes(populations, depth)) + geom_boxplot((aes(fill = populations))) + ggtitle("depth comparisions (K = 17, All trimmed data pooled)") + labs(x = "Populations", y = "Depth [X]")+ scale_fill_manual(values = levels(final_table$colours), labels= levels(final_table$populations)) + guides(fill=FALSE)+theme_bw()
boxplot1









ALL_R<-read.table(file='/proj/b2014034/private/repeats/ALL_R.txt', header=FALSE)


library(ggplot2)
library(ggfortify)
library(cluster)


samples <- ALL_R$V1
LTR <- ALL_R$V2
LINES <- ALL_R$V3
SINES <- ALL_R$V4
DNA <- ALL_R$V5
Helitron <- ALL_R$V6

total <- ALL_R$V12
population <- ALL_R$V13

LTR_fractions <- LTR/total
LINES_fractions <- LINES/total
SINES_fractions <- SINES/total
DNA_fractions <- DNA/total
Helitron_fractions <- Helitron/total

calculation_table = data.frame(LTR_fractions, LINES_fractions, SINES_fractions, DNA_fractions, Helitron_fractions, samples, population)

final_table = merge(calculation_table, final_data, )




pca.data <- data.frame( final_table$LTR_fractions, final_table$LINES_fractions, final_table$SINES_fractions , final_table$DNA_fractions,  final_table$Helitron_fractions,final_table$population, final_table$samples)
pca.data

juvernica <- subset(pca.data, pca.data$final_table.population=="IRISH_JUVERNICA"  | pca.data$final_table.population=="KAZAK_JUVERNICA")
juvernica



boxplot <- ggplot(final_table, aes(populations, genome_sizes_scaled)) + geom_boxplot((aes(fill = populations))) + ggtitle("") + labs(x = "", y = "Genome size (Mb)")+ scale_fill_manual(values = levels(final_table$colours), labels= levels(final_table$populations)) + guides(fill=FALSE)+theme_bw(base_size = 18)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


boxplot_LTR <- ggplot(final_table, aes(populations, LTR_fractions)) + geom_boxplot((aes(fill = (populations)))) + ggtitle("") + labs(x = "", y = "LTR fraction")+       scale_fill_manual(values = levels(final_table$colours), labels= levels(final_table$populations)) + guides(fill=FALSE)+theme_bw(base_size = 18)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
boxplot_LINE <- ggplot(final_table, aes(populations, LINES_fractions)) + geom_boxplot((aes(fill = (populations)))) + ggtitle("") + labs(x = "", y = "LINES fraction") + scale_fill_manual(values = levels(final_table$colours), labels= levels(final_table$populations)) + guides(fill=FALSE)+theme_bw(base_size = 18)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
boxplot_SINE <- ggplot(final_table, aes(populations, SINES_fractions)) + geom_boxplot((aes(fill = (populations)))) + ggtitle("") + labs(x = "", y = "SINES fraction")+  scale_fill_manual(values = levels(final_table$colours), labels= levels(final_table$populations)) + guides(fill=FALSE)+theme_bw(base_size = 18)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
boxplot_DNA <- ggplot(final_table, aes(populations, DNA_fractions)) + geom_boxplot((aes(fill = (populations)))) + ggtitle("") + labs(x = "", y = "DNA fraction") + scale_fill_manual(values = levels(final_table$colours), labels= levels(final_table$populations)) + guides(fill=FALSE)+theme_bw(base_size = 18)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
boxplot_Helitron <- ggplot(final_table, aes(populations, Helitron_fractions)) + geom_boxplot((aes(fill = (populations)))) + ggtitle("") + labs(x = "", y = "Helitron fraction")+  scale_fill_manual(values = levels(final_table$colours), labels= levels(final_table$populations)) + guides(fill=FALSE)+theme_bw(base_size = 18)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())



corelation_LTR = ggplot(final_table, aes(x=genome_sizes_scaled, y=LTR_fractions )) + geom_point(aes(colour = populations))+ geom_smooth(method=lm,formula=y ~ x) + ggtitle("")    +labs(x = "Genome size (Mb)", y = "LTR fraction")+ scale_colour_manual(values = levels(final_table$colours)) +theme_bw(base_size = 18)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+theme(legend.title = element_blank(), legend.key = element_rect(colour = "white"))
corelation_LINE = ggplot(final_table, aes(x=genome_sizes_scaled, y=LINES_fractions )) + geom_point(aes(colour = populations))+ geom_smooth(method=lm,formula=y ~ x) + ggtitle("") +labs(x = "Genome size (Mb)", y = "LINES fraction")+ scale_colour_manual(values = levels(final_table$colours)) +theme_bw(base_size = 18)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(legend.position='none')
corelation_SINE = ggplot(final_table, aes(x=genome_sizes_scaled, y=SINES_fractions )) + geom_point(aes(colour = populations))+ geom_smooth(method=lm,formula=y ~ x) + ggtitle("") +labs(x = "Genome size (Mb)", y = "SINES fraction")+ scale_colour_manual(values = levels(final_table$colours)) +theme_bw(base_size = 18)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(legend.position='none')
corelation_DNA = ggplot(final_table, aes(x=genome_sizes_scaled, y=DNA_fractions )) + geom_point(aes(colour = populations))+ geom_smooth(method=lm,formula=y ~ x) + ggtitle("") +labs(x = "Genome size (Mb)", y = "DNA fraction")+ scale_colour_manual(values = levels(final_table$colours)) +theme_bw(base_size = 18)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(legend.position='none')
corelation_Helitron = ggplot(final_table, aes(x=genome_sizes_scaled, y=Helitron_fractions )) + geom_point(aes(colour = populations))+ geom_smooth(method=lm,formula=y ~ x) + ggtitle("") +labs(x = "Genome size (Mb)", y = "Helitron fraction")+ scale_colour_manual(values = levels(final_table$colours)) +theme_bw(base_size = 18)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(legend.position='none')



#corelation_SINEvsLINE = ggplot(final_table, aes(x=LINES_fractions, y=SINES_fractions )) + geom_point(aes(colour = population))+ geom_smooth(method=lm,formula=y ~ x) + ggtitle("SINES vs LINES") 
#autoplot(prcomp(log(pca.data[, 1:3])), data = pca.data, colour = "final_table.population", frame = TRUE, frame.type='norm')+ggtitle("") + scale_colour_manual(values = c( "#C8C800","#006400","#FF8C00","#0000FF","#FF0000","#C04000")) +scale_fill_manual(values = c( "#C8C800","#006400","#FF8C00","#0000FF","#FF0000","#C04000")) + guides(colour=FALSE, fill=FALSE)+ theme_bw(base_size = 18)
pca_guide=autoplot(prcomp(log(pca.data[, 1:5])), data = pca.data, colour = "final_table.population", loadings = TRUE,loadings.label = TRUE,frame = TRUE, shape = FALSE, label.size = 3)+ggtitle("") + scale_colour_manual(values = c("#C8C800","#006400","#FF8C00","#0000FF","#FF0000","#C04000")) + scale_fill_manual(values = c( "#C8C800","#006400","#FF8C00","#0000FF","#FF0000","#C04000")) + guides(colour=FALSE ,fill=FALSE)+ theme_bw(base_size = 18)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
pca=autoplot(prcomp(log(pca.data[, 1:5])), data = pca.data, colour = "final_table.population",frame = TRUE)+ggtitle("") + scale_colour_manual(values = c( "#C8C800","#006400","#FF8C00","#0000FF","#FF0000","#C04000")) + scale_fill_manual(values = c( "#C8C800","#006400","#FF8C00","#0000FF","#FF0000","#C04000")) + guides(colour=FALSE ,fill=FALSE)+ theme_bw(base_size = 18)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
framed_juvernica=autoplot(prcomp(log(juvernica[, 1:5])), data = juvernica, colour = "final_table.population",frame = TRUE)+ggtitle("") + scale_colour_manual(values = c( "#C8C800","#006400")) + scale_fill_manual(values = c( "#C8C800","#006400")) + guides(colour=FALSE ,fill=FALSE)+ theme_bw(base_size = 18)+ theme_bw(base_size = 18)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

write.csv(final_table ,'/home/venkat/glob/genome_size/final_table')
write.csv(juvernica ,'/home/venkat/glob/genome_size/juvernica')



boxplot_LTR
boxplot_LINE
boxplot_SINE
corelation_LTR
corelation_LINE
corelation_SINE



#multiplot(boxplot, boxplot_LINE, boxplot_SINE, corelation_LINE, corelation_SINE, framed, cols=3)






#autoplot(prcomp(log(juvernica[, 1:3])), data = juvernica, colour = "final_table.population", frame = TRUE, frame.type = 'norm')
#framed_juvernica = autoplot(prcomp(log(juvernica[, 1:3])), data = juvernica, colour = "final_table.population", frame = TRUE)+ggtitle("PCA for Juvernica populations ")
#autoplot(prcomp(log(pca.data[, 1:3])), data =framed_juvernica , colour = "final_table.population", frame = TRUE, frame.type='norm')

framed


pdf(file= 'myOut_land.pdf' ,onefile=T,paper='A4', )
for (i in 1){
  plot(boxplot)
  plot(boxplot_LTR)
  plot(boxplot_LINE)
  plot(boxplot_SINE)
  plot(boxplot_DNA)
  plot(boxplot_Helitron)
  plot(corelation_LTR)
  plot(corelation_LINE)
  plot(corelation_SINE)
  plot(corelation_DNA)
  plot(corelation_Helitron)
  plot(pca_guide)
  plot(pca)
  plot(framed_juvernica)
}
dev.off()

