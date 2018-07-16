
library(Hmisc)
library(Rgraphviz)
library(graph)
library(ggplot2)
library(gtable)
library(grid)
library(gridExtra)
library(dplyr)
library(tidyr)

source("kendrick.mass.filter_new_all_datamatrix_terminologyupdate.R")

data_matrix <- read.csv("all17K/tab17PEG.csv", header=T, check.names=F)
head(data_matrix)

vec <- do.call(rbind,lapply(strsplit(colnames(data_matrix)[2:ncol(data_matrix)],";"),matrix,ncol=2,byrow=TRUE))
class(vec) <- "numeric"
colnames(vec)<- c("mz","rt")
head(vec)

output <- Kendrick.mass.filter(
    data_matrix, 
    vec,
    polymer = "polyethylene_glycol_other_C2H4O1", 
    KMD = 0.01, 
    RT = 0.8, 
    NOS = 2, 
    connection_filter = TRUE)

output$Kendrickfiltered_MS1features_wgraph

kendrickmassfilterinfo <- as.data.frame(output$kendrickmassfilterinfo_original)
kendrickmassfilterinfo_filtered_wgraph <- as.data.frame(output$kendrickmassfilterinfo_filtered_wgraph)
kendrickmassfilterinfo_diff <- kendrickmassfilterinfo[((kendrickmassfilterinfo[,1] %in% 
                                                        kendrickmassfilterinfo_filtered_wgraph[,1]) & 
                                                       (kendrickmassfilterinfo[,2] %in% 
                                                        kendrickmassfilterinfo_filtered_wgraph[,2])) != TRUE,]

Kendrickplot <- 
    ggplot()+
    geom_point(data=as.data.frame(output$kendrickmassfilterinfo_filtered_wgraph), aes(x=nom_kend, y=msdefect), pch=16, 
               size=2.5, alpha=0.9, col="black")+
    geom_point(data=kendrickmassfilterinfo_diff,aes(x=nom_kend,y=msdefect), pch=16, 
               size=2.5, alpha=0.9,  col="#e41a1c")+
    scale_y_continuous(breaks = seq(-0.5,0.5,0.1))+
    scale_x_continuous(breaks = seq(0,1000,100))+
    xlab("Kendrick nominal m/z") +
    ylab("Kendrick mass defect (KMD)") +
    theme_minimal() +
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
               axis.line=element_line(colour ="black",size=0.5, linetype="solid"),
               axis.ticks=element_line(colour ="black",size=0.5, linetype="solid"),
               panel.background = element_rect(fill="white", colour="white", size=0.5, linetype="solid"),
               legend.position="bottom",
               legend.background = element_rect(fill="white",size=0.25,linetype="blank"),
               legend.title = element_text("Plasma + Swab Extract Features Filtered"),
               legend.text=element_text(size=6)) +
    theme(aspect.ratio=1)
print(Kendrickplot)
ggsave("Kendrickmassplot.pdf")

#colnames(output$Kendrickfiltered_MS1features_wgraph)
KMF_vec <- do.call(rbind,lapply(strsplit(colnames(output$Kendrickfiltered_MS1features_wgraph)[2:ncol(output$Kendrickfiltered_MS1features_wgraph)],";"),matrix,ncol=2,byrow=TRUE))
class(KMF_vec) <- "numeric"
colnames(KMF_vec)<- c("mz","rt")
head(KMF_vec)

Filtered_matrix_KMF_wgraph_MS1_plot_vec <- as.data.frame(vec[((vec[,1] %in% KMF_vec[,1]) & (vec[,2] %in% KMF_vec[,2])) != TRUE,])

MS1featureplot <- 
    ggplot()+
    geom_point(data=as.data.frame(vec),aes(x=rt, y=mz), pch=16, 
              size=2.5, alpha=0.9, color="black")+
    geom_point(data=Filtered_matrix_KMF_wgraph_MS1_plot_vec, aes(x=rt, y=mz), pch=16, 
               size=2.5, alpha=0.9, col="#e41a1c")+
    scale_y_continuous(breaks = seq(0,1000,100))+
    scale_x_continuous(breaks = seq(0,10,1))+
        xlab("retention time (min)") +
        ylab("m/z") +
    theme_minimal() +
        theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
               axis.line=element_line(colour ="black",size=0.5, linetype="solid"),
               axis.ticks=element_line(colour ="black",size=0.5, linetype="solid"),
               panel.background = element_rect(fill="white", colour="white", size=0.5, linetype="solid"),
               legend.position="bottom",
               legend.background = element_rect(fill="white",size=0.25,linetype="blank"),
               legend.title = element_text("Plasma + Swab Extract Features Filtered"),
               legend.text=element_text(size=6)) +
    theme(aspect.ratio=1)
print(MS1featureplot)
ggsave("MS1featureplot.pdf")

original_spectrum <- as.data.frame(cbind(vec,as.numeric(t(data_matrix)[2:nrow(t(data_matrix)),2])))
label <- format(round(vec[,1], 2), nsmall = 4)
original_spectrum <- cbind(original_spectrum, label)
colnames(original_spectrum) <- c("mz","rt","abundance","mz_label")
head(original_spectrum)

# Original spectrum
Original_Spectrum_plot <- 
    ggplot()+
    geom_bar(data=original_spectrum,aes(x=mz, y=abundance), stat="identity",  
              size=0.5, alpha=1, color="black")+
    geom_text(data=original_spectrum,aes(x=mz, y=abundance,label=mz_label), cex=3, check_overlap = TRUE, vjust=0) +
   scale_x_continuous(breaks = seq(0,1000,100))+
        xlab("m/z") +
        ylab("abundance") +
    theme_minimal() +
        theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
               axis.line=element_line(colour ="black",size=0.5, linetype="solid"),
               axis.ticks=element_line(colour ="black",size=0.5, linetype="solid"),
               panel.background = element_rect(fill="white", colour="white", size=0.5, linetype="solid"),
               legend.position="bottom",
               legend.background = element_rect(fill="white",size=0.25,linetype="blank"),
               legend.title = element_text("Plasma + Swab Extract Features Filtered"),
               legend.text=element_text(size=6)) +
    theme(aspect.ratio=1)
print(Original_Spectrum_plot)

KMF_spectrum <- as.data.frame(cbind(KMF_vec,as.numeric(t(output$Kendrickfiltered_MS1features_wgraph)[2:nrow(t(output$Kendrickfiltered_MS1features_wgraph)),2])))
label <- format(round(KMF_vec[,1], 2), nsmall = 4)
KMF_spectrum <- cbind(KMF_spectrum, label)
colnames(KMF_spectrum) <- c("mz","rt","abundance","mz_label")
head(KMF_spectrum)

KMF_spectrum_plot <- 
    ggplot()+
    geom_bar(data=KMF_spectrum,aes(x=mz, y=abundance), stat="identity",  
              size=0.5, alpha=1, color="black")+
    geom_text(data=KMF_spectrum,aes(x=mz, y=abundance,label=mz_label), cex=3, check_overlap = TRUE, vjust=0) +
   scale_x_continuous(breaks = seq(0,1000,100))+
        xlab("m/z") +
        ylab("abundance") +
    theme_minimal() +
        theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
               axis.line=element_line(colour ="black",size=0.5, linetype="solid"),
               axis.ticks=element_line(colour ="black",size=0.5, linetype="solid"),
               panel.background = element_rect(fill="white", colour="white", size=0.5, linetype="solid"),
               legend.position="bottom",
               legend.background = element_rect(fill="white",size=0.25,linetype="blank"),
               legend.title = element_text("Plasma + Swab Extract Features Filtered"),
               legend.text=element_text(size=6)) +
    theme(aspect.ratio=1)
print(KMF_spectrum_plot)

# MS1 Features (calculate MS1 Features that were removed rather than retained)
Filtered_matrix_KMF_wgraph_MS1_plot <- 
   t(data_matrix)[((vec[,1] %in% KMF_vec[,1]) & 
                  (vec[,2] %in% KMF_vec[,2])) != TRUE,]

KMF_filtered_features_spectrum <- as.data.frame(cbind(Filtered_matrix_KMF_wgraph_MS1_plot_vec, 
                                                      as.numeric(Filtered_matrix_KMF_wgraph_MS1_plot[2:nrow(Filtered_matrix_KMF_wgraph_MS1_plot),2])))
label <- format(round(Filtered_matrix_KMF_wgraph_MS1_plot_vec[,1], 2), nsmall = 4)
KMF_filtered_features_spectrum <- cbind(KMF_filtered_features_spectrum, label)
colnames(KMF_filtered_features_spectrum) <- c("mz","rt","abundance","mz_label")
head(KMF_filtered_features_spectrum)

KMF_filtered_features_spectrum_plot <- ggplot()+
    geom_bar(data=KMF_filtered_features_spectrum,aes(x=mz, y=abundance), stat="identity",size=0.5, alpha=1, color="black")+
    geom_text(data=KMF_filtered_features_spectrum,aes(x=mz, y=abundance,label=mz_label), cex=3, check_overlap = TRUE, vjust=0) +
    scale_x_continuous(breaks = seq(0,1000,100))+
    xlab("m/z") +
    ylab("abundance") +
    theme_minimal() +
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
               axis.line=element_line(colour ="black",size=0.5, linetype="solid"),
               axis.ticks=element_line(colour ="black",size=0.5, linetype="solid"),
               panel.background = element_rect(fill="white", colour="white", size=0.5, linetype="solid"),
               legend.position="bottom",
               legend.background = element_rect(fill="white",size=0.25,linetype="blank"),
               legend.title = element_text("Plasma + Swab Extract Features Filtered"),
               legend.text=element_text(size=6)) +
    theme(aspect.ratio=1)
print(KMF_filtered_features_spectrum_plot)

Kendrickplot <- 
    ggplot()+
    geom_point(data=as.data.frame(output$kendrickmassfilterinfo_filtered_wgraph), aes(x=nom_kend, y=msdefect), pch=16, 
               size=1, alpha=0.75, col="black")+
    geom_point(data=kendrickmassfilterinfo_diff,aes(x=nom_kend,y=msdefect), pch=16, 
               size=1, alpha=0.75,  col="#e41a1c")+
    scale_y_continuous(breaks = seq(-0.5,0.5,0.1))+
    scale_x_continuous(breaks = seq(0,1000,100))+
    xlab("Kendrick nominal m/z") +
    ylab("Kendrick mass defect (KMD)") +
    theme_minimal() +
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
               axis.line=element_line(colour ="black",size=0.5, linetype="solid"),
               axis.ticks=element_line(colour ="black",size=0.5, linetype="solid"),
               axis.title.x=element_text(size=6),
               axis.title.y=element_text(size=6),
               axis.text=element_text(size=6),
               panel.background = element_rect(fill="white", colour="white", size=0.5, linetype="solid"),
               legend.position="none",
               legend.background = element_rect(fill="white",size=0.25,linetype="blank"),
               legend.title = element_blank(),
               legend.text=element_text(size=6)) +
    theme(aspect.ratio=1)
    
MS1featureplot <- 
    ggplot()+
    geom_point(data=as.data.frame(vec),aes(x=rt, y=mz), pch=16, 
              size=1, alpha=0.75, color="black")+
    geom_point(data=Filtered_matrix_KMF_wgraph_MS1_plot_vec, aes(x=rt, y=mz), pch=16, 
               size=1, alpha=0.75, col="#e41a1c")+
    scale_y_continuous(breaks = seq(0,1000,100))+
    scale_x_continuous(breaks = seq(0,10,1))+
        xlab("retention time (min)") +
        ylab("m/z") +
    theme_minimal() +
        theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
               axis.line=element_line(colour ="black",size=0.5, linetype="solid"),
               axis.ticks=element_line(colour ="black",size=0.5, linetype="solid"),
               axis.title.x=element_text(size=6),
               axis.title.y=element_text(size=6),
               axis.text=element_text(size=6),
               panel.background = element_rect(fill="white", colour="white", size=0.5, linetype="solid"),
               legend.position="none",
               legend.background = element_rect(fill="white",size=0.25,linetype="blank"),
               legend.title = element_blank(),
               legend.text=element_text(size=6)) +
    theme(aspect.ratio=1)

Original_Spectrum_plot <- 
    ggplot()+
    geom_bar(data=original_spectrum,aes(x=mz, y=abundance), stat="identity",  
              size=0.5, alpha=1, color="black")+
    geom_text(data=original_spectrum,aes(x=mz, y=abundance,label=mz_label), size=2, check_overlap = TRUE, vjust=0) +
   scale_x_continuous(breaks = seq(0,1000,100))+
        xlab("m/z") +
        ylab("abundance") +
    theme_minimal() +
        theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
               axis.line=element_line(colour ="black",size=0.5, linetype="solid"),
               axis.ticks=element_line(colour ="black",size=0.5, linetype="solid"),
               axis.title.x=element_text(size=6),
               axis.title.y=element_text(size=6),
               axis.text=element_text(size=6),
               panel.background = element_rect(fill="white", colour="white", size=0.5, linetype="solid"),
               legend.position="none",
               legend.background = element_rect(fill="white",size=0.25,linetype="blank"),
               legend.title = element_blank(),
               legend.text=element_text(size=6)) 

KMF_spectrum_plot <- 
    ggplot()+
    geom_bar(data=KMF_spectrum,aes(x=mz, y=abundance), stat="identity",  
              size=0.5, alpha=1, color="black")+
    geom_text(data=KMF_spectrum,aes(x=mz, y=abundance,label=mz_label), size=2, check_overlap = TRUE, vjust=0) +
   scale_x_continuous(breaks = seq(0,1000,100))+
        xlab("m/z") +
        ylab("abundance") +
    theme_minimal() +
        theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
               axis.line=element_line(colour ="black",size=0.5, linetype="solid"),
               axis.ticks=element_line(colour ="black",size=0.5, linetype="solid"),
               axis.title.x=element_text(size=6),
               axis.title.y=element_text(size=6),
               axis.text=element_text(size=6),
               panel.background = element_rect(fill="white", colour="white", size=0.5, linetype="solid"),
               legend.position="none",
               legend.background = element_rect(fill="white",size=0.25,linetype="blank"),
               legend.title = element_blank(),
               legend.text=element_text(size=6)) 

KMF_filtered_features_spectrum_plot <- ggplot()+
    geom_bar(data=KMF_filtered_features_spectrum,aes(x=mz, y=abundance), stat="identity",size=0.5, alpha=1, color="red")+
    geom_text(data=KMF_filtered_features_spectrum,aes(x=mz, y=abundance,label=mz_label), size=2, check_overlap = TRUE, vjust=0) +
    scale_x_continuous(breaks = seq(0,1000,100))+
    xlab("m/z") +
    ylab("abundance") +
    theme_minimal() +
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
               axis.line=element_line(colour ="black",size=0.5, linetype="solid"),
               axis.ticks=element_line(colour ="black",size=0.5, linetype="solid"),
               axis.title.x=element_text(size=6),
               axis.title.y=element_text(size=6),
               axis.text=element_text(size=6),
               panel.background = element_rect(fill="white", colour="white", size=0.5, linetype="solid"),
               legend.position="none",
               legend.background = element_rect(fill="white",size=0.25,linetype="blank"),
               legend.title = element_blank(),
               legend.text=element_text(size=6))  

lay <- rbind(c(1,3),
             c(1,3),
             c(1,4),
             c(2,4),
             c(2,5),
             c(2,5))
plot <- grid.arrange(Kendrickplot,MS1featureplot,Original_Spectrum_plot,KMF_spectrum_plot,KMF_filtered_features_spectrum_plot, layout_matrix = lay)

ggsave("20180611_KMF_Figure_PEG400.pdf", plot, scale=1, width=7.1, units="in")

# write .csv table containing the KMF info matrix
write.csv(kendrickmassfilterinfo_diff, "20160611_PEG_KMF_info.csv", row.names=FALSE)

# write .csv table containing MS1 features remaining after Kendrick mass filtering
write.csv(output$Kendrickfiltered_MS1features_wgraph, "20160611_PEG_KMF_resultingdatamatrix.csv",row.names=FALSE)

# write .csv table of MS1 features removed from the data via the Kendrick mass filter
write.csv(Filtered_matrix_KMF_wgraph_MS1_plot_vec, "20160611_PEG_KMF_featuresfiltered.csv", row.names=FALSE)
