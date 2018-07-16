
library(plotly)
library(htmlwidgets)
library(grid)
library(gridExtra)
library(viridis)
library(VennDiagram)
library(doMC)
library(IRdisplay)

source("kendrick.mass.filter_new_all_datamatrix_terminologyupdate.R")

kend_mass <- c(0.001, 0.0015, 0.002, 0.0025, 0.005, 0.0067, 0.0075, 0.01,
                   0.0125, 0.015, 0.02, 0.025, 0.033, 0.05, 0.067, 0.1)
kend_mass

rt <- c(seq(from = 1, to = 2, by = 0.25), seq(from = 0.1, to = 0.9, by = 0.1), seq(from=0.15, to=0.25, by = 0.05))
rt

conn_nun <- 1:5
conn_nun

all_comb <- expand.grid(kend_mass, conn_nun, rt)
head(all_comb)

plasma_fn = "all17K/tab17plasma.csv"
swabextract_fn = "all17K/tab17swabextract.csv"
peg_fn = "all17K/tab17PEG.csv"
plasmaspikedswab_fn = "all17K/tab17plasmaspikedswab.csv"

# Iteratively Perform KMF on Data (time parallel computation)
registerDoMC()
all_comb_ans3 <- list()
system.time(all_comb_ans3 <- foreach(i=1:nrow(all_comb)) %dopar%{
            composition.analysis(feature_matrix = plasmaspikedswab_fn, 
                                                       polymer = c("polyethylene_glycol_other_C2H4O1"),
                                                       mass_defect_parameter = all_comb[i,1],
                                                       retention_time = all_comb[i,3],
                                                       connections_min_number = all_comb[i,2],
                                                       comp_list = c(plasma_fn, 
                                                                     peg_fn, 
                                                                     swabextract_fn, 
                                                                     plasmaspikedswab_fn))})
# Origin of the extracted features to all parameter combinations
all_comb_tab <- lapply(all_comb_ans3, function(x) x$summary$polyethylene_glycol_other_C2H4O1)
all_comb_tab <- do.call(rbind, all_comb_tab)
all_comb_tab <- as.data.frame(all_comb_tab)
head(all_comb_tab)

write.csv(all_comb_tab,"20180611_Output_Table_KendrickMassFiter_TestingMultipleParameters.csv", row.names=FALSE)

# Load the matrices to perform comparison
plasma = read.csv(plasma_fn, check.names=FALSE)
swabextract = read.csv(swabextract_fn, check.names=FALSE)
peg = read.csv(peg_fn, check.names=FALSE)
plasmaspikedswab = read.csv(plasmaspikedswab_fn, check.names=FALSE)

# Origin of the extracted features to all parameter combinations
prop_polymer <- lapply(all_comb_ans3, function(x)
                                        c(length(x$data$polyethylene_glycol_other_C2H4O1),
                                          sum(colnames(plasmaspikedswab)[-1][x$data$polyethylene_glycol_other_C2H4O1] %in% colnames(plasma)),
                                          sum(colnames(plasmaspikedswab)[-1][x$data$polyethylene_glycol_other_C2H4O1] %in% colnames(swabextract)),
                                          sum(colnames(plasmaspikedswab)[-1][x$data$polyethylene_glycol_other_C2H4O1] %in% colnames(peg)),
                                          sum(colnames(plasmaspikedswab)[-1][x$data$polyethylene_glycol_other_C2H4O1] %in% colnames(swabextract) & colnames(plasmaspikedswab)[-1][x$data$polyethylene_glycol_other_C2H4O1] %in% colnames(peg))
                                          )
)
    prop_polymer <- do.call(rbind, prop_polymer)
    colnames(prop_polymer) <- c('Plasma_spiked_SwabExtract', 'plasma', 'swabextract', 'peg', 'plasma_spiked_peg')
    prop_polymer <- as.data.frame(prop_polymer)
    head(prop_polymer)

# Compute Scoring (plasma/PEG features) for plot
prop_ratio <- prop_polymer[,"plasma"] / prop_polymer[,"peg"]
PEG_optimization <- cbind(prop_polymer,prop_ratio)
    head(PEG_optimization)
write.csv(PEG_optimization, "20180611_Output_Table_PEGoptimization.csv", row.names=FALSE)

static_plot <- ggplot(prop_polymer, aes(x=prop_ratio, y=peg, colour=plasma_spiked_peg))+
    geom_point()+    
    scale_y_continuous(breaks=seq(0,100,10)) +
        xlab("Plasma Features Filtered / PEG400 Features Filtered") +
        ylab("Number of PEG400 Features Filtered") +
    scale_color_viridis() +
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
print(static_plot)
ggsave("ParameterExplorationPlot.pdf")

p <- plot_ly(type="scatter",mode="markers", data = prop_polymer, 
                              x = ~prop_ratio, y = ~peg, color=~plasma_spiked_peg)
saveWidget(p, "plot_ly.html")
display_html("<iframe width='100%' height='550' src='plot_ly.html'></iframe>")

prop_polymer[which(round(prop_ratio, 1)==0.4 & prop_polymer[,'peg']==82),]

all_comb_tab[which(round(prop_ratio, 1)==0.4 & prop_polymer[,'peg']==82),]

nun_peaks_extracted = cbind(prop_polymer[all_comb_tab[,'retention_time']== 0.8,],
                            all_comb_tab[all_comb_tab[,'retention_time']== 0.8,]
                            )
colnames(nun_peaks_extracted)[1] <- 'feat_extracted'
nun_peaks_extracted$connections_min_number <- as.factor(nun_peaks_extracted$connections_min_number)
nun_peaks_extracted <- cbind(nun_peaks_extracted$plasma/nun_peaks_extracted$peg,nun_peaks_extracted)
colnames(nun_peaks_extracted)[1] <- 'prop_ratio'
head(nun_peaks_extracted)

write.csv(nun_peaks_extracted,"20180611_Output_Table_KMFEvaluatedwithPEGparameters.csv", row.names=FALSE)

a <- ggplot(nun_peaks_extracted, aes(x=mass_defect_parameter, y=peg, colour=connections_min_number)) + 
    scale_y_continuous(breaks=seq(0,100,10)) +
        xlab("Kendrick mass defect") +
        ylab("# PEG400 MS1 Features Filtered") +
    geom_line() +          
        geom_point(aes(fill=connections_min_number), alpha=0.8, colour="black", pch=21, size=2.5) +
    
    theme_minimal() +
        theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
               axis.line=element_line(colour ="black",size=0.5, linetype="solid"),
               axis.ticks=element_line(colour ="black",size=0.5, linetype="solid"),
               panel.background = element_rect(fill="white", colour="white", size=0.5, linetype="solid"),
               legend.position="none",
               legend.background = element_rect(fill="white",size=0.25,linetype="blank"),
               legend.title = element_blank(),
               legend.text=element_text(size=8))


b <- ggplot(nun_peaks_extracted, aes(x=mass_defect_parameter, y=plasma, colour=connections_min_number)) +   
    scale_y_continuous(breaks=seq(0,100,10)) +
        xlab("Kendrick mass defect") +
        ylab("# Plasma MS1 Features Filtered") +
    geom_line() +          
        geom_point(aes(fill=connections_min_number), alpha=0.8, colour="black", pch=24, size=2.5) +
    
    theme_minimal() +
        theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
               axis.line=element_line(colour ="black",size=0.5, linetype="solid"),
               axis.ticks=element_line(colour ="black",size=0.5, linetype="solid"),
               panel.background = element_rect(fill="white", colour="white", size=0.5, linetype="solid"),
               legend.position="none",
               legend.background = element_rect(fill="white",size=0.25,linetype="blank"),
               legend.title = element_blank(),
               legend.text=element_text(size=8))

c <- ggplot(nun_peaks_extracted, aes(x=mass_defect_parameter, y=plasma_spiked_peg, colour=connections_min_number)) + 
    scale_y_continuous(breaks=seq(0,200,10)) +
        xlab("Kendrick mass defect") +
        ylab("# Plasma + PEG400 MS1 Features Filtered") +
    geom_line() +          
        geom_point(aes(fill=connections_min_number), alpha=0.8, colour="black", pch=22, size=2.5) +
    
    theme_minimal() +
        theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
               axis.line=element_line(colour ="black",size=0.5, linetype="solid"),
               axis.ticks=element_line(colour ="black",size=0.5, linetype="solid"),
               panel.background = element_rect(fill="white", colour="white", size=0.5, linetype="solid"),
               legend.position="none",
               legend.background = element_rect(fill="white",size=0.25,linetype="blank"),
               legend.title = element_blank(),
               legend.text=element_text(size=8))

d <- ggplot(nun_peaks_extracted, aes(x=mass_defect_parameter, y=feat_extracted, colour=connections_min_number)) +  
    scale_y_continuous(breaks=seq(0,1200,100)) +
        xlab("Kendrick mass defect") +
        ylab("# Plasma + Swab Extract MS1 Features Filtered") +
    geom_line() +          
        geom_point(aes(fill=connections_min_number), alpha=0.8, colour="black", pch=23, size=2.5) +
    
    theme_minimal() +
        theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
               axis.line=element_line(colour ="black",size=0.5, linetype="solid"),
               axis.ticks=element_line(colour ="black",size=0.5, linetype="solid"),
               panel.background = element_rect(fill="white", colour="white", size=0.5, linetype="solid"),
               legend.position="none",
               legend.background = element_rect(fill="white",size=0.25,linetype="blank"),
               legend.title = element_blank(),
               legend.text=element_text(size=8))

grid.arrange(a,b,c,d,ncol=2, nrow=2)
plot_all <- arrangeGrob(a,b,c,d,ncol=2, nrow=2)
ggsave("20180611_KMF_Figure_KMDandNOS.pdf", plot_all, scale=1.5)

nun_peaks_extracted_RT = cbind(prop_polymer[all_comb_tab[,'mass_defect_parameter']==0.01,],
                            all_comb_tab[all_comb_tab[,'mass_defect_parameter']==0.01,]
                            )
colnames(nun_peaks_extracted_RT)[1] <- 'feat_extracted'
nun_peaks_extracted_RT$connections_min_number <- as.factor(nun_peaks_extracted_RT$connections_min_number)
nun_peaks_extracted_RT <- cbind(nun_peaks_extracted_RT$plasma/nun_peaks_extracted_RT$peg,nun_peaks_extracted_RT)
colnames(nun_peaks_extracted_RT)[1] <- 'prop_ratio'
head(nun_peaks_extracted_RT)

a <- ggplot(nun_peaks_extracted_RT, aes(x=retention_time, y=peg, colour=connections_min_number)) +  
    scale_y_continuous(breaks=seq(0,100,10)) +
        xlab("Retention time (min)") +
        ylab("# PEG400 MS1 Features Filtered") +
    geom_line() +          
        geom_point(aes(fill=connections_min_number), alpha=0.8, colour="black", pch=21, size=2.5) +
    
    theme_minimal() +
        theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
               axis.line=element_line(colour ="black",size=0.5, linetype="solid"),
               axis.ticks=element_line(colour ="black",size=0.5, linetype="solid"),
               panel.background = element_rect(fill="white", colour="white", size=0.5, linetype="solid"),
               legend.position="none",
               legend.background = element_rect(fill="white",size=0.25,linetype="blank"),
               legend.title = element_blank(),
               legend.text=element_text(size=8))


b <- ggplot(nun_peaks_extracted_RT, aes(x=retention_time, y=plasma, colour=connections_min_number)) +
    scale_y_continuous(breaks=seq(0,100,10)) +
        xlab("Retention time (min)") +
        ylab("# Plasma MS1 Features Filtered") +
    geom_line() +          
        geom_point(aes(fill=connections_min_number), alpha=0.8, colour="black", pch=24, size=2.5) +
    
    theme_minimal() +
        theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
               axis.line=element_line(colour ="black",size=0.5, linetype="solid"),
               axis.ticks=element_line(colour ="black",size=0.5, linetype="solid"),
               panel.background = element_rect(fill="white", colour="white", size=0.5, linetype="solid"),
               legend.position="none",
               legend.background = element_rect(fill="white",size=0.25,linetype="blank"),
               legend.title = element_blank(),
               legend.text=element_text(size=8))

c <- ggplot(nun_peaks_extracted_RT, aes(x=retention_time, y=plasma_spiked_peg, colour=connections_min_number)) +  
    scale_y_continuous(breaks=seq(0,200,10)) +
        xlab("Retention time (min)") +
        ylab("# Plasma + PEG400 MS1 Features Filtered") +
    geom_line() +          
        geom_point(aes(fill=connections_min_number), alpha=0.8, colour="black", pch=22, size=2.5) +
    
    theme_minimal() +
        theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
               axis.line=element_line(colour ="black",size=0.5, linetype="solid"),
               axis.ticks=element_line(colour ="black",size=0.5, linetype="solid"),
               panel.background = element_rect(fill="white", colour="white", size=0.5, linetype="solid"),
               legend.position="none",
               legend.background = element_rect(fill="white",size=0.25,linetype="blank"),
               legend.title = element_blank(),
               legend.text=element_text(size=8))

d <- ggplot(nun_peaks_extracted_RT, aes(x=retention_time, y=feat_extracted, colour=connections_min_number)) +
    scale_y_continuous(breaks=seq(0,1200,100)) +
        xlab("Retention time (min)") +
        ylab("# Plasma + Swab Extract MS1 Features Filtered") +
    geom_line() +          
        geom_point(aes(fill=connections_min_number), alpha=0.8, colour="black", pch=23, size=2.5) +
    
    theme_minimal() +
        theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
               axis.line=element_line(colour ="black",size=0.5, linetype="solid"),
               axis.ticks=element_line(colour ="black",size=0.5, linetype="solid"),
               panel.background = element_rect(fill="white", colour="white", size=0.5, linetype="solid"),
               legend.position="none",
               legend.background = element_rect(fill="white",size=0.25,linetype="blank"),
               legend.title = element_blank(),
               legend.text=element_text(size=8))

grid.arrange(a,b,c,d,ncol=2, nrow=2)
plot_all <- arrangeGrob(a,b,c,d,ncol=2, nrow=2)
ggsave("20180611_KMF_Figure_RTandNOS.pdf", plot_all, scale=1.5)

# RATIO (Optimation Plot) w/ Retention Time Parameter of 0.2 min
a <- ggplot(nun_peaks_extracted, aes(x=mass_defect_parameter, y=prop_ratio, colour=connections_min_number)) +
    scale_x_continuous(breaks=nun_peaks_extracted$mass_defect_parameter)+
        xlab("Kendrick mass defect") +
        ylab("# Plasma MS1 Features Filtered / # PEG400 MS1 Features Filtered") +
    geom_line() +          
        geom_point(aes(fill=connections_min_number), alpha=1, size=2) +
    
    theme_minimal() +
        theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
               axis.line=element_line(colour ="black",size=0.5, linetype="solid"),
               axis.ticks=element_line(colour ="black",size=0.5, linetype="solid"),
               panel.background = element_rect(fill="white", colour="white", size=0.5, linetype="solid"),
               legend.position="bottom",
                  axis.text.x = element_text(angle = 45, hjust = 1),
               legend.background = element_rect(fill="white",size=0.25,linetype="blank"),
               legend.title = element_blank(),
               legend.text=element_text(size=8))
print(a)
ggsave("20180611_KMF_Figure_KMFScoring_KMD.pdf", a, scale=1.5)

# RATIO (Optimation Plot) w/ Retention Time Parameter of 0.2 min
a <- ggplot(nun_peaks_extracted_RT, aes(x=retention_time, y=prop_ratio, colour=connections_min_number)) +
    #scale_y_continuous(breaks=seq(0,5,0.1)) +
    scale_x_continuous(breaks=nun_peaks_extracted_RT$retention_time)+
        xlab("Retention time (min)") +
        ylab("# Plasma MS1 Features Filtered / # PEG400 MS1 Features Filtered") +
    geom_line() +          
        geom_point(aes(fill=connections_min_number), alpha=1, size=2) +
    
    theme_minimal() +
        theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
               axis.line=element_line(colour ="black",size=0.5, linetype="solid"),
               axis.ticks=element_line(colour ="black",size=0.5, linetype="solid"),
               panel.background = element_rect(fill="white", colour="white", size=0.5, linetype="solid"),
               legend.position="bottom",
               legend.background = element_rect(fill="white",size=0.25,linetype="blank"),
               legend.title = element_blank(),
               legend.text=element_text(size=8))
print(a)
ggsave("20180611_KMF_Figure_KMFScoring_RT.pdf", a, scale=1.5)
