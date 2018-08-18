#Kendrick Mass Filter and Composition Analysis
    #Authors: Ricardo R. da Silva and Alan K. Jarmusch
    #Verion: 1.0 (Prior to Submission of Manuscript)
    #Date of Last Revision: 03-22-2018 

library(Hmisc)
library("Rgraphviz")
library(graph)

        #STILL REVELVANT??? <
            # PERFORM: change path to .csv
                #    Load in data matrix created in MZmine
            #    Transpose matrix to get into the correct format for subsequent processes
                #    Only keeps the (mass-to-charge;retention time) information
        # >

#####################################
### Kendrick Mass Filter Function ###                                                                         
#####################################

isPeak <- function(x) {
    mx <- max(x)
    peaks <- c()
    for(i in 2:(length(x)-1)) {
        if (x[i] < x[i+1] & x[i]!=mx){
           if(x[i] > x[i-1] & x[i] < x[i+1]) {
               peaks <- c(peaks, i-1, i, i+1)
               mx <- x[i]
           }
        } else {
           if((x[i] < x[i-1] & x[i] > x[i+1]) | x[i]==mx) {
               peaks <- c(peaks, i-1, i, i+1)
           }
        }
    }
    return(sum(1:length(x) %in% unique(peaks))==length(x))
}


getPeakShape <- function(vec, dimension='mz', minpeaks=3, mxtol=0.6) {
    if(nrow(vec)<minpeaks) {
        return(NULL)
    } 
    od <- order(vec[,dimension])
    x <- vec[od, 'int']
    peakmat <- lapply(minpeaks:length(x), function(y) blockparts(rep(1,length(x)),y))
    peakmat <- do.call(cbind, peakmat)

    peaklist <- list()
    idx <- 1

    for(k in 1:ncol(peakmat)) {
        deltaint <- c()
        x <- vec[od, 'int'][peakmat[,k]==1]
        for(i in 1:(length(x)-1)) {
            deltaint <- c(deltaint, abs(x[i]-x[i+1])/max(c(x[i],x[i+1])))
        }
        deltaint[is.na(deltaint)] <- 1
        if (!sum(deltaint>mxtol)) {
            peaklist[[idx]] <- peakmat[,k] 
            idx <- idx+1
        }
    }
    if(!length(peaklist)) {
        return(NULL)
    }
    ispeak <-  unlist(lapply(peaklist, function(x) isPeak(vec[od, 'int'][x==1])))
    peaklist <- peaklist[ispeak]
    selid <- od[peaklist[[which.max(unlist(lapply(peaklist, sum)))]]==1]

    return(list(peaklist=peaklist, selid=selid))
}



Kendrick.mass.filter <- function(data_matrix, 
                       # vec now includes a third column with intensity
                       vec, 
                       polymer="polyethylene_glycol_other_C2H4O1",
		       repeat_unit_mass=44.02622,
                       KMD=0.01, 
                       RT=0.2,
                       connection_filter=TRUE, 
                       NOS = 2,
		       filterShape=FALSE,
		       fractionBase=1) 

{
    mass_defect_parameter <- KMD
    retention_time <- RT
    connections_min_number <- NOS - 1
 
   if(!is.null(polymer)) {
        polymer_list <- 
            list(
                'alkane_other_CH2' = list('const' = 14.01565, 'const_mz' = 14.0),
                'oxidation' = list('const' = 15.99492, 'const_mz' = 16.0),
                'water_cluster' = list('const' = 18.01057, 'const_mz' =18.0),
                'alkane_C2H4' = list('const' = 28.03130, 'const_mz' =28.0),
                'methanol_cluster' = list('const' = 32.02622, 'const_mz' =32.0),
                'acetonitrile_cluster' = list('const' = 41.02655, 'const_mz' =41.0),
                'propylation_other_C3H6' = list('const' = 42.04695, 'const_mz' =42.0),
                'polyethylene_glycol_other_C2H4O1' = list('const' = 44.02622, 'const_mz' =44.0),
                'perfluoro_CF2' = list('const' = 49.99681, 'const_mz' =50.0),
                'ammoniumchloride_cluster' = list('const' = 53.00323, 'const_mz' =53.0),
                'butylation_other_C4H8' = list('const' = 56.06260, 'const_mz' =56.0),
                'sodiumchloride_cluster' = list('const' = 57.95862, 'const_mz' =58.0),
                'polypropylene_glycol_other_C3H6O1' = list('const' = 58.04187, 'const_mz' =58.0),
                'ammoniumformate_cluster' = list('const' = 63.03203, 'const_mz' =63.0),
                'sodiumformate_cluster' = list('const' = 67.98742, 'const_mz' =68.0),
                'potassiumchloride_cluster' = list('const' = 73.93256, 'const_mz' =74.0),
                'polysiloxane' = list('const' = 74.01879, 'const_mz' =74.0),
                'sodiumacetate_cluster' = list('const' = 82.00307, 'const_mz' =82.0)
                )
        
        const <- polymer_list[[which(names(polymer_list)==polymer)]]$const
    
        if(fractionBase==1) {
            const_mz <- round(const) 
            const <- round(const)/const
            kend0 <- vec[,1]*const
        } else {
            const_mz <- round(const) 
            #const <- round(const)/const
            #kend0 <- vec[,1]*const
            #const_mz <- round(const/fractionBase) 
            const <- round((const/fractionBase))/(const/fractionBase)
            #mass_defect_parameter <- mass_defect_parameter*const
        }
        #const_mz <- round(const/fractionBase) 
        #const <- round(const/fractionBase)/(const/fractionBase) 
        #const_mz <- polymer_list[[which(names(polymer_list)==polymer)]]$const_mz
    } else {
        const <- repeat_unit_mass
    
        if(fractionBase==1) {
            const_mz <- round(const) 
            const <- round(const)/const
            kend0 <- vec[,1]*const
        } else {
            const_mz <- round(const) 
            #const <- round(const)/const
            #kend0 <- vec[,1]*const
            #const_mz <- round(const/fractionBase) 
            const <- round((const/fractionBase))/(const/fractionBase)
            #mass_defect_parameter <- mass_defect_parameter*const
        }
        #const_mz <- round(const/fractionBase) 
        #const <- round(const/fractionBase)/(const/fractionBase) 
        #const_mz <- polymer_list[[which(names(polymer_list)==polymer)]]$const_mz

    }
    
    cb <- combn(1:(length(vec[,1])), 2)

    kend <- vec[,1]*const
    msdefect <- round(kend)-kend

    if(ncol(vec)==2) {
        mat <- cbind(vec, msdefect, kend, round(kend))
    } else {
        mat <- cbind(vec[,-3], msdefect, kend, round(kend))
    }
    colnames(mat) <- c("m/z","RT","msdefect","kend","nom_kend")
    
    v1 <- abs(mat[cb[1,],3]-mat[cb[2,],3])<= mass_defect_parameter
    if(!sum(v1)) return(NULL)
    cb <- matrix(cb[,v1], nrow=2)

    v2 <- round(abs(mat[cb[1,],4]-mat[cb[2,],4]))==const_mz
    #v2 <- round(abs(kend0[cb[1,]]-kend0[cb[2,]]))==const_mz
    if(!sum(v2)) return(NULL)
    cb <- matrix(cb[,v2], nrow=2)

    v4 <- abs(mat[cb[1,],2]-mat[cb[2,],2])<=retention_time
    if(!sum(v4)) return(NULL)
    cb2 <- matrix(cb[,v4], nrow=2)

    vpos <- unique(as.vector(cb2))
    if(connection_filter) {
        if(filterShape) {
            # create a matrix with intensity
            # and decide what to export
            gr <- ftM2graphNEL(as.matrix(t(cb)), edgemode = "undirected")
            conn <- connComp(gr)
            filt <- list()
            for (j in 1:length(conn)){
                ans <- getPeakShape(vec[as.numeric(conn[[j]]),])
                if(!is.null(ans)) {
                    filt[[j]] <- ans$selid
                } else {
                    filt[[j]] <- 0 
                }
    
            }	
	    return(list(conn=conn, filt=filt))
        }
        gr <- ftM2graphNEL(as.matrix(t(cb2)), edgemode = "undirected")
        conn <- connComp(gr)

        connLength <- unlist(lapply(conn, length))
    
        fpos <- cbind(names(degree(gr)), degree(gr))
        sz <- sapply(fpos[,1], function(y) connLength[which(unlist(lapply(conn, function(x) sum(x==y)))==1)])
        fpos <- cbind(fpos, sz)
        ffpos <- as.numeric(fpos[as.numeric(fpos[,3])> connections_min_number, 1])
    } else {
    ffpos <- NULL
    }
    mat2 <- cbind(mat[cb2[1,], 3:4], mat[cb2[2,], 3:4])
    mat2 <- cbind(mat2, mat2[,1]-mat2[,3], mat2[,2]-mat2[,4])
    colnames(mat2) <- c("msdefect","kend","msdefect","kend","massdefect_difference","mass_difference")
    
    kendrickmassfilterinfo_filtered_wograph <- mat[-c(vpos+1),]
    kendrickmassfilterinfo_filtered_wgraph <- mat[-c(ffpos+1),]                                    
                                                                          
    Kendrickfiltered_MS1features_wograph <- data_matrix[,-c(vpos+1)]
    Kendrickfiltered_MS1features_wgraph <- data_matrix[,-c(ffpos+1)]                                                             
                                                                          
    return(list(wograph=vpos, kendrickmassfilterinfo_original=mat, wgraph=ffpos,
                kendrickmassfilterinfo_filtered_wograph=kendrickmassfilterinfo_filtered_wograph,
                kendrickmassfilterinfo_filtered_wgraph=kendrickmassfilterinfo_filtered_wgraph,
                Kendrickfiltered_MS1features_wograph=Kendrickfiltered_MS1features_wograph,
                Kendrickfiltered_MS1features_wgraph=Kendrickfiltered_MS1features_wgraph))
}

#####################################
### necessary for compositional analysis ###                                                                         
##################################### 
                                                                          
get.series <- function(vec, const=44/44.02622, const_mz=44.0, 
                      mass_defect_parameter=0.005, retention_time=1.0, 
		      connection_filter=TRUE, connections_min_number = 3) 

{
    
    cb <- combn(1:(length(vec[,1])), 2)

    kend <- vec[,1]*const
    msdefect <- round(kend)-kend

    mat <- cbind(vec, msdefect, kend, round(kend))
    colnames(mat) <- c("m/z","RT","msdefect","kend","nom_kend")
    
    v1 <- abs(mat[cb[1,],3]-mat[cb[2,],3])<= mass_defect_parameter
    if(!sum(v1)) return(NULL)
    cb <- matrix(cb[,v1], nrow=2)

    v2 <- round(abs(mat[cb[1,],4]-mat[cb[2,],4]))==const_mz
    if(!sum(v2)) return(NULL)
    cb <- matrix(cb[,v2], nrow=2)

    v4 <- abs(mat[cb[1,],2]-mat[cb[2,],2])<=retention_time
    if(!sum(v4)) return(NULL)
    cb2 <- matrix(cb[,v4], nrow=2)

    #v3 <- v1 & v2 & v4

    #cb2 <- cb[,v3]
    vpos <- unique(as.vector(cb2))
    if(connection_filter) {
        gr <- ftM2graphNEL(as.matrix(t(cb2)), edgemode = "undirected")
        conn <- connComp(gr)
    	connLength <- unlist(lapply(conn, length))
    
        fpos <- cbind(names(degree(gr)), degree(gr))
        sz <- sapply(fpos[,1], function(y) connLength[which(unlist(lapply(conn, function(x) sum(x==y)))==1)])
        fpos <- cbind(fpos, sz)
    	ffpos <- as.numeric(fpos[as.numeric(fpos[,3])> connections_min_number, 1])
    } else {
	ffpos <- NULL
    }
   ### Questions here... need to output mat2 with the kendrick mass and kendrick mass defect for plotting.
    mat2 <- cbind(mat[cb2[1,], 3:4], mat[cb2[2,], 3:4])
    mat2 <- cbind(mat2, mat2[,1]-mat2[,3], mat2[,2]-mat2[,4])
    colnames(mat2) <- c("msdefect","kend","msdefect","kend","massdefect_difference","mass_difference")
    
    kendrickmassfilterinfo_filtered_wograph <- mat[-c(vpos+1),]
    kendrickmassfilterinfo_filtered_wgraph <- mat[-c(ffpos+1),]                                    

    return(list(wograph=vpos, kendrickmassfilterinfo_original=mat, wgraph=ffpos,
                kendrickmassfilterinfo_filtered_wograph=kendrickmassfilterinfo_filtered_wograph,
                kendrickmassfilterinfo_filtered_wgraph=kendrickmassfilterinfo_filtered_wgraph))
}
                                                                          
                                                                          
                                                                          
                                                                          
#####################################
### Composition Analysis Function ###                                                                         
#####################################
                                                                          
composition.analysis  <- function(
              feature_matrix,
              polymer=c('alkane_other_CH2',
              'oxidation',
              'water_cluster',
              'alkane_C2H4',
              'methanol_cluster',
              'acetonitrile_cluster',
              'propylation_other_C3H6',
              'polyethylene_glycol_other_C2H4O1',
              'perfluoro_CF2',
              'ammoniumchloride_cluster',
              'butylation_other_C4H8',
              'sodiumchloride_cluster',
              'polypropylene_glycol_other_C3H6O1',
              'ammoniumformate_cluster',
              'potassiumchloride_cluster',
              'polysiloxane',
              'sodiumacetate_cluster'), 
              mass_defect_parameter = 0.01, 
              retention_time=0.2, 
              connection_filter=TRUE,
              connections_min_number = 1,
              comp_list=NULL,
              sprint=FALSE) 

{
    polymer_list <- 
        list(
            'alkane_other_CH2' = list('const' = 14/14.01565, 'const_mz' = 14.0),
            'oxidation' = list('const' = 16/15.99492, 'const_mz' = 16.0),
            'water_cluster' = list('const' = 18/18.01057, 'const_mz' =18.0),
            'alkane_C2H4' = list('const' = 28/28.03130, 'const_mz' =28.0),
            'methanol_cluster' = list('const' = 32/32.02622, 'const_mz' =32.0),
            'acetonitrile_cluster' = list('const' = 41/41.02655, 'const_mz' =41.0),
            'propylation_other_C3H6' = list('const' = 42/42.04695, 'const_mz' =42.0),
            'polyethylene_glycol_other_C2H4O1' = list('const' = 44/44.02622, 'const_mz' =44.0),
            'perfluoro_CF2' = list('const' = 50/49.99681, 'const_mz' =50.0),
            'ammoniumchloride_cluster' = list('const' = 53/53.00323, 'const_mz' =53.0),
            'butylation_other_C4H8' = list('const' = 56/56.06260, 'const_mz' =56.0),
            'sodiumchloride_cluster' = list('const' = 58/57.95862, 'const_mz' =58.0),
            'polypropylene_glycol_other_C3H6O1' = list('const' = 58/58.04187, 'const_mz' =58.0),
            'ammoniumformate_cluster' = list('const' = 63/63.03203, 'const_mz' =63.0),
            'sodiumformate_cluster' = list('const' = 68/67.98742, 'const_mz' =68.0),
            'potassiumchloride_cluster' = list('const' = 74/73.93256, 'const_mz' =74.0),
            'polysiloxane' = list('const' = 74/74.01879, 'const_mz' =74.0),
            'sodiumacetate_cluster' = list('const' = 82/82.00307, 'const_mz' =82.0)
            )
    
    data_matrix <- read.csv(feature_matrix, header=TRUE, check.names=FALSE)
    mat1data <- colnames(data_matrix)
    mat1data <- mat1data[grep(';', mat1data)] 
    othercnames <- setdiff(colnames(data_matrix), mat1data)

    vec <- t(sapply(mat1data, function(x) strsplit(as.character(x), ";")[[1]]))
    vec <- apply(vec, 2, as.numeric)

    summList <- list()
    posList <- list()
    dataList <- list()
    for(p in polymer) {
        const <- polymer_list[[which(names(polymer_list)==p)]]$const
        const_mz <- polymer_list[[which(names(polymer_list)==p)]]$const_mz

	plist <- get.series(vec, const, const_mz, 
			   mass_defect_parameter, retention_time, 
			   connection_filter, connections_min_number)

	if(!is.null(plist)) {
	    vpos <- plist[[1]]
	    mat2 <- plist[[2]]
	    if(connection_filter) {
		ffpos <- plist[[3]]
	    }
	} else {
	    return(NULL)
	}

	if(sprint) {	
	print('Summary', quote = FALSE)
	    print(paste('Polymer Filtered - ', p), quote = FALSE)
	    print(paste('Parameter: mass defect - ', mass_defect_parameter), quote = FALSE)
	    print(paste('Parameter: retention time - ', retention_time), quote = FALSE)
	    print(paste('Parameter: node degree (minimum # of connection) - ', connections_min_number), quote = FALSE)
	    print(paste('number of MS1 features before filtering - ', length(mat1data)), quote = FALSE)
	    print(paste('number of MS1 features after filtering (before reincluding based on minimum # of connections) - ', length(mat1data)-length(vpos)), quote = FALSE)
        }
	
	if(length(vpos)) {
		summ <- data.frame(p, mass_defect_parameter, retention_time, connections_min_number, 
			       features_before_filtering=length(mat1data),
			       features_after_filtering=length(mat1data)-length(vpos)
		       ) 
	} else {
		summ <- data.frame(p, mass_defect_parameter, retention_time, connections_min_number, 
			       features_before_filtering=length(mat1data),
			       features_after_filtering=length(mat1data)
		       ) 

	}
	if(!is.null(comp_list)) {
	    for(k in comp_list) {
		    cname <- sub(".csv", "", strsplit(k, "/")[[1]][2]) 
                    # CHANGE COLUMN NAME HERE INSIDE paste FUNCTION
		    cname <- paste0("inputAfterFilterVS", cname) 
		    ctab <- read.csv(k,  check.names=FALSE, stringsAsFactors=F, nrow=1)
		    if(length(vpos)) {
			    summ <- cbind(summ, sum(mat1data[-vpos] %in% colnames(ctab)))
			    colnames(summ)[ncol(summ)] <- cname
		    } else {
			    summ <- cbind(summ, sum(mat1data %in% colnames(ctab)))
			    colnames(summ)[ncol(summ)] <- cname
		    }
	    }
	}


	if(connection_filter) {
	    if(sprint) print(paste('number of MS features after filtering (with reincluded based on minimum # of connections) - ', length(mat1data)-length(ffpos)), quote = FALSE)
	    if(length(ffpos)) {
		summ$features_after_filtering_wgraph <- length(mat1data)-length(ffpos)
	    } else {
		summ$features_after_filtering_wgraph <- length(mat1data)
	    }

	    if(!is.null(comp_list)) {
		for(k in comp_list) {
			cname <- sub(".csv", "", strsplit(k, "/")[[1]][2]) 
                        # CHANGE COLUMN NAME HERE INSIDE paste FUNCTION
			cname <- paste0("inputAfterFilterAndGraphVS", cname) 
			ctab <- read.csv(k,  check.names=FALSE, stringsAsFactors=F, nrow=1)
			#assign(cname,  sum(colnames(mat2data)[-c(vpos+1)] %in% colnames(ctab)))
			if(length(ffpos)) {
				summ <- cbind(summ, sum(mat1data[-ffpos] %in% colnames(ctab)))
				colnames(summ)[ncol(summ)] <- cname
			} else {
				summ <- cbind(summ, sum(mat1data %in% colnames(ctab)))
				colnames(summ)[ncol(summ)] <- cname
			}
		}
	    }

	    #return(list(data=ffpos, summary=summ))
	    posList[[p]] <- ffpos
	    summList[[p]] <- summ
            if(length(ffpos)) {
                dataList[[p]] <- data_matrix[,c(othercnames, mat1data[-ffpos])]
            } else {
                dataList[[p]] <- data_matrix[,c(othercnames, mat1data[-ffpos])]
            }
	} else {
	    #return(list(data=vpos, summary=summ))
	    posList[[p]] <- vpos
	    summList[[p]] <- summ
            if(length(vpos)) {
                dataList[[p]] <- data_matrix[,c(othercnames, mat1data[-vpos])]
            } else {
                dataList[[p]] <- data_matrix[,c(othercnames, mat1data[-vpos])]
            }
	}
    }

    return(list(data=posList, summary=summList, dataList=dataList))

}
