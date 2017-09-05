
library(Hmisc)
library("Rgraphviz")
library(graph)

# PERFORM: change path to .csv
#    Load in data matrix created in MZmine
#    Transpose matrix to get into the correct format for subsequent processes
#    Only keeps the (mass-to-charge;retention time) information

kendrick.mass.filter <- function(feature_matrix, polymer="PEG", const = 44/44.0262, const_mz= 44.0, 
			mass_defect_parameter = 0.005, retention_time=2.0, connection_filter=TRUE, 
			connections_min_number = 3, comp_list=NULL, sprint=FALSE) {


    data_matrix <- read.csv(feature_matrix, header=FALSE)
    mat1data <- data_matrix
    mat1data <- mat1data[1,]
    mat1data <- as.matrix(mat1data)
    mat1data <- t(mat1data)

    mat2data <- read.csv(feature_matrix, check.names=FALSE, stringsAsFactors=F)

    vec <- t(sapply(mat1data[-1,1], function(x) strsplit(as.character(x), ";")[[1]]))
    vec <- apply(vec, 2, as.numeric)


    cb <- combn(1:(length(vec[,1])), 2)

    kend <- vec[,1]*const
    msdefect <- round(kend)-kend

    mat <- cbind(vec, msdefect, kend, round(kend))
    colnames(mat) <- c("m/z","RT","msdefect","kend","nom_kend")

    v1 <- abs(mat[cb[1,],3]-mat[cb[2,],3])<= mass_defect_parameter
    v2 <- round(abs(mat[cb[1,],4]-mat[cb[2,],4]))==const_mz
    v4 <- abs(mat[cb[1,],2]-mat[cb[2,],2])<=retention_time

    v3 <- v1 & v2 & v4
    if(!sum(v3)) return(NULL)

    vpos <- unique(as.vector(cb[,v3]))
    matpos <- cbind(mat[vpos,], vpos)

    if(connection_filter) {
        gr <- ftM2graphNEL(as.matrix(t(cb[,v3])), edgemode = "undirected")
        conn <- connComp(gr)
    	connLength <- unlist(lapply(conn, length))
    
        fpos <- cbind(names(degree(gr)), degree(gr))
        sz <- sapply(fpos[,1], function(y) connLength[which(unlist(lapply(conn, function(x) sum(x==y)))==1)])
        fpos <- cbind(fpos, sz)
    	ffpos <- as.numeric(fpos[as.numeric(fpos[,3])> connections_min_number, 1])
    }
    mat2 <- cbind(mat[cb[1,v3], 3:4], mat[cb[2,v3], 3:4])
    mat2 <- cbind(mat2, mat2[,1]-mat2[,3], mat2[,2]-mat2[,4])
    colnames(mat2) <- c("msdefect","kend","msdefect","kend","massdefect_difference","mass_difference")

    if(sprint) {	
    print('Summary', quote = FALSE)
        print(paste('Polymer Filtered - ', polymer), quote = FALSE)
        print(paste('Parameter: mass defect - ', mass_defect_parameter), quote = FALSE)
        print(paste('Parameter: retention time - ', retention_time), quote = FALSE)
        print(paste('Parameter: node degree (minimum # of connection) - ', connections_min_number), quote = FALSE)
        print(paste('number of MS1 features before filtering - ', length(mat2data[1,])), quote = FALSE)
        print(paste('number of MS1 features after filtering (before reincluding based on minimum # of connections) - ', length(mat2data[2,-c(vpos+1)])), quote = FALSE)
   }
    
    if(length(vpos)) {
	    summ <- data.frame(polymer, mass_defect_parameter, retention_time, connections_min_number, 
			   features_before_filtering=length(mat2data[1,]),
			   features_after_filtering=length(mat2data[1,-c(vpos+1)])
		   ) 
    } else {
	    summ <- data.frame(polymer, mass_defect_parameter, retention_time, connections_min_number, 
			   features_before_filtering=length(mat2data[1,]),
			   features_after_filtering=length(mat2data[1,])
		   ) 

    }
    if(!is.null(comp_list)) {
	for(k in comp_list) {
		cname <- sub(".csv", "", strsplit(k, "/")[[1]][2]) 
		cname <- paste0("inputBfiltVS", cname) 
		ctab <- read.csv(k,  check.names=FALSE, stringsAsFactors=F, nrow=1)
    		if(length(vpos)) {
			summ <- cbind(summ, sum(colnames(mat2data)[-c(vpos+1)] %in% colnames(ctab)))
			colnames(summ)[ncol(summ)] <- cname
		} else {
			summ <- cbind(summ, sum(colnames(mat2data) %in% colnames(ctab)))
			colnames(summ)[ncol(summ)] <- cname
		}
	}
    }


    if(connection_filter) {
     if(sprint) print(paste('number of MS features after filtering (with reincluded based on minimum # of connections) - ', length(mat2data[2,-c(ffpos+1)])), quote = FALSE)
    if(length(ffpos)) {
	summ$features_after_filtering_wgraph <- length(mat2data[2,-c(ffpos+1)])
    } else {
	summ$features_after_filtering_wgraph <- length(mat2data[2,])
    }

    if(!is.null(comp_list)) {
	for(k in comp_list) {
		cname <- sub(".csv", "", strsplit(k, "/")[[1]][2]) 
		cname <- paste0("inputBfiltGraphVS", cname) 
		ctab <- read.csv(k,  check.names=FALSE, stringsAsFactors=F, nrow=1)
		#assign(cname,  sum(colnames(mat2data)[-c(vpos+1)] %in% colnames(ctab)))
		if(length(ffpos)) {
			summ <- cbind(summ, sum(colnames(mat2data)[-c(ffpos+1)] %in% colnames(ctab)))
			colnames(summ)[ncol(summ)] <- cname
		} else {
			summ <- cbind(summ, sum(colnames(mat2data) %in% colnames(ctab)))
			colnames(summ)[ncol(summ)] <- cname
		}
	}
    }

        return(list(data=mat2data[,-c(ffpos+1)], summary=summ))
    } else {
        return(list(data=mat2data[,-c(vpos+1)], summary=summ))
    }
}
