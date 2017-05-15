library(data.table)
library(dplyr)

## FUNCTION Transpose a datatable and fix column names retirn a dataframe
transp_DF <- function(df_to_transp){
  df_to_transp.T <- t(df_to_transp)
  colnames(df_to_transp.T) <- df_to_transp.T[1,]
  df_to_transp.T <- df_to_transp.T[2:nrow(df_to_transp.T), ]
  df_to_transp.T <- as.data.frame((df_to_transp.T))
}

## FUNCTION take a vector of permutated loci count the number of MLST types
## (called on columns in a matrix of loci permutations with apply)
get_mlst_groups_count <- function(mlst_genes){
      #cat(mlst_genes, ":\t")
      mlst_string <- paste0(sort(mlst_genes), sep="", collapse=":")  
     # mlst_string <- gsub("\\s", "", mlst_string)
      #cat(mlst_string, ":\t")
       mlst_genes <- as.vector(mlst_genes)
       ## select loci typing  info only for current permutation 
      temp_mlst_selection <- dplyr::filter(loci_info.T, gene %in% mlst_genes)  
      #View(temp_mlst_selection)
      
      ## transpose dataframe again and tidy up a bit
       temp_mlst_selectionT2  <-  transp_DF(temp_mlst_selection)
       temp_mlst_selectionT2  <-  tibble::rownames_to_column(temp_mlst_selectionT2, var = "isolate") %>% 
                                  select(-1) %>% 
                                  mutate_each(funs(as.numeric(as.character(.)))) 
      ## get current loci names                    
       temp_mlst_genenames <- names(temp_mlst_selectionT2)
       ## group by overall MLST types
       temp_mlst_summary   <- temp_mlst_selectionT2  %>% dplyr::group_by_(.dots = temp_mlst_genenames) %>% 
                                                         dplyr::tally()
      ## get count of MLST types for current permutation 
       temp_mlst_summary_count <- nrow(temp_mlst_summary)
       
       ## print MLST loci for MLST groups above a threshold  ## comment out if not wanted
     #  if(temp_mlst_summary_count >= 74){
  #     fileConn<-file("output.txt")

 #      writeLines(c("Hello","World"), fileConn)

  #     close(fileConn)
       cat(mlst_string, "\t", sep = "", file = "MLST_bs_out.txt", append = TRUE)
       cat(temp_mlst_summary_count, "\n", sep = "", file = "MLST_bs_out.txt", append = TRUE)

     #  }
       ## return count of MLST types for current permutation 
       ## with apply usage return value is appended to a vector
       return(temp_mlst_summary_count)
}

## Load / transpose / tidy up original data
## NB: columns extracted hard coded by index
loci_info <- fread("./loci_info.csv") %>% dplyr::select(isolate, 10:21)
## transpose data columns to rows
loci_info.T <- transp_DF(loci_info)
## NB: rows made numeric hard coded by row index
loci_info.T <- tibble::rownames_to_column(loci_info.T, var = "gene") %>% mutate_each(funs(as.numeric(as.character(.))), 2:108)




# loop through a count of 2:12 sampled loci
#for(i in 2:12){
i <- 7
  cat("loci_number = ", i, "\n")
  cat("loci\tcount", "\n", file = "MLST_bs_out.txt", append = FALSE)
  ## combn - generates a matrix of permutations for a given number of chosen loci
  mlst_permuationM <-  combn(loci_info.T$gene, i)
  ## apply loops through (the matrix of gene loci permutations, Columns, mlst counting function)
  


  temp_mlst_perm <- apply(mlst_permuationM, 2, get_mlst_groups_count)
   cat( "loci permutations = ", length(temp_mlst_perm),
        "\tmin mlst groups = ",  min(temp_mlst_perm),
        "\tmean mlst groups  = ", round(mean(temp_mlst_perm), 1),
        "\tmax  mlst groups = ",  max(temp_mlst_perm),
        "\n")
   

   temp_mlst_counts <- tibble(loci_count = i, min = min(temp_mlst_perm), 
                          mean = round(mean(temp_mlst_perm), 1),
                          max = max(temp_mlst_perm))
#  }

hist(temp_mlst_perm,
     main="7 Loci MLST Histogram", 
     xlab="MLST types", 
     border="blue", 
     col="green",
     las=1, 
     breaks=12)