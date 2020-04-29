#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 
if (length(args)<=1) {
  # contrasts
  args[2] = "pairwise"
}
if (length(args)<=2) {
  # default control condition
  args[3] = ""
}
if (length(args)<=3) {
  # default output prefix
  args[4] = "out"
}

# load the MSstats library 
require(MSstats)

# read dataframe into MSstats
data <- read.csv(args[1])
quant <- OpenMStoMSstatsFormat(data,
                               removeProtein_with1Feature = FALSE)

# process data
processed.quant <- dataProcess(quant, censoredInt = 'NA')

lvls <- levels(as.factor(data$Condition))
if (length(lvls) == 1)
{
  print("Only one condition found. No contrasts to be tested. If this is not the case, please check your experimental design.")
} else {
  if (args[2] == "pairwise")
  {
    if (args[3] == "")
    {
      l <- length(lvls)
      contrast_mat <- matrix(nrow = l * (l-1) / 2, ncol = l)
      rownames(contrast_mat) <- rep(NA, l * (l-1) / 2)
      colnames(contrast_mat) <- lvls
      c <- 1
      for (i in 1:(l-1))
      {
        for (j in (i+1):l)
        {
          comparison <- rep(0,l)
          comparison[i] <- -1
          comparison[j] <- 1
          contrast_mat[c,] <- comparison
          rownames(contrast_mat)[c] <- paste0(lvls[i],"-",lvls[j])
          c <- c+1
        }
      }
    } else {
      control <- which(as.character(lvls) == args[3])
      if (length(control) == 0)
      {
        stop("Control condition not part of found levels.n", call.=FALSE)
      }
      
      l <- length(lvls)
      contrast_mat <- matrix(nrow = l-1, ncol = l)
      rownames(contrast_mat) <- rep(NA, l-1)
      colnames(contrast_mat) <- lvls
      c <- 1
      for (j in setdiff(1:l,control))
      {
        comparison <- rep(0,l)
        comparison[i] <- -1
        comparison[j] <- 1
        contrast_mat[c,] <- comparison
        rownames(contrast_mat)[c] <- paste0(lvls[i],"-",lvls[j])
        c <- c+1
      }
    }
  }
  
  print ("Contrasts to be tested:")
  print (contrast_mat)
  #TODO allow for user specified contrasts
  test.MSstats <- groupComparison(contrast.matrix=contrast_mat, data=processed.quant)
  
  #TODO allow manual input (e.g. proteins of interest)
  write.csv(test.MSstats$ComparisonResult, "msstats_results.csv")
  
  groupComparisonPlots(data=test.MSstats$ComparisonResult, type="ComparisonPlot",
                       width=12, height=12,dot.size = 2,ylimUp = 7)
  
  test.MSstats$Volcano = test.MSstats$ComparisonResult[!is.na(test.MSstats$ComparisonResult$pvalue),]
  groupComparisonPlots(data=test.MSstats$Volcano, type="VolcanoPlot",
                       width=12, height=12,dot.size = 2,ylimUp = 7)
  
  if (nrow(constrast_mat) > 1)
  {
    groupComparisonPlots(data=test.MSstats$ComparisonResult, type="Heatmap",
                         width=12, height=12,dot.size = 2,ylimUp = 7)
  }
  
  # annotate how often the protein was quantified in each condition (NA values introduced by merge of completely missing are set to 1.0)
  ############ also calculate missingness on condition level

	# input: ProcessedData matrix of MSstats
	# output: 
	#   calculate fraction of na in condition (per protein)
	# Groups:   PROTEIN [762]
	#   PROTEIN                 `1`   `2`
	#   <fct>                 <dbl> <dbl>
	# 1 sp|A1ANS1|HTPG_PELPD   0    0.5  
	# 2 sp|A2I7N3|SPA37_BOVIN  0    0.5  
	# 3 sp|A2VDF0|FUCM_HUMAN   0    0.5  
	# 4 sp|A6ND91|ASPD_HUMAN   0.5  0.5  
	# 5 sp|A7E3W2|LG3BP_BOVIN  0.5  0.5  
	# 6 sp|B8FGT4|ATPB_DESAA   0    0.5

	getMissingInCondition <- function(processedData)
	{
		p <- processedData

		# count number of samples per condition
		n_samples = p %>% group_by(GROUP) %>% summarize(n_samples = length(unique((as.numeric(SUBJECT))))) 

		p <- p %>% 
		 filter(!is.na(INTENSITY)) %>% # remove rows with INTENSITY=NA
		 select(PROTEIN, GROUP, SUBJECT) %>%
		 distinct() %>% 
		 group_by(PROTEIN, GROUP) %>% 
		 summarize(non_na = n())  # count non-NA values for this protein and condition

		p <- left_join(p, n_samples) %>% 
				 mutate(missingInCondition = 1 - non_na/n_samples) # calculate fraction of missing values in condition

		# create one column for every condition containing the missingness
		p <- spread(data = p[,c("PROTEIN", "GROUP", "missingInCondition")], key = GROUP, value = missingInCondition)
		return(p)
	}

		mic <- getMissingInCondition(processed.quant$ProcessedData)

		test.MSstats$ComparisonResult <- merge(x=test.MSstats$ComparisonResult, y=mic, by.x="Protein", by.y="PROTEIN")
		commoncols <- intersect(colnames(mic), colnames(test.MSstats$ComparisonResult))
		test.MSstats$ComparisonResult[, commoncols]<-test.MSstats$ComparisonResult %>% select(commoncols) %>% mutate_all(funs(replace(., is.na(.), 1))) 

		#write comparison to CSV (one CSV per contrast)							 
		writeComparisonToCSV <- function(DF) 
		{
			write.table(DF, file=paste0("comparison_",unique(DF$Label),".csv"), quote=FALSE, sep='\t', row.names = FALSE)
			return(DF)
		}
		test.MSstats$ComparisonResult %>% group_by(Label) %>% do(writeComparisonToCSV(as.data.frame(.)))  

		#for (comp in rownames(contrast_mat))
		#{
		#  groupComparisonPlots(data=test.MSstats$ComparisonResult, type="ComparisonPlot",
		#                       width=12, height=12,dot.size = 2,ylimUp = 7, sig=1)#,
		#                       which.Comparison = comp,
		#                       address=F)
		#  # try to plot all comparisons
		#}
	}
