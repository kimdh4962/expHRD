#############################
# Overview
# R packages for ssGSEA
# ssGSEA is one of GSVA options
#############################

#############################
# Required packages
# GSVA: For ssGSEA running and gene set determination
# ggplot2: For graph visualization
#############################

#############################
# Required data
# Gene set: including 356 genes which represent HRD
# Input Data: DESEQ2 normalized RNAseq data, CSV format
#############################

#############################
# Parmeters 
# Data: RNAseq Data frame
# gene_set: Gene set for GSVA
# ssGSEA.es (may be changed): result data frame of ssGSEA
#############################


library('GSVA') # sever must be installed GSVA
library(ggplot2) # graph plotting
library('ggpubr')

########running time
args = c('0','./../data/TCGA_TNBC_DEG_test.ENSG.drop.csv', './../output/OUT_TABLE.csv', './../output/OUT_plot_DEG.png', 'confidence', 0.95)

# Firstly, I have to set only 2 args 
# args[1] = option 1: DESEQ2 int
# args[2] = input data ({INPUT_GENE_EXPRESSION.csv}/ split by '.csv' ) input query RNA expression, default data is TCGA_TNBC_DEG
# args[3] = output table ({OUTPUT_TABLE.csv} )
# args[4] = plot save ({OUTPUT_FIGURE.png} )
# args[5] = "confidence" or "prediction" interval
# args[6] = percentage = 0.95 or etc.
# all args follow Model of DATABASE
#print(args)


######################################
# input data read : input -> RNA seq
read_RNA <- function(option){
  #print(option)
  if (option == "0"){ # deseq2
    data <- read.csv('./../data/TCGA_OV_DEG_test.ENSG.drop.csv', header = TRUE,sep =',',  row.names = 1, check.names = FALSE)
  }
  #print(head(data))
  return(data)
}

#!#!#!#!!# real data and plot
# background RNA data
bacK_data <- read_RNA(args[1]) # following option args
##!#!#!##!#!#!#!##

# to check lm equation
#bacK_data <- read_RNA(1)
#data <- read_RNA(1)
#data <- back_data #as.matrix(back_data)

# data from input query -> use after drawing and regression of background RNA
#!#!#!##!#
data <- read.csv(args[2], header = TRUE,sep =',',  row.names = 1, check.names = FALSE)
#!#!#!#!#!

###########################
# data : RNA expression normalized by TPM or DESEQ2
##print(colnames(data)[1:10])
data <- as.matrix(data)
back_data <- as.matrix(bacK_data)

#########################
# import gene sets : 356 genes 
gene_set_before <- readLines('./../data/FINAL_98_GENE_ENSG_ID_repo.csv')
#print(gene_set_before)

##################
# make gene set
gene_set  <- list()
bigname  <- c()
for (m in 1:length(gene_set_before)){
  names  <- c()
  temp <- c()
  j=1
  for (i in unlist(strsplit(gene_set_before[m],','))){
    if (i ==''){
      break
    }
    if (j==1){
      names <- append(names, i)
      j <- j+1
      next
    }
    temp <- append(temp,i)
  }
  bigname <- append(bigname, names)
  gene_set  <- append(gene_set,list(names= temp ))
}
names(gene_set)  <- bigname
#print(gene_set)
#############################


###################
# GSVA -> ssGSEA option -> calculate ssGSEA score
ssGSEA <- gsva(back_data, gene_set, verbose=F, method='ssgsea', ssgsea.norm=F)
ssGSEA[3,] <- ssGSEA[1,] - ssGSEA[2,] # adjustment of ssGSEA score
# ssGSEA positive - ssGSEA negative

##################
# import HRD score from TCGA OV test sets (n=58)
hrd_score <- read.csv('./../data/OV_HRD.csv', row.names = 1, check.names = F) # pre-calculated HRD score
# do not changed

new_df <- cbind(hrd_score$HRD, ssGSEA[3,])
colnames(new_df) <- c('HRD', 'expHRD')
new_df=as.data.frame(new_df)
#print(head(new_df))

####################### 
# regression : HRD score ~ expHRD
model <- lm(HRD ~ expHRD, data=new_df)
#confint(model)
##################
# HRD prediction function from linear regression model
# cal_HRD <- function(expHRD, model){
#   low_HRD <- confint(model)[1] + expHRD * confint(model)[4] # low b0 + high b1
#   HRD <- model$coefficients[1] + expHRD * model$coefficients[2] # b0 + b1
#   high_HRD <- confint(model)[3] + expHRD * confint(model)[2] # high b0 + low b1
#   
#   df <- data.frame(expHRD,low_HRD,HRD,high_HRD )
#   colnames(df) <- c('expHRD', 'LOWER_HRD', 'p_HRD','HIGHER_HRD')
#   
#   for (x in seq(1,length(rownames(df)))) {
#     if (df[x,4] < df[x,2]) {
#       temp <- df[x,2]
#       df[x,2] <- df[x,4]
#       df[x,4] <- temp
#     }
#   }
#   
#   return(df)
# }
##################

cal_HRD <- function(expHRD, model){
  temp_df <- data.frame(expHRD)
  colnames(temp_df) <- c("expHRD")
  results <- predict(model, newdata = temp_df, interval = args[5], level=as.numeric(args[6]) ) # user can choose confidence interval and rate
  #return(results)
  df <- data.frame(expHRD, results[,'lwr'], results[,'fit'], results[,'upr'])
  colnames(df) <- c('expHRD', 'LOWER_HRD', 'p_HRD','HIGHER_HRD')
  return(df)
  
}

# return new data frame with predicted (low, mid, high HRD by regression model)
########################
# validation of our regression models
# probability whether expHRD correlate with true HRD (by scarHRD)
validation_of_model <- cbind(cal_HRD(new_df['expHRD'], model),new_df['HRD'])

#####################
# calculate probability and the number of outlier

# condition.1
# -> it may overly estimate than true HRD
# (ex. lower expHRD > true HRD) => outlier
# condition.2
# -> it may less estimate than true HRD
# (ex. higher expHRD < true HRD) => outlier
# condition.3
# -> the others : the range of expHRD located in true HRD (CI=95%)

# result
# Probability = 1 - (# of outlier/ # total)
# Originally, our TCGA-OV test set represents 0.6552 probability of regression

#######################

# probability calculation only for base table (TCGA-OV)
probability <- function(df){
  df <- df
  out_range <- 0
  bins <- c()
  for (x in seq(1,length(rownames(df))) ){
    #print(df[x,])
    if (df[x,5] > max(df[x,2:4]) ) {
      out_range <- out_range + 1
      bins <- append(bins, 'YES')
      #print('HRD is over CI95')
    }
    else if (df[x,5] < min(df[x,2:4])){
      #print('HRD is under CI05')
      out_range <- out_range + 1
      bins <- append(bins, 'YES')
    }
    else{
      bins <- append(bins, 'NO')
      #print('HRD is within in prediction')
    }
  }
  prob <- round(x=(1-(out_range/length( rownames(df) ) ) ),digits = 4 )
  #print(bins)
  df <- cbind(df,bins)
  df$prob <- prob
  
  colnames(df) <- c('expHRD', 'CI05_HRD', 'p_HRD','CI95_HRD',
                    'true_HRD',paste0('outlier(n=',out_range,')'),
                    paste0('prob=(',prob,')'))
  return(df)
}

# result of base table
base_table <- probability(validation_of_model)
# do not express anywhere. No body catch the table
# Base table is not important to our true result
# We have to show another query and its prediction table, but we have already made
# model (TCGA-OV with 58 samples) -> which determines our new query with 0.6552 probability
# make one-pipeline to calculate expHRD score with our query

#######################
# start calculate

ssGSEA.qu <- gsva(data, gene_set, verbose=F, method='ssgsea', ssgsea.norm=F)
ssGSEA.qu[3,] <- ssGSEA.qu[1,] - ssGSEA.qu[2,] # adjustment of ssGSEA score

# making result table and predictions 
result_table <- cal_HRD(ssGSEA.qu[3,],model)
temp_result <- data.frame(result_table$expHRD)
colnames(temp_result) <- c('expHRD')

# optional value for interval
interval <- predict(model, newdata=temp_result, interval=args[5],
                    level = as.numeric(args[6]))
rownames(interval) <- rownames(result_table)
interval <- cbind(interval, temp_result)


##########################
# time to plot scatter plot
# color = gray
# regression line = black
#p <- ggplot(new_df, aes(x=expHRD, y=HRD), fig(10,10))  +
#  stat_smooth(method = lm,formula = y~x, color='gray')+ theme_bw() + 
#  theme(panel.grid.major = element_blank(),
#        panel.grid.minor = element_blank(),
#        panel.background = element_rect(colour = "black", size=2),
#        axis.text.x = element_text(color="black", face='bold',size=12),
#        axis.text.y = element_text(color="black", face='bold',size=12),
#        axis.title.x = element_text(color='black', size=15,face='bold' ),
#        axis.title.y = element_text(color='black', size=15,face='bold'),
#        axis.ticks.length =  unit(0.2,'cm'),
#        plot.title =element_text(hjust=0.5, face='bold'),)+
#  stat_cor(label.x=0)+
#  geom_point(data=result_table,aes(x=expHRD,y=p_HRD), color='red') + # optional 
  
#  geom_line(data=interval, aes(x=expHRD, y=lwr), col='blue',linetype = "dashed")+
#  geom_line(data=interval, aes(x=expHRD, y=upr), col='blue',linetype = "dashed")+ylim(-20, 120) #expand_limits(x=0, y=c(0,100))


#p
### test ###
#ggsave(args[4], device = "png", units='px', width=1500, height=1500)
###      ###
#ggsave(paste0('./figures/','OUT_',substr(args[5],1,1),substr(args[6],3,5),'.png'),  device = "png", units='px', width=1500, height=1500 )
#pdf(args[4])
#print(p)
#dev.off()

#######################
# save data frame result -> after ssGSEA and after regression equation adjustment
write.csv(result_table, args[3])


