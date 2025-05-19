####################################################################################################################################################################################################################################


#Running Occupancy Models on Bat Activity# 


####################################################################################################################################################################################################################################

#data is a dataframe with all possible recording sites, and the presence of each species there#
data<- read.csv('merged_df.csv')

library(unmarked)
library(rlist)
library(MuMIn)
library(knitr)
library(flextable)
library(ggplot2)
library(spOccupancy)
library(tidyverse)
library(dplyr)
library(tidyr)
library(gridExtra)

#Filtering unsurveyed times and site 2019
data <- data %>%
  filter((time >= 170000) | (time <= 074500))

data <- data %>%
  filter(site != "2019")

data <- data %>%
  filter(!Cluster_32PCs %in% c("2", "4"))


#Split the data to do single species occupancy so only cluster 0 is present

#Convert date time to format
data$time_padded <- sprintf("%06d", as.numeric(data$time))
data$datetime <- as.POSIXct(paste(data$date, data$time_padded), format="%Y-%m-%d %H%M%S")

#Unique site_day combinations#
data$site_day <- paste(data$site, as.Date(data$datetime))

head(data)

#Changed 32PCs for Guild

effort <- data %>%
  mutate(value = ifelse(guild %in% c('Edge', 'Clutter', 'Open'), 1, 1)) %>%
  select(datetime, site, value) %>%
  pivot_wider(names_from = site, values_from = value, values_fill = list(value = NA))


#order alphabetically
site_cols <- colnames(effort)[-1] 
effort <- effort %>%
  select(datetime, sort(site_cols))


print(effort$datetime)


#################################################################################################################

        #Now to create one day Matrices

#################################################################################################################


#####PREPARE MATRICES####

data <- rename(data, Site = site) #to match function below
data <- rename(data, Species = guild) #change to match function
data <- rename(data, DateTime = datetime) #change to match function
data$DateTime <- as.POSIXct(data$DateTime)
unique(sort(data$Species)) #check



startDate <- as.POSIXct("2019-09-27 07:00:00", format="%Y-%m-%d %H:%M:%S") 
endDate <- as.POSIXct("2019-11-11 10:45:00", format="%Y-%m-%d %H:%M:%S") 
nrow(effort) #2597


all_cams<-effort #so matches function
str(data)
calcOcc <-
  function(species, # species name - in dataframe - that the function is to be run for
           d = d, # dataframe with species, site, and each date it was seen at that site - must have a columns called Species, Site and DateTime
           all_cams = all_cams, # matrix with all the survey dates, 1s for dates when a camera was working/deployed and NAs for when it wasn't
           startDate = startDate,#start date in date format
           endDate = endDate) {
    # Make a vector of breaks
    brks <-seq(startDate, endDate, by = "15 min")   #makes a sequence of all the days from start to end
    
    # Create an empty matrix of dim sites x time periods
    occ <-matrix(0, ncol = length(unique(d$Site)), nrow = length(brks))
    colnames(occ) <- sort(unique(d$Site))
    rownames(occ) <- strftime(brks, format = "%Y-%m-%d %H%M%S")
    
    for (s in unique(d$Site)) {
      #this loops through each site and inserts 1's on days which there were captures
      seen <- NA
      captures <-na.omit(d$DateTime[d$Species == species & d$Site == s])
      # Were animals seen at the site
      seen <- which(brks %in% captures)
      # If the species was seen, occ = 1
      col_i <- which(colnames(occ) == s)
      occ[seen, col_i] <- 1
    
    # Trim occ to match the number of rows in all_cams
    #if (nrow(occ) > nrow(all_cams)) {
      #print(paste("Trimming occ from", nrow(occ), "to", nrow(all_cams), "rows"))
      #occ <- occ[1:nrow(all_cams), ]
    }
    
    occ <- occ * all_cams[, 2:ncol(all_cams)]
    print(paste0(species, " done!"))
    species_name <- gsub(" ", "", species)
    
    #row.names(occ) <- brks
    
    write.csv(occ, paste0("1d_matrix_", species_name, ".csv"))
    return(occ)
    }

# applying function to each species (label)
lapply(
  X = unique(data$Species),
  FUN = calcOcc,
  d = data,
  all_cams=all_cams,
  startDate = startDate,
  endDate = endDate) # this will save CSVs of spp matrices in the working directory

###check there are 0s and 1s in the matrix###
filenames <- list.files("Matrices", pattern="*.csv", full.names=TRUE)
ldf <- lapply(filenames, read.csv) # reading all CSVs into a list
head(ldf)
str(ldf)
ldf[[1]] # just checking spp91 #shoat ##didnt do this
label <- basename(filenames)
label <- sub(pattern=".csv", replace="",label)
label <- sub(pattern="1d_matrix_", replace="",label)
names(ldf) <- label
new_list <- lapply(ldf, function(x) x%>% select(-X)) 
head(new_list)
names(new_list)
new_list[[1]]  #check 1st matrix
#

##############################
#THIS WORKS NOW AND HAS 1, this doesn't seem right

#Matrix 0 = 64 detections
#Matrix 1 = 377 detections
#Matrix 2 = 275 detections
#Matrix 3 = 25 detections
#Matrix 4 = 10 detections
#Matrix 5 = 69 detections
#Matrix 6 = 131 detections
species_counts <- data %>% 
  filter(Species %in% c("Clutter", "Edge", "Open")) %>%
  group_by(Species) %>%
  summarize(count = n())
print(species_counts)

#Raw species counts:
#0 = 892
#1 = 303
#2 = 99
#3 = 168
#4 = 33
#5 = 1591
#6 = 274 

#Clutter = 307
#Edge = 1154
#Open = 1496


#### aggregating days into survey occasions ####
##timestepper - creates matrices of a given timestep, can choose to include or exclude NAs

timestepper <- function(occ_in, timestep, na_mode = "include") {
  if (na_mode == "include") {
    occ_in[is.na(occ_in)] <-0 #replacing NAs with 0s if we want to include them in analysis.
  }
  if (timestep > nrow(occ_in) / 2) {
    print(paste(
      "Time step is too large! Please reduce to",
      nrow(occ_in) / 2 ,
      "or less."
    ))
  } else {
    start <- seq(1, nrow(occ_in), by = timestep)
    end <- seq(timestep, nrow(occ_in), by = timestep)
    if (length(start) > length(end)) {
      start <- start[-length(start)]
    }
    timesteps <- matrix(nrow = length(start), ncol = ncol(occ_in))
    colnames(timesteps) <- colnames(occ_in)
    rownames(timesteps) <-
      paste(rownames(occ_in)[start], rownames(occ_in)[end], sep = ":")
    for (i in 1:length(start)) {
      timestep_out <- colSums(occ_in[start[i]:end[i], ])
      timesteps[i, ] <- timestep_out
      timesteps[timesteps > 0] <- 1
    }
    timesteps <- t(timesteps)
    return(timesteps)
  }
}


matrix_1d <- lapply(X = new_list,
                    FUN = timestepper,
                    timestep = 4, #Grouping by 15 minutes so here one hour
                    na_mode = "exclude") # must be exclude otherwise will consider days when CT wasn't working as actual survey days
matrix_1d[[1]]
names(matrix_1d) <- label # adding species names


##############################





#Select covariates
covariates <- read.csv("AM_covs.csv")

##############################
#Variables are:
#Site - '1'
#Porportion of open veg (500m) - '5'
#Distance to human settlements - '19'
#Distance to water - '14'
#Grazing - as from camera traps - '34/35'
#SAVI (mean) - '36'


covs_numeric<-covariates[ , c(5, 14, 19, 34, 35, 36)]

#Removing NAs
covs_numeric$cattle_30min_event_rate[is.na(covs_numeric$cattle_30min_event_rate)] <- 0
covs_numeric$shoat_30min_event_rate[is.na(covs_numeric$shoat_30min_event_rate)] <- 0


covs_numeric <- as.data.frame(scale(covs_numeric))

#Extracting site seperately#
col_1 <- covariates[, 1, drop = FALSE]  
covs_numeric <- cbind(col_1, covs_numeric)


#Is it possible to add lunar illumination when it is based on the date of the survey?#


head(covs_numeric)





#Changing to guild
cluster0 <- list.extract(matrix_1d, "Clutter")
head(cluster0)
sum(cluster0, na.rm = TRUE)

cluster1 <- list.extract(matrix_1d, "Edge")
head(cluster1)
sum(cluster1, na.rm = TRUE)

cluster2 <- list.extract(matrix_1d, "Open")
head(cluster2)
sum(cluster2, na.rm = TRUE)


all_bats <- cluster0 + cluster1 + cluster2
all_bats[all_bats > 0] <- 1
sum(all_bats, na.rm = TRUE)



        #Extract matrix for each cluster
        cluster0 <- list.extract(matrix_1d, 1)
        head(cluster0)
        sum(cluster0, na.rm = TRUE)
        
        cluster1 <- list.extract(matrix_1d, 2)
        head(cluster1)
        sum(cluster1, na.rm = TRUE)
        
        cluster2 <- list.extract(matrix_1d, 3)
        head(cluster2)
        sum(cluster2, na.rm = TRUE)
        
        cluster3 <- list.extract(matrix_1d, 4)
        head(cluster3)
        sum(cluster3, na.rm = TRUE)
        
        cluster4 <- list.extract(matrix_1d, 5)
        head(cluster4)
        sum(cluster4, na.rm = TRUE)
        
        cluster5 <- list.extract(matrix_1d, 6)
        head(cluster5)
        sum(cluster5, na.rm = TRUE)
        
        cluster6 <- list.extract(matrix_1d, 7)
        head(cluster6)
        sum(cluster6, na.rm = TRUE)
        
        
        
        all_bats <- cluster0 + cluster1 + cluster2 + cluster3 + cluster4 + cluster5 + cluster6
        all_bats[all_bats > 0] <- 1
        sum(all_bats, na.rm = TRUE)
        

print(table(all_bats))  # Should show counts of 0s and 1s


#####################################
#Need to find the covariates now online


covs_sites <- unique(covariates$site)
my_sites <- colnames(effort)

non_CT_sites <- setdiff(my_sites, covs_sites)
#now there are no sites missing 
##############################





##############################
#Filtering by sites in my subset
cluster1_sites <- setdiff(my_sites, non_CT_sites)


covs_numeric <- covs_numeric[covs_numeric$site %in% cluster1_sites, ]
cluster0 <- cluster0[rownames(cluster0) %in% cluster1_sites, , drop = FALSE]
cluster1 <- cluster1[rownames(cluster1) %in% cluster1_sites, , drop = FALSE]
cluster2 <- cluster2[rownames(cluster2) %in% cluster1_sites, , drop = FALSE]
cluster3 <- cluster3[rownames(cluster3) %in% cluster1_sites, , drop = FALSE]
cluster4 <- cluster4[rownames(cluster4) %in% cluster1_sites, , drop = FALSE]
cluster5 <- cluster5[rownames(cluster5) %in% cluster1_sites, , drop = FALSE]
cluster6 <- cluster6[rownames(cluster6) %in% cluster1_sites, , drop = FALSE]

all_bats <- all_bats[rownames(all_bats) %in% cluster1_sites, , drop = FALSE]
dim(covs_numeric)
dim(all_bats)



#Creating dataframe for unmarked fuction for each cluster
umf.cluster0<- unmarkedFrameOccu(y = cluster0, siteCovs = covs_numeric)
summary(umf.cluster0)

  umf.cluster1<- unmarkedFrameOccu(y = cluster1, siteCovs = covs_numeric)
  summary(umf.cluster1)

  umf.cluster2<- unmarkedFrameOccu(y = cluster2, siteCovs = covs_numeric)
  summary(umf.cluster2)

  umf.cluster3<- unmarkedFrameOccu(y = cluster3, siteCovs = covs_numeric)
  summary(umf.cluster3)
  
  umf.cluster4<- unmarkedFrameOccu(y = cluster4, siteCovs = covs_numeric)
  summary(umf.cluster4)
  
  umf.cluster5<- unmarkedFrameOccu(y = cluster5, siteCovs = covs_numeric)
  summary(umf.cluster5)
  
  umf.cluster6<- unmarkedFrameOccu(y = cluster6, siteCovs = covs_numeric)
  summary(umf.cluster6)
  
  
umf.all_bats<- unmarkedFrameOccu(y = all_bats, siteCovs = covs_numeric)
summary(umf.all_bats)
  
#Debugging HERE!!!!!! - HOW CAN I FIX THE MASSIVE SE Values
##############################

#Null model:
#P is the detection
#psi is the occupancy 

m0BF <- occu(~1 ~ 1, data=umf.cluster0) #null #After first tilda in the model is detection after second is occupancy predictors
                    #Instead of 1 would change to propopen500m + waterdist etc
                    #When maybe I have hourly precipitation can change the first 1 to precipitation or wind as it may alter the pr(detection)
m0BF

any(is.na(umf.cluster0))


      m1BF <- occu(~1 ~ 1, data=umf.cluster1)
      m1BF
      
      m2BF <- occu(~1 ~ 1, data=umf.cluster2)
      m2BF
      
      m3BF <- occu(~1 ~ 1, data=umf.cluster3)
      m3BF
      
      m4BF <- occu(~1 ~ 1, data=umf.cluster4)
      m4BF
      
      m5BF <- occu(~1 ~ 1, data=umf.cluster5)
      m5BF
      
      m6BF <- occu(~1 ~ 1, data=umf.cluster6)
      m6BF



mall_batsBF <- occu( ~1 ~ 1, data=umf.all_bats)
mall_batsBF


#detection doesn't vary with any covariate
#summary(m1WB)

    m0BF.p<-backTransform(mall_batsBF, type="det") #detection #chance of observation
    m0BF.psi<-backTransform(mall_batsBF, type="state") #occupancy #prop. of observation
    
    
    m0BFCI.psi<-confint(m0BF.psi)
    m0BFCI.p<-confint(m0BF.p)



##############################

#Holly's Plotting Code:
#https://github.com/hollypringle/occupancy/blob/master/allspeciesploteasy.R 


##############################


#global model used for dredging 
mfull.all_bats <- occu( ~ 1 ~ propopen500m + cattle_30min_event_rate + shoat_30min_event_rate + waterdist_short + humdist_short + Mean.savi, data=umf.cluster2)
summary(mfull.all_bats)


# dredge function - use ir wisely! Think hard about variables included in the global model and don't let the computer do the thinking...
drg.all_bats<- dredge(mfull.all_bats, rank = AIC) #, fixed = ('psi(pa_type)')) #dredging models
drg.all_bats


################

#Preference is for the first model

#model-average estimates considering models with delta AIC <2
avg1.cluster_0<-model.avg(drg.all_bats, subset=delta<=2, fit=T) #model-average estimates
summary(avg1.cluster_0)


      #Printing nice table
            summ <- summary(avg1.cluster_0)
            
            # Full model-averaged coefficients (full average)
            coef_table_full <- as.data.frame(summ$coefmat.full)
            
            # Conditional model-averaged coefficients (conditional average)
            coef_table_conditional <- as.data.frame(summ$coefmat.subset)
            
            # View them nicely
            
            kable(coef_table_full, digits = 3)
            kable(coef_table_conditional, digits = 3)
           
            coef_table_conditional$Variable <- rownames(coef_table_conditional)
            coef_table_conditional <- coef_table_conditional[, c("Variable", "Estimate", "Std. Error", "z value", "Pr(>|z|)")]
            
            coef_table_conditional$Estimate <- round(coef_table_conditional$Estimate, 3)
            coef_table_conditional$`Std. Error` <- round(coef_table_conditional$`Std. Error`, 3)
            coef_table_conditional$`z value` <- round(coef_table_conditional$`z value`, 3)
            coef_table_conditional$`Pr(>|z|)` <- round(coef_table_conditional$`Pr(>|z|)`, 3)
            
            # Add significance asterisks based on p-value thresholds
            coef_table_conditional$Significant <- ""
            coef_table_conditional$Significant[coef_table_conditional$`Pr(>|z|)` < 0.05] <- "*"
            coef_table_conditional$Significant[coef_table_conditional$`Pr(>|z|)` < 0.01] <- "**"
            coef_table_conditional$Significant[coef_table_conditional$`Pr(>|z|)` < 0.001] <- "***"
            
            flextable(coef_table_conditional)
            
            
            m0BF.p<-backTransform(mall_batsBF, type="det") #detection #chance of observation
            m0BF.psi<-backTransform(mall_batsBF, type="state")
            
            m0BFCI.psi<-confint(m0BF.psi)
            m0BFCI.p<-confint(m0BF.p)

#model-average estimates considering models up to cumulative AIC weight of .95

#avg1.cluster_0 <-model.avg(drg.all_bats, subset = cumsum(weight) <= .95, fit=T) 
#summary(avg1.cluster_0)
      #p(INT) - highly significant the rest are not? 



#predicting occupancy based on continuous variable, and a fixed management regime

#HUMAN
range(covs_numeric$humdist_short)
newDatahuman=data.frame(humdist_short=seq(-1.316735, 3.071820,by=0.1), propopen500m=0, cattle_30min_event_rate=0, shoat_30min_event_rate=0, waterdist_short=0, Mean.savi=0)
pred.avgfull.human.cluster_0<-predict(avg1.cluster_0, type="state", se.fit=TRUE, full=T, newDatahuman)
prdOccucluster_0human <- as.data.frame(pred.avgfull.human.cluster_0)
prdOccucluster_0human$humdist_short <- newDatahuman$humdist_short
prdOccucluster_0human

#WATER - NUMBERS ARE NOVEL
range(covs_numeric$waterdist_short)
newDatawater=data.frame(waterdist_short=seq(-1.294211,  3.552385,by=0.1), propopen500m=0, cattle_30min_event_rate=0, shoat_30min_event_rate=0, humdist_short=0, Mean.savi=0)
pred.avgfull.water.cluster_0<-predict(avg1.cluster_0, type="state", se.fit=TRUE, full=T, newDatawater)
prdOccucluster_0water <- as.data.frame(pred.avgfull.water.cluster_0)
prdOccucluster_0water$waterdist_short <- newDatawater$waterdist_short
prdOccucluster_0water

#OPEN - NUMBERS ARE NOVEL
  #large number of warning messages arise on this block...
range(covs_numeric$propopen500m)
newDataopen=data.frame(propopen500m=seq(-2.05379, 1.16234,by=0.1), waterdist_short=0, cattle_30min_event_rate=0, shoat_30min_event_rate=0, humdist_short=0, Mean.savi=0)
pred.avgfull.open.cluster_0<-predict(avg1.cluster_0, type="state", se.fit=TRUE, full=T, newDataopen)
prdOccucluster_0open <- as.data.frame(pred.avgfull.open.cluster_0)
prdOccucluster_0open$propopen500m <- newDataopen$propopen500m
prdOccucluster_0open

#SAVI - NUMBERS ARE NOVEL
  #large number of warning messages arise on this block...
range(covs_numeric$Mean.savi)
newDatasavi=data.frame(Mean.savi=seq(-2.249580,  2.452621,by=0.1), waterdist_short=0, cattle_30min_event_rate=0, shoat_30min_event_rate=0, humdist_short=0, propopen500m=0)
pred.avgfull.savi.cluster_0<-predict(avg1.cluster_0, type="state", se.fit=TRUE, full=T, newDatasavi)
prdOccucluster_0savi <- as.data.frame(pred.avgfull.savi.cluster_0)
prdOccucluster_0savi$Mean.savi <- newDatasavi$Mean.savi
prdOccucluster_0savi

#SHOAT
range(covs_numeric$shoat_30min_event_rate)
newDatashoat=data.frame(shoat_30min_event_rate=seq(-0.3659014,  5.4965962,by=0.1), propopen500m=0, cattle_30min_event_rate=0, humdist_short=0, waterdist_short=0, Mean.savi=0)
pred.avgfull.shoat.cluster_0<-predict(avg1.cluster_0, type="state", se.fit=TRUE, full=T, newDatashoat)
prdOccucluster_0shoat <- as.data.frame(pred.avgfull.shoat.cluster_0)
prdOccucluster_0shoat$shoat_30min_event_rate <- newDatashoat$shoat_30min_event_rate
prdOccucluster_0shoat

#CATTLE
range(covs_numeric$cattle_30min_event_rate)
newDatacattle=data.frame(cattle_30min_event_rate=seq(-0.4140332,  6.6212613,by=0.1), propopen500m=0, humdist_short=0, shoat_30min_event_rate=0, waterdist_short=0, Mean.savi=0)
pred.avgfull.cattle.cluster_0<-predict(avg1.cluster_0, type="state", se.fit=TRUE, full=T, newDatacattle)
prdOccucluster_0cattle <- as.data.frame(pred.avgfull.cattle.cluster_0)
prdOccucluster_0cattle$cattle_30min_event_rate <- newDatacattle$cattle_30min_event_rate
prdOccucluster_0cattle


#PLOT
cluster_0humanplot<-ggplot(prdOccucluster_0human, aes(x=humdist_short, y=fit)) +
  theme_bw()+ # white background (as opposed to the default grey)
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_ribbon(aes(ymin= fit-(se.fit*1.96), ymax=fit+(se.fit*1.96)), alpha=0.2)+
  geom_line(linewidth=1,colour="blue")+
  coord_cartesian(ylim = c(0, 1), xlim = c(-1.224601, 2.832376))+
  theme(
    panel.background = element_blank(),      
    panel.grid.major = element_blank(),      
    panel.grid.minor = element_blank(),     
    panel.border = element_blank(),          
    axis.line = element_line(colour = "black"), 
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    plot.title = element_text(size = 18)
  ) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(legend.position="none")+
  ylab("Probability of Occupancy")+
  xlab("Distance to human structure")+
  ggtitle("")

cluster_0cattleplot<-ggplot(prdOccucluster_0cattle, aes(x=cattle_30min_event_rate, y=fit)) +
  theme_bw()+ # white background (as opposed to the default grey)
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_ribbon(aes(ymin= fit-(se.fit*1.96), ymax=fit+(se.fit*1.96)), alpha=0.2)+
  coord_cartesian(ylim = c(0, 1), xlim = c(-0.5149533, 4.6834225))+
  geom_line(linewidth=1,colour="blue")+
  theme(
    panel.background = element_blank(),      
    panel.grid.major = element_blank(),      
    panel.grid.minor = element_blank(),     
    panel.border = element_blank(),          
    axis.line = element_line(colour = "black"), 
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    plot.title = element_text(size = 18)
  ) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  #scale_fill_manual(values=c("yellow4", "forestgreen", "tomato3"))+
  theme(legend.position="none")+
  ylab("Probability of Occupancy")+
  xlab("Cattle grazing")+
  ggtitle("")

cluster_0shoatplot<-ggplot(prdOccucluster_0shoat, aes(x=shoat_30min_event_rate, y=fit)) +
  theme_bw()+ # white background (as opposed to the default grey)
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_ribbon(aes(ymin= fit-(se.fit*1.96), ymax=fit+(se.fit*1.96)), alpha=0.2)+
  coord_cartesian(ylim = c(0, 1), xlim = c(-0.4168032, 4.4785513))+
  geom_line(linewidth=1,colour="blue")+
  theme(
    panel.background = element_blank(),      
    panel.grid.major = element_blank(),      
    panel.grid.minor = element_blank(),     
    panel.border = element_blank(),          
    axis.line = element_line(colour = "black"), 
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    plot.title = element_text(size = 18)
  ) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  #scale_fill_manual(values=c("yellow4", "forestgreen", "tomato3"))+
  theme(legend.position="none")+
  ylab("")+
  xlab("Shoat grazing")+
  ggtitle("")

cluster_0waterplot<-ggplot(prdOccucluster_0water, aes(x=waterdist_short, y=fit)) +
  theme_bw()+ # white background (as opposed to the default grey)
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_ribbon(aes(ymin= fit-(se.fit*1.96), ymax=fit+(se.fit*1.96)), alpha=0.2)+
  coord_cartesian(ylim = c(0, 1), xlim = c(-1.327513, 3.395769))+
  geom_line(linewidth=1,colour="blue")+
  theme(
    panel.background = element_blank(),      
    panel.grid.major = element_blank(),      
    panel.grid.minor = element_blank(),     
    panel.border = element_blank(),          
    axis.line = element_line(colour = "black"), 
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    plot.title = element_text(size = 18)
  ) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  #scale_fill_manual(values=c("yellow4", "forestgreen", "tomato3"))+
  theme(legend.position="none")+
  ylab("")+
  xlab("Distance to water")+
  ggtitle("All bats (12 hours)")

cluster_0saviplot<-ggplot(prdOccucluster_0savi, aes(x=Mean.savi, y=fit)) +
  theme_bw()+ # white background (as opposed to the default grey)
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_ribbon(aes(ymin= fit-(se.fit*1.96), ymax=fit+(se.fit*1.96)), alpha=0.2)+
  coord_cartesian(ylim = c(0, 1), xlim = c(-2.143278, 2.102924))+
  geom_line(linewidth=1,colour="blue")+
  theme(
    panel.background = element_blank(),      
    panel.grid.major = element_blank(),      
    panel.grid.minor = element_blank(),     
    panel.border = element_blank(),          
    axis.line = element_line(colour = "black"), 
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    plot.title = element_text(size = 18)
  ) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  #scale_fill_manual(values=c("yellow4", "forestgreen", "tomato3"))+
  theme(legend.position="none")+
  ylab("")+
  xlab("SAVI")+
  ggtitle("")

cluster_0openplot<-ggplot(prdOccucluster_0open, aes(x=propopen500m, y=fit)) +
  theme_bw()+ # white background (as opposed to the default grey)
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_ribbon(aes(ymin= fit-(se.fit*1.96), ymax=fit+(se.fit*1.96)), alpha=0.2)+
  coord_cartesian(ylim = c(0, 1), xlim = c(-2.332393, 1.034469))+
  geom_line(linewidth=1,colour="blue")+
  theme(
    panel.background = element_blank(),      
    panel.grid.major = element_blank(),      
    panel.grid.minor = element_blank(),     
    panel.border = element_blank(),          
    axis.line = element_line(colour = "black"), 
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    plot.title = element_text(size = 18)
  ) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  #scale_fill_manual(values=c("yellow4", "forestgreen", "tomato3"))+
  theme(legend.position="none")+
  ylab("")+
  xlab("Proportion of open habitat")+
  ggtitle("")


cluster_0plot<-grid.arrange(cluster_0humanplot, cluster_0waterplot, cluster_0openplot, cluster_0cattleplot, cluster_0shoatplot, cluster_0saviplot, ncol = 3)

######################################################################################################################################################