# Load packages ####
listepackages<-c("glue","ggplot2","readr","lubridate","tidyverse","knitr","DT","rmarkdown","cowplot","TSclust","fpc","ape","dendextend","cowplot")
for (pack in listepackages) {
  if (!is.element(pack, installed.packages()[,1])){
    install.packages(pack, dep = TRUE)
  }
  eval(parse(text=paste0("library(",pack,")")))
}
rm(pack)

# Import Data ####
data_cim <- readRDS(paste0(rept2,"datacim.rds"))
data_pat <- readRDS(paste0(rept2,"datapat.rds"))
data_nsej <- readRDS(paste0(rept2,"datasej.rds"))

# Create directories ####
rept2 <- paste0(getwd(),"/") 

if(!dir.exists(paste0(rept2,"Reports/"))){
  dir.create(paste0(rept2,"Reports/"),showWarnings = F)
}
if(!dir.exists(paste0(rept2,"Graphs/"))){
  dir.create(paste0(rept2,"Graphs/"),showWarnings = F)
}
if(!dir.exists(paste0(rept2,"DataProfil/"))){
  dir.create(paste0(rept2,"DataProfil/"),showWarnings = F)
}

# Data management ####

data_pat$date <- ymd(paste0(data_pat$annee,"-",data_pat$mois,"-01"))
data_nsej$date <- ymd(paste0(data_nsej$annee,"-",data_nsej$mois,"-01"))
# Merging data and Derived variables
data_all <- left_join(data_cim,data_pat,by=c("mois","annee"))
data_all <- data_all %>% mutate(freq=N/Npatients)
data_all <- data_all %>% mutate(jour="1")
data_all <- data_all %>% unite(temps,jour,mois,annee, sep="/")
data_all <- data_all %>% mutate(temps_num=as.numeric(dmy(temps)))
data_all <- data_all %>% mutate(mois=month(dmy(temps)),annee=year(dmy(temps)))

# Preprocessing (generate plots, compute moving mean) ####

# 3-character ICD codes frequencies over time ####

if(!file.exists(paste0(rept2,"Reports/cim3movingmean.csv"))){
  liste_codes <- unique(data_all$diag2)
  # Output for amplitude computation
  col_bp <- c("ICD","Times","RelativeAmplitude","MinFreqs","MaxFreqs")
  res_bp <- data.frame(matrix(NA,nrow=length(liste_codes), ncol=length(col_bp) ))
  res_bp[,1] <- liste_codes
  colnames(res_bp) <- col_bp
  # Output moving mean for clustering
  dtt_all_mov <- data.frame()
  for(k in 1:length(liste_codes)){
    
    print(paste0(k,": ",liste_codes[k]))
    dtt2 <- data_all %>% filter(diag2==liste_codes[k])
    
    # Amplitudes computation ####
    if(length(unique(dtt2$temps_num))>15){
      dtt2$freqs <- NA
      i <- 6
      for(j in (1+i):(length(unique(dtt2$temps_num))-i)){
        temp <- 0
        for(w in -i:i){
          temp <- temp + dtt2[j+w,"freq"]
        }
        dtt2[j,"freqs"] <- temp/(2*i+1)
      }
      dtt_all_mov <- rbind(dtt_all_mov,dtt2[,c("diag2", "temps", "freqs")])
      res_bp[k,"MinFreqs"] <- min(dtt2[,"freqs"], na.rm=T)
      res_bp[k,"MaxFreqs"] <- max(dtt2[,"freqs"], na.rm=T)
      res_bp[k,"RelativeAmplitude"] <- abs(res_bp[k,"MinFreqs"]-res_bp[k,"MaxFreqs"])/mean(res_bp[k,"MinFreqs"],res_bp[k,"MaxFreqs"])
    }else{
      res_bp[k,"MinFreqs"] <- NA
      res_bp[k,"MaxFreqs"] <- NA
    }
    
    # Plots ####
    obj <- ggplot(dtt2) + geom_point(aes(x=as.Date(temps_num,origin="1970-01-01"),y=freq)) + theme_classic() + ylab("Frequency") + xlab(paste0("ICD10 code:",liste_codes[k]))
    
    if(!is.null(dtt2$freqs)){
      obj <- obj + geom_point(aes(x=as.Date(temps_num,origin="1970-01-01"),y=freqs),shape=3)
    }
    
    ggsave(obj,filename = paste0(paste0(rept2,"Graphs/"),liste_codes[k],".png"),device = "png",dpi = 300)
  }
  
  # Save outputs ####
  write.csv2(res_bp, paste0(rept2,"Reports/cim3freq.csv"), row.names = F)
  write.csv2(dtt_all_mov,paste0(rept2,"Reports/cim3movingmean.csv"), row.names = F)
}else{
  print("Files already generated, delete the Reports directory to relaunch computing")
}

# Clustering SAX ####

# Data preparation ####
dtt_all_mov <- read.csv2(paste0(rept2,"Reports/cim3movingmean.csv"),stringsAsFactors = F)
ICDdata <- dtt_all_mov
ICDdata <- ICDdata[!is.na(ICDdata$freqs),]
ICDdata <- ICDdata %>%
  mutate(temps = temps %>% as.Date(format = "%d/%m/%Y")) 

dataKLM_wide <- ICDdata %>% select(diag2, temps, freqs) %>% spread(temps, freqs)
dataKLM_wide <- as.data.frame(dataKLM_wide)  
# dataKLM_wide[is.na(dataKLM_wide)] <- 0
colmax <- length(dataKLM_wide)
vectmed <- apply(dataKLM_wide[,-1],1,function(x) mean(x,na.rm=T))
vectsd <- apply(dataKLM_wide[,-1],1,function(x) sd(x,na.rm=T))
dataKLM_long <- (dataKLM_wide[,2:colmax]-vectmed)/vectsd
row.names(dataKLM_long) <- dataKLM_wide[,1]
dataKLM_long <- as.matrix(dataKLM_long)
dimnames(dataKLM_long)[[2]] <- NULL
dataKLM_long <- drop_na(as.data.frame(dataKLM_long))
dataKLM_long = as.data.frame(t(dataKLM_long))

# Descriptive statistics ####

allres <- read.csv2(paste0(rept2,"Reports/cim3freq.csv"))
allres <- allres[,c("ICD","RelativeAmplitude")]
allres <- allres %>% rename(RAfreq=RelativeAmplitude)
# 3.0 Number of encounters ####
sum(data_nsej$Nsejours)
# 3.1	How many ICD codes have temporal instability? ####
sum(allres$RAfreq<0.5, na.rm = T)
sum(allres$RAfreq>0.5 & allres$RAfreq<=1, na.rm = T)
sum(allres$RAfreq>1, na.rm = T)


# Clustering SAX based on the Symbolic Aggregate approXimation measure

# Compute the distances ####
distSAX_w20a6 <- diss(dataKLM_long, METHOD = "MINDIST.SAX", w=20, alpha=6)

# "SumSquare","Silhouette" and "Dunn" by number of clusters ####
nclu <- c(2,3,4,5,6,7,8,9,10)
comparecls2 <- as.data.frame(matrix(NA, ncol=4,nrow=length(nclu)))
comparecls2[,1] <- nclu
colnames(comparecls2) <- c("Nclust","SumSquare","Silhouette","Dunn")
for(k in 1:length(nclu)){
  evalcls_20c <- cluster.stats(distSAX_w20a6,cutree(hclust(distSAX_w20a6),nclu[k] ))
  comparecls2[k,2] <- evalcls_20c$within.cluster.ss
  comparecls2[k,3] <- evalcls_20c$avg.silwidth
  comparecls2[k,4] <- evalcls_20c$dunn2
}

# Number of clusters ####

# par(mfrow = c(1,1))  
# plot(comparecls2[,1],comparecls2[,2])
# plot(comparecls2[,1],comparecls2[,3])
# plot(comparecls2[,1],comparecls2[,4])


# Figure 1: Dendogramm ####

colorarticle <- c("#E8BF28","#F8D488","#FF4901","#B6EDD9","#318CE7",
                  "#80D0D0")
distSAX <- distSAX_w20a6
hclustSAX <- hclust(distSAX)

d <- mutate(dataKLM_long, Time=ymd(colnames(dataKLM_wide)[-1]))
dg <- gather(d,"code","freq",1:ncol(dataKLM_long))
catClust <- cutree(hclustSAX,4)
catClust <- data.frame(catClust)
catClust$code <- rownames(catClust)
d <- left_join(dg,catClust, by="code")
d$catClust <- as.factor(d$catClust)
levels(d$catClust) <- paste0("Cluster ",levels(d$catClust))

# Clusters Description ####
cluster.stats(distSAX_w20a6,cutree(hclust(distSAX_w20a6),4 ))$cluster.size
cluster.stats(distSAX_w20a6,cutree(hclust(distSAX_w20a6),4 ))$clus.avg.silwidths
unique(d[which(d$code %in% c("A15","A37","B05","B18","E43","E44","E66","E46","C54","C55","C53","C67")),c("code","catClust")])
datclusterscat <- unique(d[,c("code","catClust")])
# Figure 2 ####
plot_50b10 <- ggplot(as.data.frame(d), aes(x=Time,y=freq)) + 
  geom_path(aes(group=code,col=as.character(catClust)),alpha=0.15) + 
  facet_wrap(~ catClust,ncol = 2, scales="free") + theme_classic() + 
  ylab("Normalized Moving Mean Frequency")+guides(colour=FALSE) +
  ylim(c(-4,4)) + geom_hline(yintercept = 0,alpha=1)
ggsave(plot_50b10,filename = paste0(rept2,"Figure2.png"),device = "png",dpi = 300)

#009681 Turquoire Cluster 3
#767F00 verdatre Cluster 2
#7866D8 Mauve Cluster 4
#CC476B Rouge Cluster 1
# Cluster 1 : "#F8766D"
# Cluster 2 : "#7CAE00"
# Cluster 3 : "#00BFC4"
# Cluster 4 : "#C77CFF"
# Cluster 2 Cluster 4 Cluster 1 Cluster 3
colorarticle <- c("#F8766D","#7CAE00","#00BFC4","#C77CFF")
# colorarticle <- colorarticle[c(2,4,3,1)]
evalcls_20c <- cluster.stats(distSAX,cutree(hclust(distSAX),4 ))
listecodesclusters <- cutree(hclust(distSAX),4)
listecodesclusters[names(listecodesclusters)=="A37"]
#
dend <- hclust(distSAX) %>% as.dendrogram %>% set("branches_k_color", k=4) %>% set("branches_lwd", 1) %>%
  set("labels_colors","white") %>% set("labels_cex", c(.1)) %>%
  set("leaves_pch", 8) %>% set("leaves_col", "black")
dend <- color_branches(dend,k=4, col = colorarticle[c(2,4,1,3)])
dend <- remove_leaves_nodePar(dend)
dend <- remove_nodes_nodePar(dend)
dend <- set(dend,"labels",value = "")
# plot(dend)
# unique(d[which(d$code %in% c("G32")),])
# labels_cex(dend) <- 0
# labels_cex(dend)[300] <- 2 
ggd1 <- ggplot(dend) + annotate("text", x = 115, y = -0.5, label = "Cluster 2",fontface =2)+
  annotate("text", x = 600, y = -0.5, label = "Cluster 4",fontface =2)+
  annotate("text", x = 1100, y = -0.5, label = "Cluster 1",fontface =2)+
  annotate("text", x = 1450, y = -0.5, label = "Cluster 3",fontface =2)
ggsave(ggd1,filename = paste0(rept2,"Figure1.png"),device = "png",dpi = 800)


# Figure 3: Mean curves ####

datash <- as.data.frame(d)
datash$freq
datash <- datash %>% select(catClust,Time,freq)%>% group_by(catClust,Time) %>% mutate(freq2=mean(freq)) %>%  ungroup()
datash <- unique(datash[,c("catClust", "Time", "freq2")])
ggplot_mean <- ggplot(datash, aes(x=Time,y=freq2, group=catClust)) + geom_path() +
  facet_wrap(~ catClust,ncol = 2, scales="free") + theme_classic() + 
  ylab("Normalized Moving Mean Frequency")+guides(colour=FALSE) +
  ylim(c(-2,2)) + geom_hline(yintercept = 0,alpha=1, lty=2)
ggsave(ggplot_mean ,filename = paste0(rept2,"Figure3.png"),device = "png",dpi = 300)

# Figures 4 & 5 ####

options(scipen=999)

fun_graph <- function(x,y){
  code <- x
  codelab <- y
  dtt2 <- data_all %>% filter(diag2==code)
  dtt2$freqs <- NA
  i <- 6
  for(j in (1+i):(length(unique(dtt2$temps_num))-i)){
    temp <- 0
    for(w in -i:i){
      temp <- temp + dtt2[j+w,"freq"]
    }
    dtt2[j,"freqs"] <- temp/(2*i+1)
  }
  obj <- ggplot(dtt2) + geom_point(aes(x=as.Date(temps_num,origin="1970-01-01"),y=freq)) + 
    theme_bw() + 
    geom_point(aes(x=as.Date(temps_num,origin="1970-01-01"),y=freqs),shape=3) + 
    ylab("Frequency") + xlab(codelab) + xlim(c(as.Date(13879,origin="1970-01-01"),as.Date(17167,origin="1970-01-01"))) +
    ylim(c(0,max(dtt2$freqs))) + theme(axis.text=element_text(size=12),
                                       axis.title=element_text(size=12,face="bold"))
  return(obj)
}
fun_graph(x = "E63",y="B05: Measles")

obj_A37 <- fun_graph(x = "A37",y="A37")
obj_B05 <- fun_graph(x = "B05",y="B05")
obj_A53 <- fun_graph(x = "A53",y="A53")
obj_A54 <- fun_graph(x = "A54",y="A54")

obj_C53 <- fun_graph(x = "C53",y="C53")
obj_C54 <- fun_graph(x = "C54",y="C54")
obj_C55 <- fun_graph(x = "C55",y="C55")

obj_E43 <- fun_graph(x = "E43",y="E43")
obj_E44 <- fun_graph(x = "E44",y="E44")
obj_E46 <- fun_graph(x = "E46",y="E46")

infectious <- plot_grid(obj_B05,obj_A53,obj_A37,obj_A54, labels=c("A", "B","C","D"), ncol = 2, nrow = 2)

extrinsic <- plot_grid(obj_C53,obj_C54,obj_C55,obj_E43,obj_E44,obj_E46, labels=c("A", "B","C","D","E","F"), ncol = 2, nrow = 3)

ggsave(infectious,filename = paste0(rept2,"Figure4.png"),device = "png",dpi = 300)
ggsave(extrinsic,filename = paste0(rept2,"Figure5.png"),device = "png",dpi = 300)
