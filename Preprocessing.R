library(tidyverse)
library(flowCore)
library(flowClust)
library(flowSet)
library(ggcyto)
library(gridExtra)
library(knitr)


# Example from flowCore package: https://github.com/RGLab/flowCore/blob/master/examples/example1.r
file.name <- system.file("extdata","0877408774.B08", package="flowCore")
b08 <- read.FCS(file.name, transformation="scale")

dim(b08)
names(b08)

range(b08@exprs[,"SSC-H"])


library(flowViz)
flowPlot(b08,plotParameters=c("FSC-H","SSC-H"))
flowPlot(b08,plotParameters=c("FL1-H","FL2-H"))

# Compensation example: https://www.bioconductor.org/packages/devel/bioc/vignettes/flowCore/inst/doc/HowTo-flowCore.pdf
library(flowStats)
fcs.dir <- system.file("extdata", "compdata", "data",package="flowCore")

frames <- lapply(dir(fcs.dir, full.names=TRUE), read.FCS)
names(frames) <- c("UNSTAINED", "FL1-H", "FL2-H", "FL4-H", "FL3-H") 
frames <- as(frames, "flowSet")
comp <- spillover(frames, unstained="UNSTAINED", patt = "-H",
                  fsc = "FSC-H", ssc = "SSC-H",
                  stain_match = "ordered")
comp

fs <- read.flowSet(path=system.file("extdata", "compdata", "data", package="flowCore"), name.keyword="SAMPLE ID")
fs.exp <- fs[[1]]
rectGate <- rectangleGate(filterId="Fluorescence Region", "FL1-H"=c(0, 12), "FL2-H"=c(0, 12))
result = filter(fs.exp,rectGate) 
result

# Example from FACS by Hellmuth: https://jchellmuth.com/posts/FACS-with-R/
fs <- read.flowSet(path = "~/Desktop/2020-07-08-FACS-data//",pattern = ".fcs",alter.names = T)
pData(fs)[1:3,]  #get sample names
colnames(fs)
autoplot(fs,"FSC.H")
range(fs[[1]]@exprs[,"FSC.A"])

colnames(fs)[colnames(fs)=="FITC.A"] <- "GFP"
colnames(fs)[colnames(fs)=="Pacific.Blue.A"] <- "BFP"

# Try compensation on different flowflames, spillover matrix is identity matrix, results stay same
fs1 <- fs[[1]]
spillover(fs1)
fs1_comp <- compensate(fs1,spillover(fs1)$SPILL)
fs1_comp

# Cleaning
library(flowAI)
fs1_comp_clean <- flow_auto_qc(fs1_comp)

# Transformation
trans <- estimateLogicle(fs1_comp_clean,colnames(fs1_comp_clean[,7:9]))
fs1_comp_clean_trans <- transform(fs1_comp_clean,trans)

autoplot(fs1_comp_clean_trans)
autoplot(fs1_comp_clean_trans,x="FITC.A",y="Pacific.Blue.A",bins=256)
autoplot(fs1_comp_clean_trans,x="Time",y="FSC.A")


library(flowWorkspace)
gs <- GatingSet(fs) #set gates

ggcyto(gs[[1]],aes(x=FSC.A,y=SSC.A),subset="root")+
  geom_hex(bins=200)+
  ggcyto_par_set(limits="instrument")


# Own data v2.0
fcs.direction <- "Z-Project Bacillus_pretest/fcs"
z_fs <- read.flowSet(path = fcs.direction,alter.names = T)
z_gate <- read.FCS("Z-Project Bacillus_pretest/Gating.fcs",alter.names = T)
markernames(z_fs[[1]])

# Plot of FSC vs SSC
ggplot(z_fs[[1]],aes(x=`PMT.1`,y=`PMT.2`))+
  geom_hex(bins=100)+
  theme_bw()+
  labs(x="FSC",y="SSC",title = "FSC vs SSC Scatter Plot",subtitle = "Inner_zone_DAPI")+
  xlim(0,2500)+ylim(0,2500)


# Compensation
y_comp <- compensate(z_fs[[1]],spillover(z_fs[[1]])$SPILL)
y_comp

ggplot(y_comp,aes(x=`PMT.1`,y=`PMT.2`))+
  geom_hex(bins=100)+
  theme_bw()+
  labs(x="FSC",y="SSC",title = "FSC vs SSC after compensation",subtitle = "Inner_zone_DAPI")+
  xlim(0,2500)+ylim(0,2500)


# Cleaning
library(flowAI)
y_comp_clean <- flow_auto_qc(y_comp)

# Transformation
trans <- estimateLogicle(y_comp_clean,colnames(y_comp_clean[,c(11,13)]))
y_comp_clean_trans <- transform(y_comp_clean,trans)
autoplot(y_comp_clean_trans)

# Using arcsinh
channels_y <- c("PMT.1","PMT.2")
trans_y <- transform(y_comp_clean,channels=channels_y,arcsinhTrans)


ggplot(y_comp_clean_trans,aes(x=`PMT.1`,y=`PMT.2`))+
  geom_hex(bins=100)+
  theme_bw()+
  labs(x="FSC",y="SSC",title = "FSC vs SSC after transformation",subtitle = "Inner_zone_DAPI")

# Gating
y_gate <- flowClust(y_comp_clean_trans,varNames=c("PMT.1","PMT.2"), K=2)


# Compensation + Transformation for flowSet
comp_matrix <- spillover(z_fs[[1]])$SPILL
z_fs_comp <- compensate(z_fs,comp_matrix)
colnames(comp_matrix) <- c("PMT.1","PMT.2","PMT.3","PMT.4","PMT.5","PMT.6",
                           "PMT.7","PMT.8","PMT.9","PMT.10","PMT.11","PMT.12",
                           "PMT.13","PMT.14","ADC.15","ADC.16")
# z_fs_comp_clean <- flow_auto_qc(z_fs_comp) #这一步耗时很久 少跑
z_fs_comp_clean <- z_fs_comp
trans_matrix <- estimateLogicle(z_fs_comp_clean[[1]],colnames(z_fs_comp_clean[,11:47]))
z_fs_comp_clean_trans <- transform(z_fs_comp_clean,trans_matrix)
autoplot(z_fs_comp_clean_trans[[1]])

# check this one
autoplot(transform(z_fs_comp_clean_trans[[1]],
                   `PMT.1`=log(`PMT.1`),`PMT.2`=log(`PMT.2`)),
         "PMT.1","PMT.2")


# Set gates
z_auto_gs <- GatingSet(z_fs_comp_clean_trans)
nodes <- gs_get_pop_paths(z_auto_gs, path = "auto")
nodes


# Cell gate
z_fs_data <- gs_pop_get_data(z_auto_gs)
z_nonDebris_gate <- fsApply(z_fs_data,function(fr) 
  openCyto::gate_flowclust_2d(fr,xChannel="PMT.1",yChannel="PMT.2",K=3))
gs_pop_add(z_auto_gs,z_nonDebris_gate,parent="root",name="z_nonDebris_gate")
recompute(z_auto_gs)
autoplot(z_auto_gs[[1]],x="PMT.1",y="PMT.2","z_nonDebris_gate",bins=256)




# Singlet gate
fs_data <- gs_pop_get_data(z_auto_gs,"z_nonDebris_gate") #get parent data
singlet_gate <- fsApply(z_fs_data, function(fr)
  flowStats::gate_singlet(fr,area="PMT.1",height="PMT.2"))
gs_pop_add(z_auto_gs,singlet_gate,parent="z_nonDebris_gate",name="singlet")
recompute(z_auto_gs)
autoplot(z_auto_gs[[1]],x="PMT.1",y="PMT.2","singlet",bins=256)


#
outlier.gate <- rectangleGate(filterId = "-outlier","PMT.1"=c(0,5000),"PMT.2"=c(0,5000))
ggcyto(z_fs[[1]],aes(x=PMT.1,y=PMT.2),subset="root")+
  geom_hex(bins=64)+geom_gate(outlier.gate)

ggcyto(z_auto_gs[[1]],aes(x=PMT.1,y=PMT.2),subset="root")+geom_hex(bins=64)
ggcyto(z_auto_gs[[1]],aes(x=PMT.1,y=PMT.2),subset="singlets")+geom_hex(bins=64)



# Own data v1.0
Inner_zone_DAPI <- read.FCS("Z-Project Bacillus_pretest/fcs/Inner_zone_DAPI.fcs", 
                            transformation = "scale",alter.names = T)
colnames(Inner_zone_DAPI)
autoplot(Inner_zone_DAPI,"PMT.9")
comp_list <- spillover(Inner_zone_DAPI) #On fluorescence part


ggcyto(Inner_zone_DAPI,aes(x=PMT.1,y=PMT.2),subset="root")+
  geom_hex(bins=100)+
  ggcyto_par_set(limits="instrument")+
  ggcyto_par_set(limits=list(x=c(0,0.6),y=c(0,0.6)))

fcs.direction <- "Z-Project Bacillus_pretest/fcs"
z.fs <- read.flowSet(path = fcs.direction,emptyValue = F,alter.names = T)

z.frames <- z.fs
all_cols <- colnames(z.frames)
keep_cols <- all_cols[!(all_cols %in% c("PMT.1","PMT.2","PMT.3","PMT.4"))]
z.frames <- z.frames[,keep_cols]

z.frames1 <- z.fs[1:4]
sampleNames(z.frames1) <- c("UNSTAINED","PMT.1","PMT.2","PMT.3") 
z.frames1 <- as(z.frames1, "flowSet")
comp <- spillover(z.frames1, unstained="UNSTAINED",patt = "PMT.",
                  fsc = "PMT.1", ssc = "PMT.2",
                  stain_match = "ordered")
comp

z.comp_match <- "Z-Project Bacillus_pretest/gate"
writeLines(readLines(z.comp_match))

#z.fs.gate <-read.FCS("Z-Project Bacillus_pretest/Gating.fcs",transformation = F,alter.names = T)
# sample1 <-read.FCS("Z-Project Bacillus_pretest/Gating.fcs")


# Graph
ggcyto(z.fs.gate[[1]],aes(x='PMT.1',y='PMT.2'),subset="root")+
  geom_hex(bins=200)+
  ggcyto_par_set(limits="instrument")+
  ggcyto_par_set(limits=list(x=c(0,3000),y=c(0,3000)))



sdf <- data.frame(name=sampleNames(z.fs)) %>%
  mutate(strain=str_split(name,"-",simplify = T)[,1]) 

phenoData(z.fs)$Filename <- fsApply(z.fs,keyword, "$FIL")
pData(phenoData(z.fs))

p1<- autoplot(z.fs[sampleNames(z.fs) %in% c("Inner_zone_DAPI.fcs","Inner_zone_filter_DAPI.fcs")],"PMT.1.Area")+
  scale_x_flowjo_biexp()
p1

# Compensation + Transformation
comp <- fsApply(z.fs, function(x) spillover(x)[[1]], simplify=FALSE) 
fs_comp <- compensate(z.fs, comp)

transList <- estimateLogicle(z.fs[[1]], colnames(z.fs[[1]]))
p21 <- autoplot(transform(z.fs[[1]], transList), "PMT.1","PMT.2") +
  ggtitle("Before")
p22 <- autoplot(transform(fs_comp[[1]], transList), "PMT.1","PMT.2") +
  ggtitle("After")
grid.arrange(as.ggplot(p21), as.ggplot(p22), ncol = 2)

FSC.H=c(1e5,1e6)
SSC.H=c(1e2,1e6)

ot<- rectangleGate(filterId = "-outlier","FSC.H"=c(1e5,1e6),"SSC.H"=c(1e2,1e6))
ggcyto(z.fs[name %in% c("Inner_zone_DAPI.fcs")],aes(x=FSC.H,y=SSC.H))+geom_gate(ot)