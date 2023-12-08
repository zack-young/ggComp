#!/usr/bin/env Rscript
is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1])
if(!is.installed("optparse")){
  warning("Detect package \"optparse\" is not installed in your R enviroment.")
  warning("Trying to install the \"optparse\" package.")
  warning("If failed, please try to install it mannually.")
  install.packages("optparse")
}
if(!is.installed("gplots")){
  warning("Detect package \"gplots\" is not installed in your R enviroment.")
  warning("Trying to install the \"gplots\" package.")
  warning("If failed, please try to install it mannually.")
  install.packages("gplots")
}

## libraries
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
# Arguments
option_list <- list(
  make_option("--sample1", dest = "sample1", default = "",
              help = "sample1 serial number"),
  make_option("--sample1_name", dest = "sample1_name", default = "",
              help = "sample1 name"),
  make_option("--sample2", dest = "sample2", default = "",
              help = "sample2 serial number"),
  make_option("--sample2_name", dest = "sample2_name", default = "",
              help = "sample2 name"),
  make_option(c("-c","--color"), dest = "color", default = "#282658,#32519D,#39BAEC",
              help = "[opt] three colors needed, in order of low,middle,high [default: blue,white,red]"),
  make_option(c("-W","--width"), dest = "figure.width", default = 16,
              help = "[opt] width of figure (inch). [default: 5]"),
  make_option(c("-H","--height"), dest = "figure.height", default = 9,
              help = "[opt] height of figure (inch). [default: 4]"),
  make_option(c("-d","--DEV"), dest = "DEV", default = ".",
              help = "[opt] file direction. [default:]"),
  make_option(c("-D","--DEV2"), dest = "DEV2", default = ".",
              help = "[opt] output file direction. [default:]")
  
)

parser <- OptionParser(usage = "",
                       option_list=option_list, description = ""

)
## check arguments
arguments <- parse_args(parser)
#
figure.width <- arguments$figure.width
figure.height <- arguments$figure.height
#
color = arguments$color
if(color == ""){ 
  color = strsplit("#282658,#32519D,#39BAEC", split = ",", fixed = T)[[1]]
  }else{
    color <- strsplit(arguments$color, split = ",", fixed = T)[[1]]
  }

#
SAMPLE1=arguments$sample1
SAMPLE1_name=arguments$sample1_name
SAMPLE2=arguments$sample2
SAMPLE2_name=arguments$sample2_name
DEV=arguments$DEV
DEV2=arguments$DEV2


options(scipen=200)
files = read.table(paste(DEV,"/plotfile",sep=""), as.is = T, header = F, comment.char = "")

pdf(paste(DEV2,"/",SAMPLE1_name,"_",SAMPLE2_name,".pdf",sep=""), height = figure.height, width = figure.width)
plotChrom <- function(xleft, ybottom, height, len, centro, binsize){
  
  #len <- len/binsize
  #centro <- centro/binsize
  
  # r vertical, r horizontal
  rv <- height/2
  rh <- len/70
  rs <- seq(0,pi,len=100)
  
  # left semi-circle
  lx <- c(xleft, xleft + rh - rh*sin(rs), xleft, xleft)
  ly <- c(ybottom, ybottom + rv - rv*cos(rs), ybottom + height, ybottom)
  polygon(lx, ly, border = "white", col = "white")
  
  # right semi-circle
  rx <- c(xleft + len,xleft + len, xleft + len - rh + rh*sin(rs), xleft + len)
  ry <- c(ybottom, ybottom + height, ybottom + rv +rv*cos(rs), ybottom)
  polygon(rx, ry, border = "white", col = "white")
  
  # top tri
  tx <- c(centro - 0.5*rh, centro, centro + 0.5*rh)
  ty <- c(ybottom + height, ybottom + rv, ybottom + height)
  polygon(tx, ty, border = "white", col = "white")
  
  # bottom tri
  bx <- c(centro - 0.5*rh, centro, centro + 0.5*rh)
  by <- c(ybottom, ybottom + rv, ybottom)
  polygon(bx, by, border = "white", col = "white")
  
  cbx <- c(xleft + rh - rh*sin(rs), centro - 0.5*rh, centro, centro + 0.5*rh, xleft + len - rh + rh*sin(rs), centro + 0.5*rh, centro, centro - 0.5*rh)
  cby <- c(ybottom + rv - rv*cos(rs), ybottom + height, ybottom + rv, ybottom + height, ybottom + rv + rv*cos(rs), ybottom, ybottom + rv, ybottom)
  polygon(cbx, cby, border = "grey", lwd = 3)
}
plot(x=0, type="n", bty="n", yaxt="n",xaxt="n",
     xlab="", ylab="", 
     xlim=c(-200, 1100), ylim=c(-2, nrow(files)+1),
     xaxs="i", yaxs="i", 
     main=paste(SAMPLE1_name,"VS",SAMPLE2_name,sep=" "),
     cex.main = 2.5)
axis(1,at=c(0,100,200,300,400,500,600,700,800),labels = c("0Mbp","","200Mbp","","400Mbp","","600Mbp","","800Mbp"))


del1 <- "#60D6A9"
del2 <- "#9FEE00"
delboth <-  "#007046"
dup1 <- "#FF8B73"
dup2 <- "#FFF073"
dupboth <- "#AA0000"
high <- "#111148"
low <- "#22BBEE"

color_pad <- c(low,high,del1,del2,delboth,dup1,dup2,dupboth)
num = 22
centromere <- c(215,239,173,337,346,268,300,346,241,300,317,184,254,206,189,286,327,211,357,287,340)
wholen <- c(594102056,
            689851870,
            495453186,
            780798557,
            801256715,
            651852609,
            750843639,
            830829764,
            615552423,
            744588157,
            673617499,
            509857067,
            709773743,
            713149757,
            566080677,
            618079260,
            720988478,
            473592718,
            736706236,
            750620385,
            638686055
)

rect(xleft = 900, xright = 930, ybottom = 16, ytop = 18,
     col = color_pad[1], border = color_pad[1]);text(x=c(940), y=c(17),c("SGR"),cex = 1,pos = 4)
rect(xleft = 900, xright = 930, ybottom = 20, ytop = 22,
   col = color_pad[3], border = color_pad[3]);text(x=c(940), y=c(21),c("PHR"),cex = 1,pos = 4)
rect(xleft = 900, xright = 930, ybottom = 22, ytop = 24,
     col = color_pad[4], border = color_pad[4]);text(x=c(940), y=c(23),c(paste(SAMPLE1_name,"deletion",sep="")),cex = 1,pos = 4)
rect(xleft = 900, xright = 930, ybottom = 24, ytop = 26,
     col = color_pad[5], border = color_pad[5]);text(x=c(940), y=c(25),c(paste(SAMPLE2_name,"deletion",sep="")),cex = 1,pos = 4)
rect(xleft = 900, xright = 930, ybottom = 26, ytop = 28,
     col = color_pad[6], border = color_pad[6]);text(x=c(940), y=c(27),c("Both_deletion"),cex = 1,pos = 4)
rect(xleft = 900, xright = 930, ybottom = 28, ytop = 30,
     col = color_pad[7], border = color_pad[7]);text(x=c(940), y=c(29),c(paste(SAMPLE1_name,"duplication",sep="")),cex = 1,pos = 4)
rect(xleft = 900, xright = 930, ybottom = 30, ytop = 32,
     col = color_pad[8], border = color_pad[8]);text(x=c(940), y=c(31),c(paste(SAMPLE2_name,"duplication",sep="")),cex = 1,pos = 4)
rect(xleft = 900, xright = 930, ybottom = 32, ytop = 34,
     col = color_pad[9], border = color_pad[9]);text(x=c(940), y=c(33),c("Both_duplication"),cex = 1,pos = 4)
for(j in 1:nrow(files)){
  if(!is.na(files[j,1])){
    data1 <- read.table(files[j,1], as.is = T, header = F, comment.char = "")
    data1 <- arrange(data1, V2)
    chro <- factor(data1[,4])

    chro1 <- factor(chro, 

                    levels = c("SGR","PHR",
                               paste(SAMPLE1,"deletion_CNV",sep=""),paste(SAMPLE2,"deletion_CNV",sep=""),"deletion_both_CNV",
                               paste(SAMPLE1,"duplication_CNV",sep=""),paste(SAMPLE2,"duplication_CNV",sep=""),"duplication_both_CNV"),

                    labels = c('1','2','3',"4","5",'6',"7","8")

                    )

    
    num <- num - 1 

    for(i in 1:nrow(data1)){
      if (as.numeric(data1[i,2])+1000000 > wholen[num]) {
        xright_use = wholen[num]/1000000
      }else{
        xright_use = as.numeric(data1[i,3])/1000000
      }
      rect(xleft = as.numeric(data1[i,2])/1000000,
           #         #ybottom = j+0.1,
           ybottom = j+0.1,
           ytop = j+0.9,
           xright=xright_use,
           #         #ytop = j+0.9,
           col = color_pad[as.numeric(as.vector(chro1[i]))],
           border = color_pad[as.numeric(as.vector(chro1[i]))] 
      )
    }
    len <- wholen[num]/1000000
    centro <- centromere[num]
    
    text(x=-100, y=j+0.5, files[j,2], cex = 1)

  }
}
num = 22
for(j in 1:nrow(files)){
  if(is.na(files[j,1])){
    j <- j-1
    num <- num - 1
    centro <- centromere[num]
    len <- wholen[num]
    len <- wholen[num]/1000000
    plotChrom(0, j+0.1, 0.8, len, centro, 1000000)


  }
}


dev.off()

