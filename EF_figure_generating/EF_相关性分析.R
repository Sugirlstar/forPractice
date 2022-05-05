# 分季节计算相关性-出图
#####
library(sp)
library(maptools)
library(rgdal)
library(raster)
library(tcltk)
library(RColorBrewer)
library(trend)
library(lattice)
library(VineCopula)
library(copula)

#01  确定各季节位置
#####
mon = substring(dimnames(TX)[[3]], 1, 6) #提取年月
spring=c("03","04","05")
summer=c("06","07","08")
fall=c("09","10","11")
winter=c("12","01","02")
springix=which( substring(mon, 5, 6) %in% spring ) #各季节位置
summerix=which( substring(mon, 5, 6) %in% summer ) #各季节位置
fallix=which( substring(mon, 5, 6) %in% fall ) #各季节位置
winterix=which( substring(mon, 5, 6) %in% winter ) #各季节位置

springcor = array( dim = c(xx, yy, 2) )
summercor = array( dim = c(xx, yy, 2) )
fallcor = array( dim = c(xx, yy, 2) )
wintercor = array( dim = c(xx, yy, 2) )

springptype=summerptype=fallptype=winterptype = array(dim=c(xx,yy))
#####

  #02 分季节计算相关性
  #####
  for(x in 1:xx)
    for(y in 1:yy)
      if (TRUE %in% (SPIarray[x, y, ] >= -999) & TRUE %in% (TX[x, y, ] >= -999))
      {
        springcor[x, y, 1] = cor.test(SPIarray[x,y,springix], TSI[x,y,springix])$estimate 
        springcor[x, y, 2] = cor.test(SPIarray[x,y,springix], TSI[x,y,springix])$p.value 
        springptype[x,y]=classtype(springcor[x, y, 2])
        
        summercor[x,y,1] = cor.test(SPIarray[x,y,summerix], TSI[x,y,summerix])$estimate 
        summercor[x,y,2] = cor.test(SPIarray[x,y,summerix], TSI[x,y,summerix])$p.value
        summerptype[x,y]=classtype(summercor[x, y, 2])
        
        fallcor[x,y,1] = cor.test(SPIarray[x,y,fallix], TSI[x,y,fallix])$estimate 
        fallcor[x,y,2] = cor.test(SPIarray[x,y,fallix], TSI[x,y,fallix])$p.value
        fallptype[x,y] = classtype(fallcor[x, y, 2])
        
        wintercor[x,y,1] = cor.test(SPIarray[x,y,winterix], TSI[x,y,winterix])$estimate 
        wintercor[x,y,2] = cor.test(SPIarray[x,y,winterix], TSI[x,y,winterix])$p.value
        winterptype[x,y] = classtype(wintercor[x, y, 2])
      }

   #####

#03 出图
#####
springcorr = raster(springcor[,,1]) #转为r层
summercorr = raster(summercor[,,1])
fallcorr = raster(fallcor[,,1])
wintercorr = raster(wintercor[,,1])
extent(springcorr) = c(72, 136, 18, 54) #确定范围
extent(summercorr) = c(72, 136, 18, 54) #确定范围
extent(fallcorr) = c(72, 136, 18, 54) #确定范围
extent(wintercorr) = c(72, 136, 18, 54) #确定范围


setwd("C:/Users/Yaoying/Desktop/2108论文收尾/211001定稿/EFfig")
colr11 = rev(colorRampPalette(brewer.pal(9,"BrBG")[c(1:4,6:9)])(100)) 
colr1 = cf(colr11,whitesite=colr11[50],c(-0.45,0.15,0),z0=0) 

#03-01 spring
#####
#画显著性水平
#位置点计算
yind=xind=pfre2ind=NULL
for(x in 1:xx)
  for(y in 1:yy)
  {
    yind=c(yind,71.75+0.5*y) #经度
    xind=c(xind,54.25-0.5*x) #纬度
    pfre2ind=c(pfre2ind,springptype[x, y])
  }


png(
  file = "spring_cor.png",
  width = 2000,
  height = 1300,
  res = 72 * 3
)
par(mar=c(1,1,1,2))
plot(
  springcorr,
  ylim=c(17,55),xlim=c(72,136),
  xaxs="i",yaxs="i",
  xaxt="n",yaxt="n",
  xlab="",ylab="",
  mgp=c(0,0,0),
  col=colr1,
  legend.width=1,
  zlim=c(-0.45, 0.15),
  axis.args=list(at=seq(-0.45, 0.15, 0.1),
                 labels=seq(-0.45, 0.15, 0.1),
                 cex.axis=1)
)

xL = seq(72, 136, 5)
xlabel = rep("",length(xL)) #设置坐标轴
yL = seq(18, 54, 5)
ylabel = rep("",length(yL))

axis(1, xL, xlabel)
axis(2, yL, ylabel,las=2)
lines(Chinasp)
lines(provincesp)

points(yind,xind,pch=pfre2ind,cex=0.5) #点显著性
points(124,19,pch=4,cex=1.5) #显著性水平标注
#text(130,19.1,"α=0.05",cex=1) #显著性水平标注

dev.off()
#####

#03-02 summer
#####
#画显著性水平
#位置点计算
yind=xind=pfre2ind=NULL
for(x in 1:xx)
  for(y in 1:yy)
  {
    yind=c(yind,71.75+0.5*y) #经度
    xind=c(xind,54.25-0.5*x) #纬度
    pfre2ind=c(pfre2ind,summerptype[x, y])
  }


png(
  file = "summer_cor.png",
  width = 2000,
  height = 1300,
  res = 72 * 3
)
par(mar=c(1,1,1,2))
plot(
  summercorr,
  ylim=c(17,55),xlim=c(72,136),
  xaxs="i",yaxs="i",
  xaxt="n",yaxt="n",
  xlab="",ylab="",
  mgp=c(0,0,0),
  col=colr1,
  legend.width=1,
  zlim=c(-0.45, 0.15),
  axis.args=list(at=seq(-0.45, 0.15, 0.1),
                 labels=seq(-0.45, 0.15, 0.1),
                 cex.axis=1)
)

xL = seq(72, 136, 5)
xlabel = rep("",length(xL)) #设置坐标轴
yL = seq(18, 54, 5)
ylabel = rep("",length(yL))

axis(1, xL, xlabel)
axis(2, yL, ylabel,las=2)
lines(Chinasp)
lines(provincesp)

points(yind,xind,pch=pfre2ind,cex=0.5) #点显著性
points(124,19,pch=4,cex=1.5) #显著性水平标注
#text(130,19.1,"α=0.05",cex=1) #显著性水平标注

dev.off()
#####

#03-03 fall
#####
#画显著性水平
#位置点计算
yind=xind=pfre2ind=NULL
for(x in 1:xx)
  for(y in 1:yy)
  {
    yind=c(yind,71.75+0.5*y) #经度
    xind=c(xind,54.25-0.5*x) #纬度
    pfre2ind=c(pfre2ind,fallptype[x, y])
  }


png(
  file = "fall_cor.png",
  width = 2000,
  height = 1300,
  res = 72 * 3
)
par(mar=c(1,1,1,2))
plot(
  fallcorr,
  ylim=c(17,55),xlim=c(72,136),
  xaxs="i",yaxs="i",
  xaxt="n",yaxt="n",
  xlab="",ylab="",
  mgp=c(0,0,0),
  col=colr1,
  legend.width=1,
  zlim=c(-0.45, 0.15),
  axis.args=list(at=seq(-0.45, 0.15, 0.1),
                 labels=seq(-0.45, 0.15, 0.1),
                 cex.axis=1)
)

xL = seq(72, 136, 5)
xlabel = rep("",length(xL)) #设置坐标轴
yL = seq(18, 54, 5)
ylabel = rep("",length(yL))

axis(1, xL, xlabel)
axis(2, yL, ylabel,las=2)
lines(Chinasp)
lines(provincesp)

points(yind,xind,pch=pfre2ind,cex=0.5) #点显著性
points(124,19,pch=4,cex=1.5) #显著性水平标注
#text(130,19.1,"α=0.05",cex=1) #显著性水平标注

dev.off()
#####

#03-03 winter
#####
#画显著性水平
#位置点计算
yind=xind=pfre2ind=NULL
for(x in 1:xx)
  for(y in 1:yy)
  {
    yind=c(yind,71.75+0.5*y) #经度
    xind=c(xind,54.25-0.5*x) #纬度
    pfre2ind=c(pfre2ind,winterptype[x, y])
  }


png(
  file = "winter_cor.png",
  width = 2000,
  height = 1300,
  res = 72 * 3
)
par(mar=c(1,1,1,2))
plot(
  wintercorr,
  ylim=c(17,55),xlim=c(72,136),
  xaxs="i",yaxs="i",
  xaxt="n",yaxt="n",
  xlab="",ylab="",
  mgp=c(0,0,0),
  col=colr1,
  legend.width=1,
  zlim=c(-0.45, 0.15),
  axis.args=list(at=seq(-0.45, 0.15, 0.1),
                 labels=seq(-0.45, 0.15, 0.1),
                 cex.axis=1)
)

xL = seq(72, 136, 5)
xlabel = rep("",length(xL)) #设置坐标轴
yL = seq(18, 54, 5)
ylabel = rep("",length(yL))

axis(1, xL, xlabel)
axis(2, yL, ylabel,las=2)
lines(Chinasp)
lines(provincesp)

points(yind,xind,pch=pfre2ind,cex=0.5) #点显著性
points(124,19,pch=4,cex=1.5) #显著性水平标注
#text(130,19.1,"α=0.05",cex=1) #显著性水平标注

dev.off()
#####


#####
 length(which(springcor[,,1] <= (-0.2))) / length(which(springcor[,,1]>=-999))
 length(which(summercor[,,1] <= (-0.2))) / length(which(summercor[,,1]>=-999))

#####
