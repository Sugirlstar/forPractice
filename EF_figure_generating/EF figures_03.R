#EF文件出图
#3 Fig3-Spatial_variation_3r2c
#####
# seperate in to 3 workspace

#3-1 DR
# Open the workspace:  .../干旱计算全流程.RData
  #00 Background setting
#####
  load("C:/Users/Yaoying/Desktop/2108论文收尾/211001定稿/SPI_TX.RData") #SPIarray
  load("C:/Users/Yaoying/Desktop/2108论文收尾/211001定稿/r层.RData") #r层 吉祥物及水印作用:)
  #设置图片存储位置及名称
  setwd("C:/Users/Yaoying/Desktop/2108论文收尾/211001定稿/EFfig")
  #加载必要包
  library(tcltk)
  library(raster)
  library(sp)
  library(maptools)
  library(rgdal)
  library(RColorBrewer) #载入调色盘
  library(trend)
  library(lattice)
  #读取shp文件
  Chinasp = readShapePoly("C:/Users/Yaoying/Desktop/hyj/china/China") #中国边界
  provincesp = readShapePoly("C:/Users/Yaoying/Desktop/hyj/china/PROVINCE_region") #省份边界
  #基本参数
  xx=72
  yy=128
  days=length(SPIarray[1,1,])
  st=90 #spi累积降水
  da = 10 #最低连续时间
  da2 = 5
  thr1 = -1 #第一阈值
  thr2 = 0 #第二阈值
  s=1961
  e=2020
  LL=e-s+1 #总年数
#####
      #01 cuculate linear trend
  #####
      objx = DRfre
      ck = DRFREh
      #计算突变点
      rst=bp1(s,e,ck)
      #计算两段变化率
      bp=rst[1]-s+1 #突变点序号
      tend2=array(dim=c(xx,yy))
      pfre2=array(dim=c(xx,yy))
      #第二段
      for(x in 1:xx)
        for(y in 1:yy)
          if (TRUE %in% (objx[x,y,]>=0))
          {
            tend2[x,y]=tendency(objx[x,y,(bp+1):(e-s+1)])[[1]]
            pfre2[x,y]=classtype(tendency(objx[x,y,(bp+1):(e-s+1)])[[2]])
          }
      #放入r层
      tend2r = raster(tend2) #转为r层
      extent(tend2r) = c(72, 136, 18, 54) #确定范围
      crs(tend2r) = crs(r) #确定投影
      #画显著性水平
      #位置点计算
      yind=xind=pfre1ind=pfre2ind=NULL
      for(x in 1:xx)
        for(y in 1:yy)
        {
          yind=c(yind,71.75+0.5*y) #经度
          xind=c(xind,54.25-0.5*x) #纬度
          pfre2ind=c(pfre2ind,pfre2[x,y])
        }
  #####
          #02 Putout_linrar trend
      #####
          max(tend2,na.rm=T)
          min(tend2,na.rm=T)
          
          colr1=rev(colorRampPalette(brewer.pal(9,"GnBu")[1:9])(150)) 
          colr2=colorRampPalette(brewer.pal(9,"OrRd")[1:9])(150)
          #tend2col = c(colr1,colr2)
          tend2col = cf(c(colr1,colr2),whitesite="#FFF7EC",
                        c(-0.060, 0.08,0),z0=0) 
      
          png(
            file = "fig3_tend_drought.png",
            width = 2000,
            height = 1300,
            res = 72 * 3
          )
          par(mar=c(1,1,1,2))
          plot(
            tend2r,
            ylim=c(17,55),xlim=c(72,136),
            xaxs="i",yaxs="i",
            xaxt="n",yaxt="n",
            xlab="",ylab="",
            mgp=c(0,0,0),
            col=tend2col,
            legend.width=1,
            zlim=c(-0.060, 0.08),
            axis.args=list(at=seq(-0.06, 0.08, 0.02),
                           labels=seq(-0.06, 0.08, 0.02),
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
          
          #Plus-01 cuculate linear trend
          #####
          objx = DRfre
          ck = DRFREh
          #计算突变点
          rst=bp1(s,e,ck)
          #计算两段变化率
          bp=rst[1]-s+1 #突变点序号
          tend2=array(dim=c(xx,yy))
          pfre2=array(dim=c(xx,yy))
          #第二段
          for(x in 1:xx)
            for(y in 1:yy)
              if (TRUE %in% (objx[x,y,]>=0))
              {
                tend2[x,y]=tendency(objx[x,y,1:bp])[[1]]
                pfre2[x,y]=classtype(tendency(objx[x,y,1:bp])[[2]])
              }
          #放入r层
          tend2r = raster(tend2) #转为r层
          extent(tend2r) = c(72, 136, 18, 54) #确定范围
          crs(tend2r) = crs(r) #确定投影
          #画显著性水平
          #位置点计算
          yind=xind=pfre1ind=pfre2ind=NULL
          for(x in 1:xx)
            for(y in 1:yy)
            {
              yind=c(yind,71.75+0.5*y) #经度
              xind=c(xind,54.25-0.5*x) #纬度
              pfre2ind=c(pfre2ind,pfre2[x,y])
            }
          #####
          
          #Plus-02 Putout_linrar trend
          #####
          max(tend2,na.rm=T)
          min(tend2,na.rm=T)
          crange = c(-0.120, 0.12)
          cgap = 0.04
          
          colr1=rev(colorRampPalette(brewer.pal(9,"GnBu")[1:9])(150)) 
          colr2=colorRampPalette(brewer.pal(9,"OrRd")[1:9])(150)
          #tend2col = c(colr1,colr2)
          tend2col = cf(c(colr1,colr2),whitesite="#FFF7EC",
                        c(crange,0),z0=0) 
          
          png(
            file = "fig3_tend1_drought.png",
            width = 2000,
            height = 1300,
            res = 72 * 3
          )
          par(mar=c(1,1,1,2))
          plot(
            tend2r,
            ylim=c(17,55),xlim=c(72,136),
            xaxs="i",yaxs="i",
            xaxt="n",yaxt="n",
            xlab="",ylab="",
            mgp=c(0,0,0),
            col=tend2col,
            legend.width=1,
            zlim=crange,
            axis.args=list(at=seq(crange[1],crange[2], cgap),
                           labels=seq(crange[1],crange[2], cgap),
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

              #03 caculate M-K trend
          #####
          zz = MK.raster(objx,"year") #逐年
          z <- raster(zz[,,1]) #转为raster层
          extent(z)<- c(72, 136, 18, 54) #确定范围
          res(z)<-res(r) #设分辩率
          projection(z)=projection(r)
          sigz=array(dim=c(xx,yy))
          for(x in 1:xx)
            for(y in 1:yy)
              if (TRUE %in% (zz[x,y,]>=0))
                sigz[x,y]=classtype(zz[x,y,2])
          #确定位置点
          yind=xind=sigind=NULL
          for(x in 1:xx)
            for(y in 1:yy)
            {
              yind=c(yind,71.75+0.5*y) #经度
              xind=c(xind,54.25-0.5*x) #纬度
              sigind=c(sigind,sigz[x,y])
            }
          
          #####
          
          
          #04 Putout_M-K trend
          #####
          #画图
          z # -0.6,0.5
          
          colr=rev(colorRampPalette(brewer.pal(11, "RdYlBu"))(100)) #颜色
          MKcol =  cf(colr,whitesite=colr[50],c(-0.5,0.5,0),z0=0) 
          
          png(
            file = "fig3_MK_drought.png",
            width = 2000,
            height = 1300,
            res = 72 * 3
          )
          par(mar=c(1,1,1,2))
          plot(
            z,
            ylim=c(17,55),xlim=c(72,136),
            xaxs="i",yaxs="i",
            xaxt="n",yaxt="n",
            xlab="",ylab="",
            mgp=c(0,0,0),
            col=MKcol,
            legend.width=1,
            zlim=c(-0.5,0.5),
            axis.args=list(at=seq(-0.5, 0.5, 0.2),
                           labels=seq(-0.5, 0.5, 0.2),
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
          
          points(yind,xind,pch=sigind,cex=0.5) #点显著性
          points(124,19,pch=4,cex=1.5) #显著性水平标注
          #text(130,19.1,"α=0.05",cex=1.2) #显著性水平标注
          
          dev.off()
          #####
          
#3-2 HW
# Open the workspace:  .../热浪计算全流程.RData
#00 Background setting
#####
load("C:/Users/Yaoying/Desktop/2108论文收尾/211001定稿/SPI_TX.RData") #SPIarray
load("C:/Users/Yaoying/Desktop/2108论文收尾/211001定稿/r层.RData") #r层 吉祥物及水印作用:)
#设置图片存储位置及名称
setwd("C:/Users/Yaoying/Desktop/2108论文收尾/211001定稿/EFfig")
#加载必要包
library(tcltk)
library(raster)
library(sp)
library(maptools)
library(rgdal)
library(RColorBrewer) #载入调色盘
library(trend)
library(lattice)
#读取shp文件
Chinasp = readShapePoly("C:/Users/Yaoying/Desktop/hyj/china/China") #中国边界
provincesp = readShapePoly("C:/Users/Yaoying/Desktop/hyj/china/PROVINCE_region") #省份边界
#基本参数
xx=72
yy=128
days=length(SPIarray[1,1,])
st=90 #spi累积降水
da = 10 #最低连续时间
da2 = 5
thr1 = -1 #第一阈值
thr2 = 0 #第二阈值
s=1961
e=2020
LL=e-s+1 #总年数
#####

#01 cuculate linear trend
#####
objx = DRfre
ck = DRFREh
#计算突变点
rst=bp1(s,e,ck)
#计算两段变化率
bp=rst[1]-s+1 #突变点序号
tend2=array(dim=c(xx,yy))
pfre2=array(dim=c(xx,yy))
#第二段
for(x in 1:xx)
  for(y in 1:yy)
    if (TRUE %in% (objx[x,y,]>=0))
    {
      tend2[x,y]=tendency(objx[x,y,(bp+1):(e-s+1)])[[1]]
      pfre2[x,y]=classtype(tendency(objx[x,y,(bp+1):(e-s+1)])[[2]])
    }
#放入r层
tend2r = raster(tend2) #转为r层
extent(tend2r) = c(72, 136, 18, 54) #确定范围
crs(tend2r) = crs(r) #确定投影
#画显著性水平
#位置点计算
yind=xind=pfre1ind=pfre2ind=NULL
for(x in 1:xx)
  for(y in 1:yy)
  {
    yind=c(yind,71.75+0.5*y) #经度
    xind=c(xind,54.25-0.5*x) #纬度
    pfre2ind=c(pfre2ind,pfre2[x,y])
  }
#####

    #02 Putout_linrar trend
    #####
    max(tend2,na.rm=T)
    min(tend2,na.rm=T)
    
    colr1=rev(colorRampPalette(brewer.pal(9,"GnBu")[1:9])(150)) 
    colr2=colorRampPalette(brewer.pal(9,"OrRd")[1:9])(150)
    #tend2col = c(colr1,colr2)
    tend2col = cf(c(colr1,colr2),whitesite="#FFF7EC",
                  c( -0.08,0.27,0),z0=0) 
    
    png(
      file = "fig3_tend_heatwave.png",
      width = 2000,
      height = 1300,
      res = 72 * 3
    )
    par(mar=c(1,1,1,2))
    plot(
      tend2r,
      ylim=c(17,55),xlim=c(72,136),
      xaxs="i",yaxs="i",
      xaxt="n",yaxt="n",
      xlab="",ylab="",
      mgp=c(0,0,0),
      col=tend2col,
      legend.width=1,
      zlim=c(-0.08,0.27),
      axis.args=list(at=seq(-0.08,0.27, 0.05),
                     labels=seq(-0.08,0.27, 0.05),
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
    
    #Plus-01 cuculate linear trend
    #####
    objx = DRfre
    ck = DRFREh
    #计算突变点
    rst=bp1(s,e,ck)
    #计算两段变化率
    bp=rst[1]-s+1 #突变点序号
    tend2=array(dim=c(xx,yy))
    pfre2=array(dim=c(xx,yy))
    #第二段
    for(x in 1:xx)
      for(y in 1:yy)
        if (TRUE %in% (objx[x,y,]>=0))
        {
          tend2[x,y]=tendency(objx[x,y,1:bp])[[1]]
          pfre2[x,y]=classtype(tendency(objx[x,y,1:bp])[[2]])
        }
    #放入r层
    tend2r = raster(tend2) #转为r层
    extent(tend2r) = c(72, 136, 18, 54) #确定范围
    crs(tend2r) = crs(r) #确定投影
    #画显著性水平
    #位置点计算
    yind=xind=pfre1ind=pfre2ind=NULL
    for(x in 1:xx)
      for(y in 1:yy)
      {
        yind=c(yind,71.75+0.5*y) #经度
        xind=c(xind,54.25-0.5*x) #纬度
        pfre2ind=c(pfre2ind,pfre2[x,y])
      }
    #####
    
    #Plus-02 Putout_linrar trend
    #####
    max(tend2,na.rm=T)
    min(tend2,na.rm=T)
    crange = c(-1.10, 0.3)
    cgap = 0.2
    seq(crange[1],crange[2],cgap)
    
    colr1=rev(colorRampPalette(brewer.pal(9,"GnBu")[1:9])(150)) 
    colr2=colorRampPalette(brewer.pal(9,"OrRd")[1:9])(150)
    #tend2col = c(colr1,colr2)
    tend2col = cf(c(colr1,colr2),whitesite="#FFF7EC",
                  c( crange,0),z0=0) 
    
    png(
      file = "fig3_tend1_heatwave.png",
      width = 2000,
      height = 1300,
      res = 72 * 3
    )
    par(mar=c(1,1,1,2))
    plot(
      tend2r,
      ylim=c(17,55),xlim=c(72,136),
      xaxs="i",yaxs="i",
      xaxt="n",yaxt="n",
      xlab="",ylab="",
      mgp=c(0,0,0),
      col=tend2col,
      legend.width=1,
      zlim=crange,
      axis.args=list(at=seq(crange[1],crange[2], cgap),
                     labels=seq(crange[1],crange[2], cgap),
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
    

    
    
        #03 caculate M-K trend
        #####
        zz = MK.raster(objx,"year") #逐年
        z <- raster(zz[,,1]) #转为raster层
        extent(z)<- c(72, 136, 18, 54) #确定范围
        res(z)<-res(r) #设分辩率
        projection(z)=projection(r)
        sigz=array(dim=c(xx,yy))
        for(x in 1:xx)
          for(y in 1:yy)
            if (TRUE %in% (zz[x,y,]>=0))
              sigz[x,y]=classtype(zz[x,y,2])
        #确定位置点
        yind=xind=sigind=NULL
        for(x in 1:xx)
          for(y in 1:yy)
          {
            yind=c(yind,71.75+0.5*y) #经度
            xind=c(xind,54.25-0.5*x) #纬度
            sigind=c(sigind,sigz[x,y])
          }
        
        #####
    
            #04 Putout_M-K trend
            #####
            #画图
            z # -0.35,0.65
            
            colr=rev(colorRampPalette(brewer.pal(11, "RdYlBu"))(100)) #颜色
            MKcol =  cf(colr,whitesite=colr[50],c(-0.35,0.65,0),z0=0) 
            
            png(
              file = "fig3_MK_heatwave.png",
              width = 2000,
              height = 1300,
              res = 72 * 3
            )
            par(mar=c(1,1,1,2))
            plot(
              z,
              ylim=c(17,55),xlim=c(72,136),
              xaxs="i",yaxs="i",
              xaxt="n",yaxt="n",
              xlab="",ylab="",
              mgp=c(0,0,0),
              col=MKcol,
              legend.width=1,
              zlim=c(-0.35,0.65),
              axis.args=list(at=seq(-0.35,0.65,0.2),
                             labels=seq(-0.35,0.65, 0.2),
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
            
            points(yind,xind,pch=sigind,cex=0.5) #点显著性
            points(124,19,pch=4,cex=1.5) #显著性水平标注
            #text(130,19.1,"α=0.05",cex=1.2) #显著性水平标注
            
            dev.off()
            #####
    
#3-3 CDH
# Open the workspace:  .../复合计算_final.RData
#00 Background setting
#####
load("C:/Users/Yaoying/Desktop/2108论文收尾/211001定稿/SPI_TX.RData") #SPIarray
load("C:/Users/Yaoying/Desktop/2108论文收尾/211001定稿/r层.RData") #r层 吉祥物及水印作用:)
#设置图片存储位置及名称
setwd("C:/Users/Yaoying/Desktop/2108论文收尾/211001定稿/EFfig")
#加载必要包
library(tcltk)
library(raster)
library(sp)
library(maptools)
library(rgdal)
library(RColorBrewer) #载入调色盘
library(trend)
library(lattice)
#读取shp文件
Chinasp = readShapePoly("C:/Users/Yaoying/Desktop/hyj/china/China") #中国边界
provincesp = readShapePoly("C:/Users/Yaoying/Desktop/hyj/china/PROVINCE_region") #省份边界
#基本参数
xx=72
yy=128
days=length(SPIarray[1,1,])
st=90 #spi累积降水
da = 10 #最低连续时间
da2 = 5
thr1 = -1 #第一阈值
thr2 = 0 #第二阈值
s=1961
e=2020
LL=e-s+1 #总年数
#####

    #01 cuculate linear trend
    #####
    objx = CDHmag
    ck = CDHMAGh
    #计算突变点
    rst=bp1(s,e,ck)
    #计算两段变化率
    bp=rst[1]-s+1 #突变点序号
    tend2=array(dim=c(xx,yy))
    pfre2=array(dim=c(xx,yy))
    #第二段
    for(x in 1:xx)
      for(y in 1:yy)
        if (TRUE %in% (objx[x,y,]>=0))
        {
          tend2[x,y]=tendency(objx[x,y,(bp+1):(e-s+1)])[[1]]
          pfre2[x,y]=classtype(tendency(objx[x,y,(bp+1):(e-s+1)])[[2]])
        }
    #放入r层
    tend2r = raster(tend2) #转为r层
    extent(tend2r) = c(72, 136, 18, 54) #确定范围
    crs(tend2r) = crs(r) #确定投影
    #画显著性水平
    #位置点计算
    yind=xind=pfre1ind=pfre2ind=NULL
    for(x in 1:xx)
      for(y in 1:yy)
      {
        yind=c(yind,71.75+0.5*y) #经度
        xind=c(xind,54.25-0.5*x) #纬度
        pfre2ind=c(pfre2ind,pfre2[x,y])
      }
    #####

        #02 Putout_linrar trend
        #####
        max(tend2,na.rm=T)
        min(tend2,na.rm=T)
        
        colr1=rev(colorRampPalette(brewer.pal(9,"GnBu")[1:9])(150)) 
        colr2=colorRampPalette(brewer.pal(9,"OrRd")[1:9])(150)
        #tend2col = c(colr1,colr2)
        tend2col = cf(c(colr1,colr2),whitesite="#FFF7EC",
                      c(-0.14,0.26,0),z0=0) 
        
        png(
          file = "fig3_tend_CDH.png",
          width = 2000,
          height = 1300,
          res = 72 * 3
        )
        par(mar=c(1,1,1,2))
        plot(
          tend2r,
          ylim=c(17,55),xlim=c(72,136),
          xaxs="i",yaxs="i",
          xaxt="n",yaxt="n",
          xlab="",ylab="",
          mgp=c(0,0,0),
          col=tend2col,
          legend.width=1,
          zlim=c(-0.14,0.26),
          axis.args=list(at=seq(-0.14,0.26,0.08),
                         labels=seq(-0.14,0.26,0.08),
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
        
        #Plus-01 cuculate linear trend
        #####
        objx = CDHmag
        ck = CDHMAGh
        #计算突变点
        rst=bp1(s,e,ck)
        #计算两段变化率
        bp=rst[1]-s+1 #突变点序号
        tend2=array(dim=c(xx,yy))
        pfre2=array(dim=c(xx,yy))
        #第二段
        for(x in 1:xx)
          for(y in 1:yy)
            if (TRUE %in% (objx[x,y,]>=0))
            {
              tend2[x,y]=tendency(objx[x,y,1:bp])[[1]]
              pfre2[x,y]=classtype(tendency(objx[x,y,1:bp])[[2]])
            }
        #放入r层
        tend2r = raster(tend2) #转为r层
        extent(tend2r) = c(72, 136, 18, 54) #确定范围
        crs(tend2r) = crs(r) #确定投影
        #画显著性水平
        #位置点计算
        yind=xind=pfre1ind=pfre2ind=NULL
        for(x in 1:xx)
          for(y in 1:yy)
          {
            yind=c(yind,71.75+0.5*y) #经度
            xind=c(xind,54.25-0.5*x) #纬度
            pfre2ind=c(pfre2ind,pfre2[x,y])
          }
        #####
        
        #Plus-02 Putout_linrar trend
        #####
        max(tend2,na.rm=T)
        min(tend2,na.rm=T)
        crange = c(-12, 0.9)
        cgap = 0.2
        seq(crange[1],crange[2],cgap)
        
        
        colr1=rev(colorRampPalette(brewer.pal(9,"GnBu")[1:9])(150)) 
        colr2=colorRampPalette(brewer.pal(9,"OrRd")[1:9])(150)
        #tend2col = c(colr1,colr2)
        tend2col = cf(c(colr1,colr2),whitesite="#FFF7EC",
                      c(-0.14,0.26,0),z0=0) 
        
        png(
          file = "fig3_tend1_CDH.png",
          width = 2000,
          height = 1300,
          res = 72 * 3
        )
        par(mar=c(1,1,1,2))
        plot(
          tend2r,
          ylim=c(17,55),xlim=c(72,136),
          xaxs="i",yaxs="i",
          xaxt="n",yaxt="n",
          xlab="",ylab="",
          mgp=c(0,0,0),
          col=tend2col,
          legend.width=1,
          zlim=c(-0.14,0.26),
          axis.args=list(at=seq(-0.14,0.26,0.08),
                         labels=seq(-0.14,0.26,0.08),
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
        

            #03 caculate M-K trend
            #####
            zz = MK.raster(objx,"year") #逐年
            z <- raster(zz[,,1]) #转为raster层
            extent(z)<- c(72, 136, 18, 54) #确定范围
            res(z)<-res(r) #设分辩率
            projection(z)=projection(r)
            sigz=array(dim=c(xx,yy))
            for(x in 1:xx)
              for(y in 1:yy)
                if (TRUE %in% (zz[x,y,]>=0))
                  sigz[x,y]=classtype(zz[x,y,2])
            #确定位置点
            yind=xind=sigind=NULL
            for(x in 1:xx)
              for(y in 1:yy)
              {
                yind=c(yind,71.75+0.5*y) #经度
                xind=c(xind,54.25-0.5*x) #纬度
                sigind=c(sigind,sigz[x,y])
              }
            
            #####
        
                #04 Putout_M-K trend
                #####
                #画图
                z # -0.35,0.65
                
                colr=rev(colorRampPalette(brewer.pal(11, "RdYlBu"))(100)) #颜色
                MKcol =  cf(colr,whitesite=colr[50],c(-0.5,0.4,0),z0=0) 
                
                png(
                  file = "fig3_MK_CDH.png",
                  width = 2000,
                  height = 1300,
                  res = 72 * 3
                )
                par(mar=c(1,1,1,2))
                plot(
                  z,
                  ylim=c(17,55),xlim=c(72,136),
                  xaxs="i",yaxs="i",
                  xaxt="n",yaxt="n",
                  xlab="",ylab="",
                  mgp=c(0,0,0),
                  col=MKcol,
                  legend.width=1,
                  zlim=c(-0.5,0.4),
                  axis.args=list(at=seq(-0.5,0.4,0.1),
                                 labels=seq(-0.5,0.4,0.1),
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
                
                points(yind,xind,pch=sigind,cex=0.5) #点显著性
                points(124,19,pch=4,cex=1.5) #显著性水平标注
                #text(130,19.1,"α=0.05",cex=1.2) #显著性水平标注
                
                dev.off()
                #####
        



