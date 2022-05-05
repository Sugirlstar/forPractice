#EF文件出图
#4 Fig4-SPI-STI
#两个工作空间，分别出两幅图
#打开:复合计算_final

#00 background
#####
load("C:/Users/Yaoying/Desktop/2108论文收尾/211001定稿/droughtFlag.RData") 
load("C:/Users/Yaoying/Desktop/2108论文收尾/211001定稿/heatwaveFlag.RData") 

library("tcltk")
library(raster)
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

#准备函数：找v中位置dot1前面最近的1的位置
findNestOne = function(v,dot1) 
{
  v[which(is.na(v)==TRUE)] = 0 #把NA值置为0，否则循环会报错
  i = dot1 - 1
  outputIndex = i
  
  if(dot1 == 1)
  {
    print("dot1 > 1")
  }else
  {
    while( v[i] != 1 )
    {
      i=i-1
      outputIndex = i
      if( i == 0 )
      {
        outputIndex = NA
        break
      }
    }
    return(outputIndex)
  }
}

######

  #01 计算干旱主导比例
  #####
  leadTime = array( dim=dim(Flag3) ) #设定矩阵放置时间差
  multiEvent = array( dim=dim(Flag3) )
  xx=dim(Flag)[1]
  yy=dim(Flag)[2]
  #计算时间差及重叠事件
  pb <- tkProgressBar(title="进度",label="已完成 %", min=0, max=100)
  for(x in 1:xx)
  {
    for(y in 1:yy)
      if(TRUE %in% (Flag[x,y,] >= -999))
      {
        CDHindex = which( Flag3[x,y,] == 1 ) #起始位置标记：复合事件开始第一天
        
        if( length(CDHindex) == 0 )
          next
        
        for( i in 1:length(CDHindex) )
        {
          dot = CDHindex[i] #复合事件起始点
          
          #找干旱起始点
          if( is.na(FD3[ x,y,dot ]) == TRUE )
            FDdotindex = findNestOne( FD3[x,y,], dot ) else
              if( FD3[ x,y,dot ] ==  0)
                FDdotindex = findNestOne( FD3[x,y,], dot ) else #如果当日为NA或0，则直接找今日前的第一个1
                  FDdotindex = dot
                
                #找热浪起始点
                if( is.na(FH3[ x,y,dot ]) == TRUE )
                  FHdotindex = findNestOne( FH3[x,y,], dot ) else
                    if( FH3[ x,y,dot ] ==  0)
                      FHdotindex = findNestOne( FH3[x,y,], dot ) else #如果当日为NA或0，则直接找今日前的第一个1
                        FHdotindex = dot
                      
                      #计算时间差
                      leadTime[x,y,dot] = FHdotindex - FDdotindex  # >0 - 干旱在前; <0 - 热浪在前; 
                      
        }
        
        #找有没有事件重叠
        if( length(CDHindex)>1 )
        {
          test = CDHindex
          
          for(i in 1:(length(test)-1) )
          {
            intervalDays = (test[i]+1):(test[i+1]-1) #间隔日
            
            t1 = FD[ x,y,intervalDays]
            t2 = FH[ x,y,intervalDays]
            
            if( length(which(t1>0)) == length(t1) ) #如果间隔干旱旗子全为1
            {
              leadTime[ x,y, test[i+1] ] = length(intervalDays)+1
              multiEvent[x,y,test[i+1] ] = "D"
            }
            if( length(which(t2>0)) == length(t2) ) #如果间隔热浪旗子全为1
            {
              leadTime[ x,y, test[i+1] ] = -(length(intervalDays)+1)
              multiEvent[x,y,test[i+1]] = "H"
            }
            
          }
        }
        
      }
    info<- sprintf("已完成 %d%%", round(x*100/xx))   #查看循环进度
    setTkProgressBar(pb, value = x*100/xx, title = sprintf("进度 (%s)",info),label = info)   
    
  }
  close(pb)  
  
  #逐格点计算正负次数
  DorHarray = array( dim = c(xx,yy) )
  Dmean = array( dim = c(xx,yy) )
  multis = array( dim = c(xx,yy) )
  multisType = array( dim = c(xx,yy) )
  multisD = array( dim = c(xx,yy) )
  multisH = array( dim = c(xx,yy) )
  
  
  pb <- tkProgressBar(title="进度",label="已完成 %", min=0, max=100)
  for(x in 1:xx)
  {
    for(y in 1:yy)
      if(TRUE %in% (Flag[x,y,] >= -999))
      {
        
        CDHindex = which( Flag3[x,y,] == 1 ) #起始位置标记：复合事件开始第一天
        if( length(CDHindex) == 0 )
        {
          DorHarray[x,y] = NA
          Dmean[x,y] = NA
          next
        }
        
        DorHarray[x,y] = length(which(leadTime[x,y,CDHindex]>0)) / length(CDHindex) #计算正值比例
        Dmean[x,y] = mean(leadTime[x,y,CDHindex],na.rm=TRUE) #计算均值
        multis[x,y] = length(which( multiEvent[x,y,CDHindex] == "D" ) ) + 
          length(which( multiEvent[x,y,CDHindex] == "H" ) )
        multisD[x,y] = length(which( multiEvent[x,y,CDHindex] == "D" ) )
        multisH[x,y] = length(which( multiEvent[x,y,CDHindex] == "H" ) )
        
        if( multis[x,y] > 0 )
          multisType[x,y] = length(which( multiEvent[x,y,CDHindex] == "D" ) )/multis[x,y]
        
      }
    
    info<- sprintf("已完成 %d%%", round(x*100/xx))   #查看循环进度
    setTkProgressBar(pb, value = x*100/xx, title = sprintf("进度 (%s)",info),label = info)   
    
  }
  close(pb)  
  
  DorHarray = DorHarray *100 
  DorHarrayr = raster(DorHarray) #正值比例
  extent(DorHarrayr) = c(72, 136, 18, 54)
  # DorHarrayr = raster(DorHarray) #正值比例
  # Dmeanr = raster(Dmean) #差值日均值
  # multisr = raster(multis) #多重事件数
  # multisDr = raster(multisD) #干旱主导多重事件数
  # multisHr = raster(multisH) #干旱主导多重事件数
  # multisTyper = raster(multisType)
  #####
  
    #02 画空间分布图
    #####
    colr = colorRampPalette(brewer.pal(11, "Set3"))(11)[c(4,6,8,11,7)] #颜色
    #MKcol =  cf(colr,whitesite=colr[50],c(-0.5,0.5,0),z0=0) 
    #colr2 = colorRampPalette(colr)(100)
    
    #手动输出为eps
    dev.new(width = 2000,height = 1300)
    par(mar=c(1,1,1,2))
    plot(
      DorHarrayr,
      ylim=c(17,55),xlim=c(72,136),
      xaxs="i",yaxs="i",
      xaxt="n",yaxt="n",
      xlab="",ylab="",
      mgp=c(0,0,0),
      col=colr,
      legend.width=1,
      zlim=c(0,100),breaks = seq(0, 100, 20),
      axis.args=list(at=seq(0, 100, 20),
                     labels=seq(0, 100, 20),
                     cex.axis=1,las=0)
    ) #,
    
    xL = seq(72, 136, 5)
    xlabel = rep("",length(xL)) #设置坐标轴
    yL = seq(18, 54, 5)
    ylabel = rep("",length(yL))
    
    axis(1, xL, xlabel)
    axis(2, yL, ylabel,las=2)
    lines(Chinasp)
    lines(provincesp)
    
    # png(
    #   file = "fig5_spatial.png",
    #   width = 2000,
    #   height = 1300,
    #   res = 72 * 3
    # )
    # 
    # dev.off()
  
    #####

      #03 画过程线01-云南
      #####
      #
      #具体事件
      x=59
      y=57
      testr=r
      testr[x,y]=10000
      plot(testr)
      
      indx = c(19000:19080)  #20130107-20130328
      dimnames(TX)[[3]][indx]
    
      #03-01 事件1
      hstime=indx
      hstime[which(FH[x,y,indx]!=1)]=NA
      dstime=indx
      dstime[which(FD[x,y,indx]!=1)]=NA
      
      dev.new(width=1300,height=800,res=72) #控制输出图片的大小  
      par(mar=c(3,4.5,0.5,4.5))
      #画SPI-左
      plot(indx,SPIarray[x,y,indx], type = "l", xaxt = "n", yaxt = "n",pch=20,
           ylab = "", xlab = "Date",col="deepskyblue", ylim=c(-4,4),
           cex=1.2,lwd=2,mgp=c(1.5,0.8,0),
           xaxs = "i",yaxs = "i")
      
      lines(indx,TSI[x,y,indx], type = "l", pch=20, col="brown1",lwd=2,cex=1.2)
      abline(h=-1,lty=2)
      abline(h=1,lty=2)
      
      points(dstime,SPIarray[x,y,dstime],pch=20,col="deepskyblue",cex=1.5)
      points(hstime,TSI[x,y,hstime],pch=18,col="brown1",cex=1.5)
      
      
      axis(at=seq(-4,4,1)
           ,label=as.character(seq(-4,4,1)),
           side = 2,col.axis="black",cex.axis=0.8,mgp=c(1,0.8,0),
           las=2)
      axis(at= seq(19000,19080,20), 
           label=c("7  Jan 2013","27 Jan 2013",
                   "16 Feb 2013","8  Mar 2013","28 Mar 2013"),
           side=1,mgp = c(2.8, 0.4, 0))
      
      mtext("SPI/-STI",
            side = 2, line = 2.2, cex.lab = 1, 
            las = 0,cex=0.8)
      
      # 画CDH-右
      par(new = TRUE)
      plot(indx,CDHI[x,y,indx], type = "o", xaxt = "n", yaxt = "n",
           ylab = "", xlab = "", pch=17, col="purple2",lwd=2,xaxs = "i",yaxs = "i",
           ylim=c(0.94,1.04))
    
      axis(at=seq(0.980,1,0.01),
           label=as.character(seq(0.980,1,0.01)),
           side = 4, mgp=c(1,0.8,0),
           col.axis="black",cex.axis=0.8,las=2)
    
      mtext("CDHId" ,
            side = 4, line = 2.2, cex.lab = 1,
            las = 0, cex=0.8)
      
      legend(x="top",
             c("SPI","STI","CDHId"),
             lty=c(1,1,1),pch=c(NA,NA,17), lwd=c(2,2,2),
             col=c("deepskyblue","brown1","purple2"),horiz=T,
             bty="n",cex=1.5, x.intersp=c(0.2,0.2,0.2), inset=c(0.01,0) )
      
      #####
      
        #04 画过程线02-东北
        #####
        #具体事件
        x=13
        y=103
        testr=r
        testr[x,y]=10000
        plot(testr)
        
        indx = c(21220:21300)  #20130107-20130328
        dimnames(TX)[[3]][indx]
        
        #03-01 事件1
        hstime=indx
        hstime[which(FH[x,y,indx]!=1)]=NA
        dstime=indx
        dstime[which(FD[x,y,indx]!=1)]=NA
        
        dev.new(width=1300,height=800,res=72) #控制输出图片的大小  
        par(mar=c(3,4.5,0.5,4.5))
        #画SPI-左
        plot(indx,SPIarray[x,y,indx], type = "l", xaxt = "n", yaxt = "n",pch=20,
             ylab = "", xlab = "Date",col="deepskyblue", ylim=c(-3,4),
             cex=1.2,lwd=2,mgp=c(1.5,0.8,0),
             xaxs = "i",yaxs = "i")
        
        lines(indx,TSI[x,y,indx], type = "l", pch=20, col="brown1",lwd=2,cex=1.2)
        abline(h=-1,lty=2)
        abline(h=1,lty=2)
        
        points(dstime,SPIarray[x,y,dstime],pch=20,col="deepskyblue",cex=1.5)
        points(hstime,TSI[x,y,hstime],pch=18,col="brown1",cex=1.5)
        
        
        axis(at=seq(-4,4,1)
             ,label=as.character(seq(-4,4,1)),
             side = 2,col.axis="black",cex.axis=0.8,mgp=c(1,0.8,0),
             las=2)
        axis(at= seq(21220,21300,20), 
             label=c("5  Feb 2019","25 Feb 2019",
                     "17 Mar 2019","6  Apr 2019","26 Apr 2013"),
             side=1,mgp = c(2.8, 0.4, 0))
        
        mtext("SPI/-STI",
              side = 2, line = 2.2, cex.lab = 1, 
              las = 0,cex=0.8)
        
        # 画CDH-右
        par(new = TRUE)
        plot(indx,CDHI[x,y,indx], type = "o", xaxt = "n", yaxt = "n",
             ylab = "", xlab = "", pch=17, col="purple2",lwd=2,xaxs = "i",yaxs = "i",
             ylim=c(0.2,1.5))
        
        axis(at=seq(0.6,1,0.1),
             label=as.character(seq(0.6,1,0.1)),
             side = 4, mgp=c(1,0.8,0),
             col.axis="black",cex.axis=0.8,las=2)
        
        mtext("CDHId" ,
              side = 4, line = 2.2, cex.lab = 1,
              las = 0, cex=0.8)
        
        legend(x="top",
               c("SPI","STI","CDHId"),
               lty=c(1,1,1),pch=c(NA,NA,17), lwd=c(2,2,2),
               col=c("deepskyblue","brown1","purple2"),horiz=T,
               bty="n",cex=1.5, x.intersp=c(0.2,0.2,0.2), inset=c(0.01,0) )
        
        
        #####
        
        #05 画过程线03-云南
        #####
        #具体事件
        x=63
        y=59
        testr=r
        testr[x,y]=10000
        plot(testr)
        
        indx = c(21250:21350)  #20130107-20130328
        dimnames(TX)[[3]][indx]
        
        #03-01 事件1
        hstime=indx
        hstime[which(FH[x,y,indx]!=1)]=NA
        dstime=indx
        dstime[which(FD[x,y,indx]!=1)]=NA
        
        dev.new(width=1300,height=800,res=72) #控制输出图片的大小  
        par(mar=c(3,4.5,0.5,4.5))
        #画SPI-左
        plot(indx,SPIarray[x,y,indx], type = "l", xaxt = "n", yaxt = "n",pch=20,
             ylab = "", xlab = "Date",col="deepskyblue", ylim=c(-4,4),
             cex=1.2,lwd=2,mgp=c(1.5,0.8,0),
             xaxs = "i",yaxs = "i")
        
        lines(indx,TSI[x,y,indx], type = "l", pch=20, col="brown1",lwd=2,cex=1.2)
        abline(h=-1,lty=2)
        abline(h=1,lty=2)
        
        points(dstime,SPIarray[x,y,dstime],pch=20,col="deepskyblue",cex=1.5)
        points(hstime,TSI[x,y,hstime],pch=18,col="brown1",cex=1.5)
        
        
        axis(at=seq(-4,4,1)
             ,label=as.character(seq(-4,4,1)),
             side = 2,col.axis="black",cex.axis=0.8,mgp=c(1,0.8,0),
             las=2)
        axis(at= seq(21250,21350,25), 
             label=c("7  Mar 2019","1  Apr 2019",
                     "26 Apr 2019","21 MAy 2019","15 June 2013"),
             side=1,mgp = c(2.8, 0.4, 0))
        
        mtext("SPI/-STI",
              side = 2, line = 2.2, cex.lab = 1, 
              las = 0,cex=0.8)
        
        # 画CDH-右
        par(new = TRUE)
        plot(indx,CDHI[x,y,indx], type = "o", xaxt = "n", yaxt = "n",
             ylab = "", xlab = "", pch=17, col="purple2",lwd=2,xaxs = "i",yaxs = "i",
             ylim=c(0,1.6))
        
        axis(at=seq(0.6,1,0.1),
             label=as.character(seq(0.6,1,0.1)),
             side = 4, mgp=c(1,0.8,0),
             col.axis="black",cex.axis=0.8,las=2)
        
        mtext("CDHId" ,
              side = 4, line = 2.2, cex.lab = 1,
              las = 0, cex=0.8)
        
        legend(x="top",
               c("SPI","STI","CDHId"),
               lty=c(1,1,1),pch=c(NA,NA,17), lwd=c(2,2,2),
               col=c("deepskyblue","brown1","purple2"),horiz=T,
               bty="n",cex=1.5, x.intersp=c(0.2,0.2,0.2), inset=c(0.01,0) )
        
        
        #####
  
  