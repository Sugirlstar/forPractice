#EF文件出图
#4 Fig4-SPI-STI
#两个工作空间，分别出两幅图

#4-1 打开空间：复合计算-过度识别版本
  #00 Background
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
  #定义函数：计算连续1的长度
  sumfun=function(v,i) #给参:v为某格点逐日spi,起算位置i
  {
    a=0
    while(v[i]==1)
    {
      if( i == length(v) ) #放在i=i+1前面
      {a=a+1; break}
      a=a+1
      i=i+1
    }
    return(a)
  }
  #定义函数：计算连续0的长度
  sumfun2=function(v,i) #给参:v为某格点逐日spi,起算位置i
  {
    a=0
    while(v[i]==0)
    {
      if( i == length(v) )
      {a=a+1; break}
      a=a+1
      i=i+1
    }
    return(a)
  }
  #函数准备工作
  fenqi=function(v,a) #a为起算位置，a为一个向量
  {
    v[is.na(v)==TRUE]=-999
    z=list()
    for(j in 1:length(a) )
    {
      i=j
      inx=a[i]
      b=NULL
      while( v[inx]==1 )
      {
        if(inx == length(v))
        {
          b=c(b,inx)
          break
        }
        b=c(b,inx)
        inx=inx+1
      }
      z[[j]]=b
    }
    return(z)
  }
  
  #设置图片存储位置及名称
  setwd("C:/Users/Yaoying/Desktop/2108论文收尾/211001定稿/EFfig")
  #####
  
    #01 CDHI分级
    #####
    x=30
    y=50
    longtt=71.75+0.5*y #经度 #96.75
    latt=54.25-0.5*x #纬度 #39.25
    cbind(longtt,latt)
  
    #分级
    ind0= which( CDHdi[x,y,] < 0.3 & CDHdi[x,y,] >= 0.2 ) # -0.5
    ind1= which( CDHdi[x,y,] < 0.2 & CDHdi[x,y,] >= 0.1 ) #-0.8
    ind2= which( CDHdi[x,y,] < 0.1 & CDHdi[x,y,] >= 0.05 ) #-1.3
    ind3= which( CDHdi[x,y,] < 0.05 & CDHdi[x,y,] >= 0.02 ) #-1.6
    ind4= which( CDHdi[x,y,] < 0.02 ) #-2
    #具体事件
    indx = c(20331:20369)  #19791113-19800103 CDHE
    indx5 = c(16828:16837) #20131129-20131211 CDHE
    # indx30=ind0[which(Flagtest3[x,y,ind0]==1)]
    
    #####
  
      #02 画分布图
    
      #02画分布图
      #####
      colr=rev(colorRampPalette(brewer.pal(9,"Blues")[1:9])(12))
      #03 Putout
      dev.new(width=1700,height=1700,res=72) #控制输出图片的大小  
      
      xxx=SPIarray[x,y,ind0]
      yyy=-TSI[x,y,ind0]
      zzz=CDHdi[x,y,ind0]
      plot(xxx,yyy,type="p",col=colr[9],
           xlab="SPI",ylab="-STI",pch=16,mgp = c(1.7, 0.7, 0),
           xlim=c(-4,4),ylim=c(-4,4), cex=1.3,
           cex.axis=1.3, cex.lab=1.3,las=1)#
      
      xxx=SPIarray[x,y,ind1]
      yyy=-TSI[x,y,ind1]
      zzz=CDHdi[x,y,ind1]
      points(xxx,yyy,type="p",col=colr[7],xlab="SPI",ylab="-STI",pch=16)
      
      xxx=SPIarray[x,y,ind2]
      yyy=-TSI[x,y,ind2]
      zzz=CDHdi[x,y,ind2]
      points(xxx,yyy,type="p",col=colr[5],xlab="SPI",ylab="-STI",pch=16)
      
      xxx=SPIarray[x,y,ind3]
      yyy=-TSI[x,y,ind3]
      zzz=CDHdi[x,y,ind3]
      points(xxx,yyy,type="p",col=colr[3],xlab="SPI",ylab="-STI",pch=16)
      
      xxx=SPIarray[x,y,ind4]
      yyy=-TSI[x,y,ind4]
      zzz=CDHdi[x,y,ind4]
      points(xxx,yyy,type="p",col=colr[1],xlab="SPI",ylab="-STI",pch=16)
      
      xxx=SPIarray[x,y,indx5]
      yyy=-TSI[x,y,indx5]
      zzz=CDHdi[x,y,indx5]
      points(xxx,yyy,type="p",col="darkorange",xlab="SPI",ylab="-STI",pch=3,cex=1.5)
      
      xxx=SPIarray[x,y,indx]
      yyy=-TSI[x,y,indx]
      zzz=CDHdi[x,y,indx]
      points(xxx,yyy,type="p",col="darkorange",xlab="SPI",ylab="-STI",pch=4, cex=1.5)
      
      #
      text(2.1,3.7,"  previous index \n ")
      
      points(2.1,3.2, pch=16,cex=1.8,col=colr[9])
      text(2.7,3.2," <0.3")
      
      points(2.1,2.9, pch=16,cex=1.8,col=colr[7])
      text(2.7,2.9, "<0.2" )
      
      points(2.1,2.6, pch=16,cex=1.8,col=colr[5])
      text(2.7,2.6, "<0.1" )
      
      points(2.1,2.3, pch=16,cex=1.8,col=colr[3])
      text(2.7,2.3, "<0.05" )
      
      points(2.1,2.0, pch=16,cex=1.8,col=colr[1])
      text(2.7,2.0, "<0.02" )
      
      abline(h=-1,lty=2)
      abline(v=-1,lty=2)
      
      
      points(1.8,-3.6, pch=4,cex=1.8,col="darkorange")
      text(3.1,-3.6,"20160830-20161007",cex=1.1)
      
      points(1.8,-4, pch=3,cex=1.8,col="darkorange")
      text(3.1,-4,"20070127-20070205",cex=1.1)
      
      #####

        #03 画过程线
        #####
        #具体事件
        indx = c(20331:20369)  #20160830-20161007 CDHE
        indx5 = c(16828:16837) #20070127-20070205 CDHE
        
        #03-01 事件1
        dev.new(width=1300,height=800,res=72) #控制输出图片的大小  
        par(mar=c(3,4.8,0.5,3))
        #画SPI-左
        plot(indx,SPIarray[x,y,indx], type = "o", xaxt = "n", yaxt = "n",pch=20,
             ylab = "", xlab = "Date",col="darkblue", ylim=c(-3.5,3.5),
             cex=0.8,lwd=2,mgp=c(1.5,0.8,0),
             xaxs = "i",yaxs = "i")
        
        lines(indx,-TSI[x,y,indx], type = "o", pch=20, col="brown",lwd=2)
        abline(h=-1,lty=2)
        
        axis(at=seq(-3.5,3.5,1)
             ,label=as.character(seq(-3.5,3.5,1)),
             side = 2,col.axis="black",cex.axis=0.8,mgp=c(1,0.8,0),
             las=2)
        axis(at=c(seq(20331,20369,10),20369), 
             label=c("27 Aug 2016","9 Sep 2016",
                     "19 Sep 2016","29 Sep 2016","7 Oct 2016"),
             side=1,mgp = c(2.8, 0.4, 0))
        
        mtext("SPI/-STI",
              side = 2, line = 2.2, cex.lab = 1, 
              las = 0, col="darkblue",cex=0.8)
        
        legend(x="bottomright",
               c("SPI","-STI"),
               lty=c(1,1),pch=c(20,20), lwd=c(2,2),col=c("darkblue","brown"),
               text.col=c("darkblue","brown"),
               bty="n",cex=1.5, x.intersp=c(0.2,0.2), inset=c(0.01,0) )
    
        # 略_画CDH-右
        # par(new = TRUE)
        # plot(indx,CDHdi[x,y,indx], type = "l", xaxt = "n", yaxt = "n", 
        #      ylab = "", xlab = "", pch=20, col="darkorange",lwd=2,xaxs = "i",yaxs = "i",
        #      ylim=c(0.00,0.16))
        # 
        # axis(at=seq(0.00,0.16,0.02), 
        #      label=as.character(seq(0.00,0.16,0.02)),
        #      side = 4, mgp=c(1,0.8,0),
        #      col.axis="black",cex.axis=0.8,las=2)
        # 
        # mtext("CDH index" ,
        #       side = 4, line = 2.2, cex.lab = 1, 
        #       las = 0, cex=0.8)
        
        #03-01 事件2
        dev.new(width=1300,height=800,res=72) #控制输出图片的大小  
        par(mar=c(3,4.8,0.5,3))
        #画SPI-左
        plot(indx5,SPIarray[x,y,indx5], type = "o", xaxt = "n", yaxt = "n",pch=20,
             ylab = "", xlab = "Date",col="darkblue", ylim=c(-3.5,3.5),
             cex=0.8,lwd=2,mgp=c(1.5,0.8,0),
             xaxs = "i",yaxs = "i")
        
        lines(indx5,-TSI[x,y,indx5], type = "o", pch=20, col="brown",lwd=2)
        abline(h=-1,lty=2)
        
        axis(at=seq(-3.5,3.5,1)
             ,label=as.character(seq(-3.5,3.5,1)),
             side = 2,col.axis="black",cex.axis=0.8,mgp=c(1,0.8,0),
             las=2)
        
        axis(at=seq(16828,16837,3), 
             label=c("27 Jan 2007","30 Jan 2007",
                     "2 Feb 2007","5 Feb 2007"),
             side=1,mgp = c(2.8, 0.4, 0))
        
        mtext("SPI/-STI",
              side = 2, line = 2.2, cex.lab = 1, 
              las = 0, col="darkblue",cex=0.8)
        
        legend(x="bottomright",
               c("SPI","-STI"),
               lty=c(1,1),pch=c(20,20), lwd=c(2,2),col=c("darkblue","brown"),
               text.col=c("darkblue","brown"),
               bty="n",cex=1.5, x.intersp=c(0.2,0.2), inset=c(0.01,0) )
        #####
        
  
#4-2 打开空间：复合计算_final  
  
  #00 Background
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
  #定义函数：计算连续1的长度
  sumfun=function(v,i) #给参:v为某格点逐日spi,起算位置i
  {
    a=0
    while(v[i]==1)
    {
      if( i == length(v) ) #放在i=i+1前面
      {a=a+1; break}
      a=a+1
      i=i+1
    }
    return(a)
  }
  #定义函数：计算连续0的长度
  sumfun2=function(v,i) #给参:v为某格点逐日spi,起算位置i
  {
    a=0
    while(v[i]==0)
    {
      if( i == length(v) )
      {a=a+1; break}
      a=a+1
      i=i+1
    }
    return(a)
  }
  #函数准备工作
  fenqi=function(v,a) #a为起算位置，a为一个向量
  {
    v[is.na(v)==TRUE]=-999
    z=list()
    for(j in 1:length(a) )
    {
      i=j
      inx=a[i]
      b=NULL
      while( v[inx]==1 )
      {
        if(inx == length(v))
        {
          b=c(b,inx)
          break
        }
        b=c(b,inx)
        inx=inx+1
      }
      z[[j]]=b
    }
    return(z)
  }
  
  #设置图片存储位置及名称
  setwd("C:/Users/Yaoying/Desktop/2108论文收尾/211001定稿/EFfig")
  #####
  
    #01 CDHI分级
    #####
    x = 30
    y = 50
    longtt=71.75+0.5*y #经度
    latt=54.25-0.5*x #纬度
    
    #分级:0.3 0.2 0.1 0.05 0.02
    ind0=which( CDHdi[x,y,] > 0 & CDHdi[x,y,]<=0.3 )
    ind1=which(CDHdi[x,y,]> 0.3 & CDHdi[x,y,]<=0.6 )
    ind2=which(CDHdi[x,y,]> 0.6 & CDHdi[x,y,]<=0.9 )
    ind3=which(CDHdi[x,y,]> 0.9 & CDHdi[x,y,]<=0.98 )
    ind4=which(CDHdi[x,y,]> 0.98)
    indxx = which( CDHdi[x,y,] > 0)
    indx=indxx[which(Flag2[x,y,indxx]==1)] #这当中被识别为复合事件的点
    #####
    
      #02-1 画分布图-01版
      #####
      colr=rev(colorRampPalette(brewer.pal(9,"Blues")[1:9])(12))
      #03 Putout
      dev.new(width=1700,height=1700,res=72) #控制输出图片的大小  
      
      
      
      xxx=SPIarray[x,y,ind0]
      yyy=-TSI[x,y,ind0]
      zzz=CDHI[x,y,ind0]
      plot(xxx,yyy,type="p",col=colr[9],
           xlab="SPI",ylab="-STI",pch=16,mgp = c(1.7, 0.7, 0),
           xlim=c(-4,4),ylim=c(-4,4), cex=1.3,
           cex.axis=1.3, cex.lab=1.3,las=1)#
      
      xxx=SPIarray[x,y,ind1]
      yyy=-TSI[x,y,ind1]
      zzz=CDHI[x,y,ind1]
      points(xxx,yyy,type="p",col=colr[7],xlab="SPI",ylab="-STI",pch=16)
      
      xxx=SPIarray[x,y,ind2]
      yyy=-TSI[x,y,ind2]
      zzz=CDHI[x,y,ind2]
      points(xxx,yyy,type="p",col=colr[5],xlab="SPI",ylab="-STI",pch=16)
      
      xxx=SPIarray[x,y,ind3]
      yyy=-TSI[x,y,ind3]
      zzz=CDHI[x,y,ind3]
      points(xxx,yyy,type="p",col=colr[3],xlab="SPI",ylab="-STI",pch=16)
      
      xxx=SPIarray[x,y,ind4]
      yyy=-TSI[x,y,ind4]
      zzz=CDHI[x,y,ind4]
      points(xxx,yyy,type="p",col=colr[1],xlab="SPI",ylab="-STI",pch=16)
      
      xxx=SPIarray[x,y,indx]
      yyy=-TSI[x,y,indx]
      zzz=CDHI[x,y,indx]
      points(xxx,yyy,type="p",col="orangered",xlab="SPI",ylab="-STI",pch=4,cex=1.5)
      
      
      text(2.1,3.7,"  CDHId \n ")
      points(2.1,3.2, pch=16,cex=1.8,col=colr[9])
      text(2.7,3.2," >0")
      points(2.1,2.9, pch=16,cex=1.8,col=colr[7])
      text(2.7,2.9, ">0.3" )
      points(2.1,2.6, pch=16,cex=1.8,col=colr[5])
      text(2.7,2.6, ">0.6" )
      points(2.1,2.3, pch=16,cex=1.8,col=colr[3])
      text(2.7,2.3, ">0.9" )
      points(2.1,2.0, pch=16,cex=1.8,col=colr[1])
      text(2.7,2.0, ">0.98" )
      
      abline(h=-1,lty=2)
      abline(v=-1,lty=2)
      
      points(1.8,-3.6, pch=4,cex=1.8,col="orangered")
      text(3.1,-3.6,"19950526-19950604",cex=1.1)
      
      #####
    
        #02-2 画分布图-02版
        #####
        colr=rev(colorRampPalette(brewer.pal(9,"Blues")[1:9])(12))
        #03 Putout
        dev.new(width=1700,height=1700,res=72) #控制输出图片的大小  
        
        xxx=SPIarray[x,y,ind0]
        yyy=-TSI[x,y,ind0]
        zzz=CDHI[x,y,ind0]
        plot(xxx,yyy,type="p",col=colr[9],
             xlab="SPI",ylab="-STI",pch=16,mgp = c(3, 0.7, 0),
             xlim=c(-4,-1),ylim=c(-4,-1), cex=1.3,
             cex.axis=1.3, cex.lab=1.3,las=1)#
        
        xxx=SPIarray[x,y,ind1]
        yyy=-TSI[x,y,ind1]
        zzz=CDHI[x,y,ind1]
        points(xxx,yyy,type="p",col=colr[7],xlab="SPI",ylab="-STI",pch=16)
        
        xxx=SPIarray[x,y,ind2]
        yyy=-TSI[x,y,ind2]
        zzz=CDHI[x,y,ind2]
        points(xxx,yyy,type="p",col=colr[5],xlab="SPI",ylab="-STI",pch=16)
        
        xxx=SPIarray[x,y,ind3]
        yyy=-TSI[x,y,ind3]
        zzz=CDHI[x,y,ind3]
        points(xxx,yyy,type="p",col=colr[3],xlab="SPI",ylab="-STI",pch=16)
        
        xxx=SPIarray[x,y,ind4]
        yyy=-TSI[x,y,ind4]
        zzz=CDHI[x,y,ind4]
        points(xxx,yyy,type="p",col=colr[1],xlab="SPI",ylab="-STI",pch=16)
        
        xxx=SPIarray[x,y,indx]
        yyy=-TSI[x,y,indx]
        zzz=CDHI[x,y,indx]
        points(xxx,yyy,type="p",col="orangered",xlab="SPI",ylab="-STI",pch=4,cex=1.5)
        
        text(-3.5,-2.95,"CDHId\n ")
        points(-3.7,-3.2, pch=16,cex=1.8,col=colr[9])
        text(-3.5,-3.2,"> 0")
        points(-3.7,-3.4, pch=16,cex=1.8,col=colr[7])
        text(-3.5,-3.4, "> 0.3" )
        points(-3.7,-3.6, pch=16,cex=1.8,col=colr[5])
        text(-3.5,-3.6, "> 0.6" )
        points(-3.7,-3.8, pch=16,cex=1.8,col=colr[3])
        text(-3.5,-3.8, "> 0.9", )
        points(-3.7,-4.0, pch=16,cex=1.6,col=colr[1])
        text(-3.5,-4.0, "> 0.98" )
    
        abline(h=-1,lty=2)
        abline(v=-1,lty=2)
        
        points(-4,-1.2, pch=4,cex=1.8,col="orangered")
        text(-3.5,-1.2,"19950526-19950604",cex=1.1)
        
        #####
        
            #03 画过程线
          #####
          #具体事件
          #03-01 事件1
          dev.new(width=1300,height=800,res=72) #控制输出图片的大小  
          par(mar=c(3,4.8,0.5,3))
          #画SPI-左
          plot(indx,SPIarray[x,y,indx], type = "o", xaxt = "n", yaxt = "n",pch=20,
               ylab = "", xlab = "Date",col="darkblue", ylim=c(-3.5,-0.5),
               cex=0.8,lwd=2,mgp=c(1.5,0.8,0),
               xaxs = "i",yaxs = "i")
          
          lines(indx,-TSI[x,y,indx], type = "o", pch=20, col="brown",lwd=2)
          abline(h=-1,lty=2)
          
          axis(at=seq(-3.5,-0.5,1)
               ,label=as.character(seq(-3.5,-0.5,1)),
               side = 2,col.axis="black",cex.axis=0.8,mgp=c(1,0.8,0),
               las=2)
          axis(at=seq(12564,12573,3), 
               label=c("26 May 1995","29 May 1995",
                       "1 June 1995","4 June 1995"),
               side=1,mgp = c(2.8, 0.4, 0))
          
          mtext("SPI/-STI",
                side = 2, line = 2.2, cex.lab = 1, 
                las = 0, col="black",cex=0.8)
          
          legend(x="bottomright",
                 c("SPI","-STI"),
                 lty=c(1,1),pch=c(20,20), lwd=c(2,2),col=c("darkblue","brown"),
                 text.col=c("darkblue","brown"),
                 bty="n",cex=1.5, x.intersp=c(0.2,0.2), inset=c(0.01,0) )
          
          #####
    
    
  
  