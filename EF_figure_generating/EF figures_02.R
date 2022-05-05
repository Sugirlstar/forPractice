#EF文件出图
#2 Fig2-Temporal_variation_4r3c
#####
# seperate in to 3 workspace

#2-1 DR
# Open the workspace:  .../干旱计算全流程.RData
#####
  #00 Background setting
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
  #画图坐标
  xL = seq(72, 136, 5)
  xlabel = paste(xL, "°", sep = "") #设置坐标轴
  yL = seq(18, 54, 5)
  ylabel = paste(yL, "°", sep = "")
  
    #01 functions
    figure2out = function(objx,ylabname,aty) 
    {
      atylabel=as.character(aty)
      ylimget = c(min(aty,na.rm=T),max(aty,na.rm=T))
      
      objx=-DRSTGh
      #计算突变点
      rst=bp1(s,e,objx)
      #02-01-01-a
      if(rst[1]!=0) #如果有突变点
      {
        FREt=c(s,rst[1,(1:(ncol(rst)-1))],e)
        FREy=FREt
        FREy[1]=FREt[1]*rst[2,1]+rst[3,1]
        for(i in 2: length(FREt) )
          FREy[i]=FREt[i]*rst[2,i-1]+rst[3,i-1]
        FREp1=rst[4,1] ###
        FREp2=rst[4,2] ###
      }
      #02-01-01-b
      if(rst[1]==0) #如果无突变点——进行线性拟合、MK检验
      {
        tend=tendency(objx) #直线拟合
        tendx=c(s:e)
        tendy=tendx*tend[[1]]+tend[[3]] #计算y坐标
        #MK趋势检验
        mktau=mk.test(objx)$estimates[3] #tau
        mkpvalue=mk.test(objx)$p.value #p值
      }
      
      
      atx = seq(1960,2020,10)
      atxlabel = as.character(atx)
      plot(c(s:e),objx,type="o",
           xlim = c(1960,2020),
           ylim = ylimget,
           ylab = ylabname,
           pch=20, xlab="Year",
           xaxs = "i",yaxs = "i", xaxt = "n", yaxt = "n",
           lwd=0.8,
           mgp=c(2.4,0.4,0),cex.lab=1.2) 
      
      axis(1,at=atx,atxlabel)
      
      axis(2,at=aty,atylabel,las = 2)
      
      
      #Linear trends
      if(rst[1]!=0)
      {
        lines(FREt, FREy, lwd=2,col="red")
        points(FREt[2],FREy[2],pch=17,cex=3)
      }
      #M-K
      if(rst[1]==0)
        lines(tendx,tendy,lwd=2,col="blue")
      
    }
    
      #02 Test for ylim (not run)
      objx=DRper*100 #objx=DRFREh
      ydur = max(objx,na.rm=T) - min(objx,na.rm=T)
      getylim1 = c(min(objx,na.rm=T), max(objx,na.rm=T)+ydur*0.2)
      getylim1;ydur/4
      #得到y刻度线
      
          #03 Putout
          dev.new(width=1300,height=3000,res=72) #控制输出图片的大小  
          par( mfrow = c(4,1),mar=c(4,4.8,0.5,1.2) )       
          figure2out(DRFREh,"Frequency (times)",seq(0.6,2.2,0.4))        
          figure2out(DRDURh,"Duration (days)",seq(20,100,20))        
          figure2out(-DRSTGh,"Severity",seq(0.6,1.6,0.2))        
          figure2out(DRper*100,"Coverage (%)",seq(40,100,10)) 
          #手动另存为eps
          
              #04 输出数值结果
              rst1=bp1(s,e,DRFREh)
              rst2=bp1(s,e,DRDURh)
              rst3=bp1(s,e,DRSTGh)
              rst4=bp1(s,e,DRper)
              csvout = cbind(rst1,rst2,rst3,rst4)
              write.table (csvout,sep=",",file ="dr_char.csv") #文件读写

#####
          
          
#2-2 HW
# Open the workspace:  .../热浪计算全流程.RData
#####
#00 Background setting
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
#画图坐标
xL = seq(72, 136, 5)
xlabel = paste(xL, "°", sep = "") #设置坐标轴
yL = seq(18, 54, 5)
ylabel = paste(yL, "°", sep = "")

      #01 functions
      figure2out = function(objx,ylabname,aty) 
      {
        atylabel=as.character(aty)
        ylimget = c(min(aty,na.rm=T),max(aty,na.rm=T))
        
        #计算突变点
        rst=bp1(s,e,objx)
        #02-01-01-a
        if(rst[1]!=0) #如果有突变点
        {
          FREt=c(s,rst[1,(1:(ncol(rst)-1))],e)
          FREy=FREt
          FREy[1]=FREt[1]*rst[2,1]+rst[3,1]
          for(i in 2: length(FREt) )
            FREy[i]=FREt[i]*rst[2,i-1]+rst[3,i-1]
          FREp1=rst[4,1] ###
          FREp2=rst[4,2] ###
        }
        #02-01-01-b
        if(rst[1]==0) #如果无突变点——进行线性拟合、MK检验
        {
          tend=tendency(objx) #直线拟合
          tendx=c(s:e)
          tendy=tendx*tend[[1]]+tend[[3]] #计算y坐标
          #MK趋势检验
          mktau=mk.test(objx)$estimates[3] #tau
          mkpvalue=mk.test(objx)$p.value #p值
        }
        
        
        atx = seq(1960,2020,10)
        atxlabel = as.character(atx)
        plot(c(s:e),objx,type="o",
             xlim = c(1960,2020),
             ylim = ylimget,
             ylab = ylabname,
             pch=20, xlab="Year",
             xaxs = "i",yaxs = "i", xaxt = "n", yaxt = "n",
             lwd=0.8,
             mgp=c(2.4,0.4,0),cex.lab=1.2) 
        
        axis(1,at=atx,atxlabel)
        
        axis(2,at=aty,atylabel,las = 2)
        
        
        #Linear trends
        if(rst[1]!=0)
        {
          lines(FREt, FREy, lwd=2,col="red")
          points(FREt[2],FREy[2],pch=17,cex=3)
        }
        #M-K
        if(rst[1]==0)
          lines(tendx,tendy,lwd=2,col="blue")
        
      }


          #02 Test for ylim (not run)
          objx=DRper*100 #objx=DRFREh
          ydur = max(objx,na.rm=T) - min(objx,na.rm=T)
          getylim1 = c(min(objx,na.rm=T), max(objx,na.rm=T)+ydur*0.2)
          getylim1;ydur/4
          #得到y刻度线
          
                #03 Putout
                dev.new(width=1300,height=3000,res=72) #控制输出图片的大小  
                par( mfrow = c(4,1),mar=c(4,4.8,0.5,1.2) )       
                figure2out(DRFREh,"Frequency (times)",seq(1,6,1))        
                figure2out(DRDURh,"Duration (days)",seq(6,41,7))        
                figure2out(DRSTGh,"Severity",seq(0.8,2.0,0.2))        
                figure2out(DRper*100,"Coverage (%)",seq(50,110,10))        
                #手动另存为eps
                
                    #04 输出数值结果
                    rst1=bp1(s,e,DRFREh)
                    rst2=bp1(s,e,DRDURh)
                    rst3=bp1(s,e,DRSTGh)
                    rst4=bp1(s,e,DRper)
                    csvout = cbind(rst1,rst2,rst3,rst4)
                    write.table (csvout,sep=",",file ="hw_char.csv") #文件读写
                    

#####
                
                
#2-3 CDH
# Open the workspace:  .../复合计算_final.RData
#####
#00 Background setting
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
#画图坐标
xL = seq(72, 136, 5)
xlabel = paste(xL, "°", sep = "") #设置坐标轴
yL = seq(18, 54, 5)
ylabel = paste(yL, "°", sep = "")

    #01 functions
    figure2out = function(objx,ylabname,aty) 
    {
      atylabel=as.character(aty)
      ylimget = c(min(aty,na.rm=T),max(aty,na.rm=T))
      
      #计算突变点
      rst=bp1(s,e,objx)
      #02-01-01-a
      if(rst[1]!=0) #如果有突变点
      {
        FREt=c(s,rst[1,(1:(ncol(rst)-1))],e)
        FREy=FREt
        FREy[1]=FREt[1]*rst[2,1]+rst[3,1]
        for(i in 2: length(FREt) )
          FREy[i]=FREt[i]*rst[2,i-1]+rst[3,i-1]
        FREp1=rst[4,1] ###
        FREp2=rst[4,2] ###
      }
      #02-01-01-b
      if(rst[1]==0) #如果无突变点——进行线性拟合、MK检验
      {
        tend=tendency(objx) #直线拟合
        tendx=c(s:e)
        tendy=tendx*tend[[1]]+tend[[3]] #计算y坐标
        #MK趋势检验
        mktau=mk.test(objx)$estimates[3] #tau
        mkpvalue=mk.test(objx)$p.value #p值
      }
      
      
      atx = seq(1960,2020,10)
      atxlabel = as.character(atx)
      plot(c(s:e),objx,type="o",
           xlim = c(1960,2020),
           ylim = ylimget,
           ylab = ylabname,
           pch=20, xlab="Year",
           xaxs = "i",yaxs = "i", xaxt = "n", yaxt = "n",
           lwd=0.8,
           mgp=c(2.4,0.4,0),cex.lab=1.2) 
      
      axis(1,at=atx,atxlabel)
      
      axis(2,at=aty,atylabel,las = 2)
      
      
      #Linear trends
      if(rst[1]!=0)
      {
        lines(FREt, FREy, lwd=2,col="red")
        points(FREt[2],FREy[2],pch=17,cex=3)
      }
      #M-K
      if(rst[1]==0)
        lines(tendx,tendy,lwd=2,col="blue")
      
    }

        #02 Test for ylim (not run)
        objx=CDHFREh #objx=DRFREh
        ydur = max(objx,na.rm=T) - min(objx,na.rm=T)
        getylim1 = c(min(objx,na.rm=T), max(objx,na.rm=T)+ydur*0.2)
        getylim1;ydur/4
        #得到y刻度线

            #03 Putout
            dev.new(width=1300,height=3000,res=72) #控制输出图片的大小  
            par( mfrow = c(4,1),mar=c(4,4.8,0.5,1.2) )       
            figure2out(CDHFREh,"Frequency (times)",seq(0,0.30,0.05))        
            figure2out(CDHDURh,"Duration (days)",seq(0.0,3.50,0.7))        
            figure2out(CDHMAGh,"Magnitude",seq(0,3,0.6))        
            figure2out(CDHper*100,"Coverage (%)",seq(0,20,4))        
            #手动另存为eps
            
            #04 输出数值结果
            rst1=bp1(s,e,CDHFREh)
            rst2=bp1(s,e,CDHDURh)
            rst3=bp1(s,e,CDHMAGh)
            rst4=bp1(s,e,CDHper)
            csvout = cbind(rst1,rst2,rst3,rst4)
            write.table (csvout,sep=",",file ="CDH_char.csv") #文件读写

#####
                
                

            
                
                
                
