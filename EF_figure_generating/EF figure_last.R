#figure last

#01 Background
#####
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

tendency=function(v) #默认:由1961开始
{
  long=length(v)
  fit=lm(v~c(1961:(1961+long-1)))
  ten=fit$coefficients[2]
  intercept=fit$coefficients[1]
  pvalue=summary(fit)$coefficients[,4][2]
  tp=list(ten,pvalue,intercept)
  return(tp)
}
#01-02FLMIP（一个突变点）函数
bp1=function(s,e,STGh) 
{
  muti=function(v)  #用于检验相邻两斜率是否符号相反
  {
    x=length(v)-1
    z=0
    for( i in 1:x )
    { p=v[i]*v[i+1]
    if(p>0)
      break
    z=z+1
    }
    return(x==z)
  }
  
  m=e-s+1  #总年数
  Y=matrix(STGh,m,1)   #逐年STG值  m行矩阵
  T=s:e   #逐年年份 
  
  SSR=10000000
  RST=0
  
  bp=(s+9):(e-9)   ##一个点的选取规则不太一样
  ##2020.10.16修改:前后至少十年
  b=bp-s+1
  for( j in 1:length(bp) )  #b[j]为转折年份
  {
    A0=matrix(0,m,2)  #1个转折点对应的矩阵A
    A0[1: b[j],1]=s:bp[j]
    A0[ (b[j]+1):m ,1 ]=bp[j]
    A0[ (b[j]+1):m, 2 ]=1:(e-bp[j])
    c=matrix(rep(1,m),m,1)  #最后一列c1值    
    A=cbind(A0,c)
    S=solve(t(A)%*%A)%*%t(A)%*%Y
    
    if(muti(S[1:2]))   #满足条件：相邻两段斜率符号相反
    {
      rst=matrix(0,4,2)  #放最终结果 
      rownames(rst)=c("year","a","c","p of a") #拟合结果：y=a*x+c
      rst[1,1]=bp[j]
      rst[2,]=S[1: 2]
      rst[3,1]=S[3]
      rst[3,2]=rst[3,1]+(S[1]-S[2])*bp[j]
      
      ##显著性水平筛选
      
      rdf1=rst[1,1]-s+1-2 #自由度(总年数-2)
      rdf2=e-rst[1,1]+1-2
      
      ssr1=0 #残差平方和
      for(u in 1:b[j])
        ssr1=ssr1+(Y[u]-(rst[2,1]*T[u]+rst[3,1]))^2
      
      ssr2=0 #残差平方和
      for( u in (b[j]+1):m )
        ssr2=ssr2+(Y[u]-(rst[2,2]*T[u]+rst[3,2]))^2  #两段之和
      
      vv1=sqrt(ssr1/rdf1) #Residual standard error
      vv2=sqrt(ssr2/rdf2) #Residual standard error
      
      stderr1=vv1/sqrt(sum( (c(s:rst[1,1])-mean(c(s:rst[1,1])))^2)) #斜率的标准误差
      stderr2=vv2/sqrt(sum( ( c((rst[1,1]+1):e)-mean(c((rst[1,1]+1):e)) )^2)) #斜率的标准误差
      
      tval1=rst[2,1]/stderr1
      tval2=rst[2,2]/stderr2
      
      pr1=2*pt(abs(tval1),rdf1,lower.tail=FALSE)
      pr2=2*pt(abs(tval2),rdf2,lower.tail=FALSE)
      
      rst[4,1]=pr1
      rst[4,2]=pr2
      
      #计算总残差平方和
      ssr=ssr1+ssr2
      
      #显著性检验：若不通过则重置ssr为极大，不可能留下
      
      if( pr1>0.05 & pr2>0.05 )
        ssr=1000000000000  
      
      #冒泡比较 留下更小的
      if(ssr<SSR)
      {SSR=ssr
      RST=rst}
    }
  }
  return(RST)
}

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
#####


#02 read data
#####
s=1961
e=2020
setwd("C:/Users/Yaoying/Desktop/2108论文收尾/211001定稿")
#
ctype = "GLOTI"
df = read.table("C:/Users/Yaoying/Desktop/2108论文收尾/211001定稿/GLOTI_1961_2020.csv", header=F, sep=",")
GLOTI6120 = df[,2]
#02
ctype = "PDO"
df = read.table("C:/Users/Yaoying/Desktop/2108论文收尾/211001定稿/PDO61_20.csv", header=F, sep=",")
PDO6120 = df[,2]

#####

#03 出图
#####
objx=PDO6120 #objx=DRFREh
ydur = max(objx,na.rm=T) - min(objx,na.rm=T)
getylim1 = c(min(objx,na.rm=T), max(objx,na.rm=T)+ydur*0.2)
getylim1;ydur/4
#得到y刻度线

#03 Putout
dev.new(width=1300,height=3000,res=72) #控制输出图片的大小  
par( mfrow = c(2,1),mar=c(4,4.8,0.5,1.2) )   
figure2out(PDO6120,"PDO",seq(-2.0,1.5,0.5))        
figure2out(GLOTI6120,"Temperature anomaly",seq(-0.3,1.2,0.3))        





