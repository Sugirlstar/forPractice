#SPI与STI相关性计算+HDId指数计算+指标计算+出草图

#00 Background setting
#####
#加载必要空间
load("C:/Users/Yaoying/Desktop/2108论文收尾/211001定稿/SPI_TX.RData") #SPIarray
load("C:/Users/Yaoying/Desktop/2108论文收尾/211001定稿/STI.RData") #SPIarray
load("C:/Users/Yaoying/Desktop/2108论文收尾/211001定稿/r层.RData") #r层 吉祥物及水印作用:)
#设置图片存储位置及名称
setwd("C:/Users/Yaoying/Desktop/2108论文收尾/211001定稿/compound_draft")
ctype = "compound_-1-1-10days"
#加载必要包
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
#读取shp文件
Chinasp = readShapePoly("C:/Users/Yaoying/Desktop/hyj/china/China") #中国边界
provincesp = readShapePoly("C:/Users/Yaoying/Desktop/hyj/china/PROVINCE_region") #省份边界
#基本参数
spithr1 = -1 #干旱第一阈值
hwthr = -1 #-STI阈值
da=10 #最低连续时间
xx=dim(SPIarray)[1]
yy=dim(SPIarray)[2]
days=dim(SPIarray)[3]
st=90 #起算日
s=1961
e=2020
LL=e-s+1
# da2=5 #最低间隔时间(复合事件无合并过程)
#画图坐标
xL = seq(72, 136, 5)
xlabel = paste(xL, "°", sep = "") #设置坐标轴
yL = seq(18, 54, 5)
ylabel = paste(yL, "°", sep = "")
#####

#01 Functions preparing
#####
#01-01变化率计算
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
#01-03分类符号函数
classtype=function(v)
{
  if(is.na(v)==TRUE)
    classtype=NA else
      if(v<=0.05)   #0.05 点
        classtype=4 else 
          classtype=NA
        return(classtype)
}
#01-04色谱定位
cf=function(colr,whitesite="white",z,z0=0) 
{
  #输入色带数值，
  #中间颜色名（或色号#）,
  #z值(raster或矩阵层),
  #z0值(必须包含在z值中),
  #返回可用的颜色序列
  
  zz=as.matrix(z)
  z1=min(zz,na.rm=TRUE) #z最小
  z2=max(zz,na.rm=TRUE) #z最大
  zL1=(z0-z1)/(z2-z1) #第一段z值比例
  zL2=(z2-z0)/(z2-z1) #第二段z值比例
  cL=length(colr)
  if(whitesite=="mid")
    c0=round(cL/2)  else
      c0=which(colr==whitesite) #白色位置
  cL1=(c0-1)/cL #第一段色带比例
  cL2=(cL-c0+1)/cL #第二段色带比例
  
  if(z0<z1) #z0比z1小
  {
    x=round((z1-z0)/(z2-z0)*(cL-c0)+c0)
    colr_result=colr[x:cL]
  }else
    if(z0>z2)
    {
      x=round((z2-z1)/(z0-z1)*c0)
      colr_result=colr[1:x]
    } else
      if(zL1==0)
        colr_result=colr[c0:cL] else
          if(zL2==0)
            colr_result=colr[c0:cL]=colr[1:c0] else
              if(zL1>cL1)
              {
                x=round(c0/zL1) #修改后的长度
                colr_result=colr[1:x]
              }else
                if(zL1<cL1)
                {
                  x=round((zL2*cL+c0-cL)/zL2) #修改后的长度
                  colr_result=colr[x:cL]
                }else
                  colr_result=colr
  
  return(colr_result)        
}
#01-05逐格点M-K检验
MK.raster <- function(xraster, type="year", month0=1) #默认为year
{
  library(Kendall)
  library(raster)
  x<-as.array(xraster) #将raster转为矩阵
  year0=1961 #开始年份
  D<-dim(x)
  MK.xraster<-array(data=NA,dim=c(D[1],D[2],3)) #MK.xraster的dim:c(xx,yy,3)
  #采用MK检验计算【逐年】趋势
  if (type == "year"){
    if( all ( subset(x,is.na(x)==FALSE) >= 0) ) #如果所有值均>=0
    {
      for (i in 1:D[1])
        for (j in 1:D[2])
          if (TRUE %in% (x[i, j, ] >= -9999))
          {
            if( length( which(x[i,j,]>0) ) > 2 ) #要求：至少有3个值
            {
              xts<-ts(x[i,j,],start=year0,frequency=1)
              z<-MannKendall(xts)
              MK.xraster[i,j,1:2]<- c(z$tau,z$sl)
            }else
              MK.xraster[i,j,1:2]=NA #否则为NA
          }
    } else
    {
      for (i in 1:D[1])
        for (j in 1:D[2])
          if (TRUE %in% (x[i, j, ] >= -9999))
          {
            if( length( which(x[i,j,]> -9999 ) ) > 2 ) #要求：至少有3个值
            {
              xts<-ts(x[i,j,],start=year0,frequency=1)
              z<-MannKendall(xts)
              MK.xraster[i,j,1:2]<- c(z$tau,z$sl)
            }else
              MK.xraster[i,j,1:2]=NA #否则为NA
          }
    }
  }
  return(MK.xraster)
}

#####

#02-01 Corresponding
#####
#此步保存-运行约30mins

# rix=c(1:dim(SPIarray)[3])
# CDHcor=pcor=array(dim=c(xx,yy))
# pb = tkProgressBar(title="进度",label="已完成 %", min=0, max=100, initial = 0, width = 300)
# for (x in 1:xx)
# {
#   for (y in 1:yy)
#     if (TRUE %in% (SPIarray[x, y, ] > -9999))
#     {
#       CDHcor[x, y] = cor.test(SPIarray[x,y,rix], (-TSI)[x,y,rix])$estimate #which(monthhw[x,y,]>0)]
#       pcor[x,y]=cor.test(SPIarray[x,y,rix], (-TSI)[x,y,rix])$p.value
#     }
#   info = sprintf("已完成 %d%%", round(x*100/xx))   #查看循环进度
#   setTkProgressBar(pb, value = x*100/xx, title = sprintf("进度 (%s)",info),label = info) 
# }
# close(pb)
# CDHcorr = raster(CDHcor) #转为r层
# extent(CDHcorr) = c(72, 136, 18, 54) #确定范围
# crs(CDHcorr) = crs(r) #确定投影
# names(CDHcorr) = names(CDHcor) #名字

#####
#save.image("C:/Users/Yaoying/Desktop/2108论文收尾/211001定稿/compound_draft/复合计算_02_相关性分析结果.RData")


#02-02 Copula fitting
#####
vindex = dimnames(TX)[[3]]
mds=levels(as.factor(substring(vindex,5,8)))
#样本选择的滑动窗口长度
flowwindows=rep(-15:15,each=LL)
#取交集日
iflagarray = array(dim=dim(SPIarray)) #用于储存交集旗子
for(i in 1:days)
  iflagarray[,,i] = (SPIarray[,,i] <= (spithr1) ) & (-TSI[,,i] <= (hwthr) )

coptype = array(dim=c(xx,yy,length(mds)))
pCDHdi = array(dim=dim(SPIarray))
qCDHdi = array(dim=dim(SPIarray))
Q = array(dim=c(xx,yy,length(mds)))

pb = tkProgressBar(title="进度",label="已完成 %", min=0, max=100, initial = 0, width = 300)
for (x in 1:xx)
{
  for (y in 1:yy)
    if (TRUE %in% (SPIarray[x, y, ] > -9999))
    {
      #周围8个点，出界自动剔除
      xs=c(x,ifelse( x+1<xx & x+1>0 , x+1 , NA),ifelse( x-1<xx & x-1>0 , x-1 , NA)   )
      xs=subset(xs,xs!="NA")
      ys=c(y,ifelse( y+1<yy & y+1>0 , y+1 , NA),ifelse( y-1<yy & y-1>0 , y-1 , NA)   )
      ys=subset(ys,ys!="NA")
      
      for(i in 1:length(mds))
      {
        nx=which(substring(dimnames(TX)[[3]],5,8)==mds[i]) #取出全部60年长度序列
        nxs=nx+flowwindows #时间滑动
        nxs=subset(nxs,nxs<=days & nxs>0) #剔除超出长度范围的日数
        valueindex=nxs 
        
        pSPI=pnorm( as.vector(SPIarray[xs,ys,valueindex]) ) 
        pTSI=pnorm( as.vector(-TSI[xs,ys,valueindex]) ) #这里已经加了负号！
        iflags = as.vector( iflagarray[xs,ys,valueindex] )
        
        iflag = which( iflags == 1  ) #样本中所有交集发生日
        if(length(iflag) < 5)
          next
        
        q0 = (length(iflag))/length(iflags)  #计算基准发生率
        Q[x,y,i] = q0 
        
        valueind = which(iflagarray[x,y,nx] == 1) #nx中交集日的位置
        
        if( length( valueind ) > 0 ) #如果60个点中有交集日
        {
          #找到待求日在arx中的位置   
          finalIndex = match( SPIarray[x,y,nx[valueind]] ,(as.vector(SPIarray[xs,ys,valueindex]))[iflag] ) #交集日在子样本中的位置
          a1 = pSPI[iflag]/pnorm(spithr1)
          a2 = pTSI[iflag]/pnorm(spithr1)
          arx = cbind( a1 , a2 ) #用于拟合的样本[pTSI本身带负号][#注意有条件！]
          
          selectedCopula = BiCopSelect(a1,a2,familyset=c(3,4,5),rotations=F) #3,4,5 : "Clayton","Gumbel","Frank"
          coptype[x,y,i]=selectedCopula[[1]] #copula函数代号
          
          frankCopula = BiCopSelect(a1,a2,familyset=5,rotations=F) #实际上直接选5
          parm = frankCopula$par #获得参数
          mycopula = frankCopula(parm, dim = 2) #得到copula函数
          
          ybpcdhi = pCopula(arx, mycopula) #计算所有有值日的pdf值用于计算经验概率
          pcdhi = ybpcdhi[finalIndex]
          
          qCDHdi[x,y,nx[valueind]] = pcdhi*q0
          pCDHdi[x,y,nx[valueind]] = pcdhi
        }
      }
    }
   
  
  info = sprintf("已完成 %d%%", round(x*100/xx))   #查看循环进度
  setTkProgressBar(pb, value = x*100/xx, title = sprintf("进度 (%s)",info),label = info) 
}
close(pb)

#y=1-x 求指数
CDHdi = 1 - pCDHdi
#####
#warings: In cor(x[(x[, 1] < 0) & (x[, 2] < 0), ]) : 标准差为零

#03 Testing
#####
#测试结果准确性
x=28;y=53
colr=colorRampPalette(brewer.pal(9,"Blues")[1:9])(12)

ind0=which(is.na(pCDHdi[x,y,])==FALSE )
ind1=which(CDHdi[x,y,]> summary(CDHdi[x,y,])[2])
ind2=which(CDHdi[x,y,]> summary(CDHdi[x,y,])[3])
ind3=which(CDHdi[x,y,]> summary(CDHdi[x,y,])[4])
ind4=which(CDHdi[x,y,]> summary(CDHdi[x,y,])[5] )

xxx=SPIarray[x,y,ind0]
yyy=-TSI[x,y,ind0]
zzz=CDHdi[x,y,ind0]
plot(xxx,yyy,type="p",col="orange",xlab="SPI",ylab="-STI",pch=16,mgp = c(1.3, 0.4, 0),main=paste0(x,",",y))

xxx=SPIarray[x,y,ind1]
yyy=-TSI[x,y,ind1]
zzz=CDHdi[x,y,ind1]
points(xxx,yyy,type="p",col=colr[2],xlab="SPI",ylab="-STI",pch=16,mgp = c(1.3, 0.4, 0))

xxx=SPIarray[x,y,ind2]
yyy=-TSI[x,y,ind2]
zzz=CDHdi[x,y,ind2]
points(xxx,yyy,type="p",col=colr[5],xlab="SPI",ylab="-STI",pch=16,mgp = c(1.3, 0.4, 0))

xxx=SPIarray[x,y,ind3]
yyy=-TSI[x,y,ind3]
zzz=CDHdi[x,y,ind3]
points(xxx,yyy,type="p",col=colr[8],xlab="SPI",ylab="-STI",pch=16,mgp = c(1.3, 0.4, 0))

xxx=SPIarray[x,y,ind4]
yyy=-TSI[x,y,ind4]
zzz=CDHdi[x,y,ind4]
points(xxx,yyy,type="p",col=colr[12],xlab="SPI",ylab="-STI",pch=16,mgp = c(1.3, 0.4, 0))

#####
#save.image("C:/Users/Yaoying/Desktop/2108论文收尾/211001定稿/compound_draft/复合计算_03_指标计算结果.RData")


#04-01 Judge Flags
#####
Flag <- array(NA,dim = c(xx,yy,days)) #第一阈值
for(i in st:days)
  Flag[,,i] = ( CDHdi[,,i] >= 0 )
#####

#04-02 Counting days
#####
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

#识别连续da天以上1
Flag[which(is.na(Flag)==TRUE)] = -9999 #空值置为-9999，最后调回空值
pb <- tkProgressBar(title="进度",label="已完成 %", min=0, max=100)
for(x in 1:xx)
{
  for(y in 1:yy)
    if (TRUE %in% (Flag[x,y,] >= -999))
    {
      Flag[ x,y, st-1 ]=0 
      for( i in st:days )
        if( (Flag[x,y,i] == 1 & Flag[x,y,i-1] != 1) )
        {
          duration=sumfun(Flag[x,y,],i)
          if(duration < da)
            Flag[x,y,i:(i+duration-1)] <- rep(0,duration) #注意(i+duration-1)括号别漏
        }
    }
  info<- sprintf("已完成 %d%%", round(x*100/xx))   #查看循环进度
  setTkProgressBar(pb, value = x*100/xx, title = sprintf("进度 (%s)",info),label = info)   
}
close(pb)

#合并相邻事件，先不考虑
Flag2=Flag

#计次
Flag2[ which( is.na(Flag2) == TRUE ) ] = -9999
Flag3= array(dim = c(xx, yy, dim(Flag2)[3]))
for (x in 1:xx)
  for (y in 1:yy)
    if (TRUE %in% (Flag2[x, y,] >= -999))
      for(i in days:st)
        if( Flag2[x,y,i] == 1 & Flag2[x,y,i-1] != 1 )
          Flag3[x,y,i] = 1

Flag[which(Flag == -9999)] = NA
Flag2[which(Flag2 == -9999)] = NA

#####

#05 Calculating indicators
#####
#计算逐年频次和、年总历时、年平均强度、年规模
CDHI=CDHdi*Flag2 #得到带值旗子   
CDHI[which(CDHI==0)]=NA #要把0值置为NA,否则强度值结果不对

Hy = substring(dimnames(TX)[[3]], 1, 4)
fhy = factor(Hy) #转化为年份因子
CDHfre = array(dim = c(xx, yy, length(levels(fhy)))) #定义年频次
CDHdur = array(dim = c(xx, yy, length(levels(fhy)))) #定义年历时
CDHstg = array(dim = c(xx, yy, length(levels(fhy)))) #定义年强度
CDHmag = array(dim = c(xx, yy, length(levels(fhy)))) #定义年规模

dimnames(CDHfre)[[3]] = levels(fhy)
dimnames(CDHdur)[[3]] = levels(fhy) 
dimnames(CDHstg)[[3]] = levels(fhy) 
dimnames(CDHmag)[[3]] = levels(fhy) 

for (x in 1:xx)
  for (y in 1:yy)
    if (TRUE %in% (SPIarray[x, y,] >= -999))
    {
      CDHfre[x, y,] = tapply(Flag3[x, y,], fhy, sum, na.rm = TRUE) #计算年频次总和
      CDHdur[x, y,] = tapply(Flag2[x, y,], fhy, sum, na.rm = TRUE) #计算年历时总和
      CDHstg[x, y,] = tapply(CDHI[x, y,], fhy, mean, na.rm = TRUE) #计算年强度均值(包含了全部的0值，即计算的是干旱总强度/365天)
      CDHmag[x, y,] = tapply(CDHI[x, y,], fhy, sum, na.rm = TRUE) #计算年规模
    }

#00-03计算某地多年平均值
CDHFRE = array(dim = c(xx, yy))
CDHDUR = array(dim = c(xx, yy))
CDHSTG = array(dim = c(xx, yy))
CDHMAG = array(dim = c(xx, yy))
for (x in 1:xx)
  for (y in 1:yy)
  {
    CDHFRE[x, y] = mean(CDHfre[x, y,], na.rm = TRUE)
    CDHDUR[x, y] = mean(CDHdur[x, y,], na.rm = TRUE)
    CDHSTG[x, y] = mean(CDHstg[x, y,], na.rm = TRUE)
    CDHMAG[x, y] = mean(CDHmag[x, y,], na.rm = TRUE)
  }

#00-04放入r层方便画图（需加载工作空间：r）
CDHFREr = raster(CDHFRE) #转为r层
CDHDURr = raster(CDHDUR)
CDHSTGr = raster(CDHSTG)
CDHMAGr = raster(CDHMAG)
extent(CDHFREr) = c(72, 136, 18, 54) #确定范围
extent(CDHDURr) = c(72, 136, 18, 54) #确定范围
extent(CDHSTGr) = c(72, 136, 18, 54) #确定范围
extent(CDHMAGr) = c(72, 136, 18, 54) #确定范围
crs(CDHFREr) = crs(r) #确定投影
crs(CDHDURr) = crs(r) #确定投影
crs(CDHSTGr) = crs(r) #确定投影
crs(CDHMAGr) = crs(r) #确定投影

#00-05 计算某年多地平均值（按权重求和）
r1 = raster(CDHfre[, , 1])   #随便取一天，取其空间范围
extent(r1) = c(72, 136, 18, 54)
crs(r1) = crs(r)
w = area(r1, weights = TRUE, na.rm = TRUE) #求各网格的面积权重
w = as.matrix(w) #转为矩阵

CDHFREh = array(dim = c(LL))
CDHDURh = array(dim = c(LL))
CDHSTGh = array(dim = c(LL))
CDHMAGh = array(dim = c(LL))

for (i in 1:LL)
{
  y1 = CDHfre[, , i] * w  #权重×数值得新的数值
  CDHFREh[i] = sum(y1, na.rm = T)  #得到均值
  y2 = CDHdur[,,i] * w
  CDHDURh[i] = sum(y2, na.rm = T)
  y3 = CDHstg[,,i] * w
  CDHSTGh[i] = sum(y3, na.rm = T)
  y4 = CDHmag[,,i] * w
  CDHMAGh[i] = sum(y4, na.rm = T)
}

#另:覆盖面积DRper计算
#DRper相当于某年多地平均值
CDH1=CDHfre[,,1]*0+1 #每个格点定为1
CDHper=NULL
for( i in 1:length(CDHfre[1,1,]) )
{
  pt=which(CDHfre[,,i]>0) #有事件发生的点位置
  CDHper=c(CDHper,sum(CDH1[pt]*w[pt],na.rm=TRUE)) #发生点位*面积权重
}
names(CDHper)=levels(fhy)
#####

#06 Defining Class and Methods
#####
#定义出图类方法
pbt = function(x,...) UseMethod("pbt")

#05-01 空间分布图
pbt.SD = function(objx,ctype)
{
  colr = rev(colorRampPalette(brewer.pal(11, "RdYlGn"))(100)) #颜色
  vname = as.character(substitute(objx)) #提取变量名
  
  png(
    file = paste("01SD_",ctype, "_", vname, ".png", sep = ""),
    width = 1800,
    height = 1400,
    res = 72 * 3
  )
  plot(
    objx,
    xaxt = "n",
    yaxt = "n",
    main = paste("SD_", ctype, "_", vname, sep = ""),
    col = colr
  )
  axis(1, xL, xlabel)
  axis(2, yL, ylabel)
  lines(Chinasp)
  lines(provincesp)
  dev.off()
}

#05-02 时间变化分析
pbt.TS = function(objx,ctype) #DRFREh
{
  vname = as.character(substitute(objx)) #提取变量名
  #02-01 计算FRE突变点
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
  
  #02-01-02
  png(
    file = paste("02TS_",ctype, "_",vname,  ".png", sep = ""),
    width = 2500,
    height = 1800,
    res = 72 * 3
  )
  
  plot(c(s:e),objx,type="o",pch=20,xlab="",main= paste("TS_",ctype,"_",vname,sep=""))
  
  legend("bottomleft", 
         paste0("avg,max,min:", 
                round(mean(objx,na.rm = TRUE),4), ",",
                round(max(objx,na.rm = TRUE),4), "(", which.max(objx)+s-1,"),",
                round(min(objx,na.rm = TRUE),4), "(", which.min(objx)+s-1,")" ) 
         ,bty="n" )
  
  #02-01-02-a
  if(rst[1]!=0)
  {
    lines(FREt, FREy, lwd=3,col="red")
    points(FREt[2],FREy[2],pch=17,cex=2)
    legend("bottom",as.character(rst[1,1]),bty="n")
    legend("topright", paste0(round(rst[2,1],4),",",round(rst[2,2],4)),bty="n")
  }
  #02-01-02-b
  if(rst[1]==0)
  {
    lines(tendx,tendy,lwd=3,col="brown")
    text(1990,max(objx),paste("linear Fitting: y=",round(tend[[1]],5),"x+",round(tend[[3]],5), ",p-value=",round(tend[[2]],5),sep=""))
    text(1990,quantile(objx,0.97),paste("M-K Trend Test: tau=",round(mktau,5),",p-value=",round(mkpvalue,5),sep=""))
  }
  dev.off()
}

#05-03 变化率空间分布
pbt.TSTD = function(objx,ck,ctype) #objx=DRfre, ck = DRFREh 
{
  vname = as.character(substitute(objx)) #提取变量名
  #计算突变点
  rst=bp1(s,e,ck)
  #如果有突变点
  if(rst[1]!=0)
  {
    #Step1 计算两段变化率
    bp=rst[1]-s+1 #突变点序号
    tend1=tend2=array(dim=c(xx,yy))
    pfre1=pfre2=array(dim=c(xx,yy))
    #第一段
    for(x in 1:xx)
      for(y in 1:yy)
        if (TRUE %in% (objx[x,y,]>=0))
        {
          tend1[x,y]=tendency(objx[x,y,1:bp])[[1]]
          pfre1[x,y]=classtype(tendency(objx[x,y,1:bp])[[2]])
        }
    #第二段
    for(x in 1:xx)
      for(y in 1:yy)
        if (TRUE %in% (objx[x,y,]>=0))
        {
          tend2[x,y]=tendency(objx[x,y,(bp+1):(e-s+1)])[[1]]
          pfre2[x,y]=classtype(tendency(objx[x,y,(bp+1):(e-s+1)])[[2]])
        }
    
    #Step2 放入r层
    tend1r = raster(tend1) #转为r层
    extent(tend1r) = c(72, 136, 18, 54) #确定范围
    crs(tend1r) = crs(r) #确定投影
    
    tend2r = raster(tend2) #转为r层
    extent(tend2r) = c(72, 136, 18, 54) #确定范围
    crs(tend2r) = crs(r) #确定投影
    
    #Step3画显著性水平
    #3-2位置点计算
    yind=xind=pfre1ind=pfre2ind=NULL
    for(x in 1:xx)
      for(y in 1:yy)
      {
        yind=c(yind,71.75+0.5*y) #经度
        xind=c(xind,54.25-0.5*x) #纬度
        pfre1ind=c(pfre1ind,pfre1[x,y])
        pfre2ind=c(pfre2ind,pfre2[x,y])
      }
    
    #Step4 画格点图
    colr1=rev(colorRampPalette(brewer.pal(7,"GnBu")[1:7])(100)) 
    colr2=colorRampPalette(brewer.pal(9,"OrRd")[1:9])(100)
    
    tend1col = cf(c(colr1,"white",colr2),whitesite="white",tend1r,z0=0) 
    tend2col = cf(c(colr1,"white",colr2),whitesite="white",tend2r,z0=0) 
    
    #tend1
    png(
      file = paste("03TSTD_tend1_",rst[1],"_",ctype,"_", vname, ".png", sep = ""),
      width = 1800,
      height = 1400,
      res = 72 * 3
    )
    plot(
      tend1r,
      xaxt = "n",
      yaxt = "n",
      main = paste("TSTD_tend1_",rst[1],"_",ctype,"_", vname,  sep = ""),
      col=tend1col
    )
    axis(1, xL, xlabel)
    axis(2, yL, ylabel)
    lines(Chinasp)
    lines(provincesp)
    #
    points(yind,xind,pch=pfre1ind,cex=0.4) #点显著性
    points(126,19,pch=4,cex=1.3) #显著性水平标注
    text(130,19.1,"α=0.05",cex=1.2) #显著性水平标注
    dev.off()
    
    #tend2
    png(
      file = paste("03TSTD_tend2_",rst[1],"_",ctype,"_", vname, ".png", sep = ""),
      width = 1800,
      height = 1400,
      res = 72 * 3
    )
    plot(
      tend2r,
      xaxt = "n",
      yaxt = "n",
      main = paste("TSTD_tend2_",rst[1],"_",ctype,"_", vname, sep = ""),
      col=tend2col
    )
    axis(1, xL, xlabel)
    axis(2, yL, ylabel)
    lines(Chinasp)
    lines(provincesp)
    #
    points(yind,xind,pch=pfre2ind,cex=0.4) #点显著性
    points(126,19,pch=4,cex=1.3) #显著性水平标注
    text(130,19.1,"α=0.05",cex=1.2) #显著性水平标注
    dev.off()
  } else
    print("No break point")
}

#05-04 M-K趋势检验
pbt.MKSD = function(objx,ctype) #objx=DRfre
{
  vname = as.character(substitute(objx)) #提取变量名
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
  
  #画图
  colr=rev(colorRampPalette(brewer.pal(11, "RdYlBu"))(100)) #颜色
  colr[50]="white"
  MKcol =  cf(colr,whitesite="white",z,z0=0) 
  
  png(
    file = paste("08MKSD_",ctype, "_",vname, ".png", sep = ""),
    width = 1800,
    height = 1400,
    res = 72 * 3
  )
  plot(
    z,
    xaxt = "n",
    yaxt = "n",
    main = paste("MKSD_", ctype, "_",vname ,sep = ""),
    col=MKcol 
  )
  axis(1, xL, xlabel)
  axis(2, yL, ylabel)
  lines(Chinasp)
  lines(provincesp)
  
  points(yind,xind,pch=sigind,cex=0.4) #点显著性
  points(126,19,pch=4,cex=1.3) #显著性水平标注
  text(130,19.1,"α=0.05",cex=1.2) #显著性水平标注
  
  dev.off()
}  


#####

#07 Draft drawing
#####
#06-01 SD
pbt.SD(CDHFREr,ctype)
pbt.SD(CDHDURr,ctype)
pbt.SD(CDHSTGr,ctype)
pbt.SD(CDHMAGr,ctype)

#06-02 TS
pbt.TS(CDHFREh,ctype)
pbt.TS(CDHDURh,ctype)
pbt.TS(CDHSTGh,ctype)
pbt.TS(CDHMAGh,ctype)
pbt.TS(CDHper,ctype)

#06-03 TSTD
pbt.TSTD(CDHfre,CDHFREh,ctype)
pbt.TSTD(CDHdur,CDHDURh,ctype)
pbt.TSTD(CDHstg,CDHSTGh,ctype)
pbt.TSTD(CDHmag,CDHMAGh,ctype)

#06-04 MKSD
pbt.MKSD(CDHfre,ctype)
pbt.MKSD(CDHdur,ctype)
pbt.MKSD(CDHstg,ctype)
pbt.MKSD(CDHmag,ctype)

#####

save.image("C:/Users/Yaoying/Desktop/2108论文收尾/211001定稿/compound_draft/复合计算_final.RData")

