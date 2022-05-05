# background setting ----
setwd("C:/Users/Yaoying/Desktop/2108论文收尾/211001定稿/EF_fig1_program")
library("sp")
library("maptools")
library("rgdal")
library("raster")
library("tcltk")
library("RColorBrewer") 
Chinasp = readShapePoly("C:/Users/Yaoying/Desktop/hyj/china/China") #中国边界
provincesp = readShapePoly("C:/Users/Yaoying/Desktop/hyj/china/PROVINCE_region") #省份边界
load("C:/Users/Yaoying/Desktop/2108论文收尾/211001定稿/r层.RData") 
load("C:/Users/Yaoying/Desktop/2108论文收尾/211001定稿/SPI_TX.RData") #SPIarray
xx=72
yy=128
Longtirange = seq(72,136,0.5)
Lattirange = seq(18,54,0.5)

xL = seq(72, 136, 5)
xlabel = paste(xL, "°", sep = "")
yL = seq(18, 54, 5)
ylabel = paste(yL, "°", sep = "")


# read information ----
stationInfo = read.table("SURF_CHN_MUL_HOR_STATION.csv", header=F, sep=",")
stationID = stationInfo[,2]
stationLatti = stationInfo[,4]
stationLongti = stationInfo[,5]
snum = length(stationLatti)

# convert into real value ----
conv = function(x)
{
  xc = as.character(x)
  clen = nchar(xc)
  if (clen == 5)
  {
    x_min = as.numeric(substr( xc, 1, 3))
    x_sec = as.numeric(substr( xc, 4, 5))
  }else
  {
    x_min = as.numeric(substr( xc, 1, 2))
    x_sec = as.numeric(substr( xc, 3, 4))
  }
  x_out = x_min+x_sec/60  
  return(x_out)  
}

# calculating ----
LattiValue = rep(NA,snum)
LongtiValue = rep(NA,snum)
for ( i in 1:snum )
{
  LattiValue[i] = conv(stationLatti[i])
  LongtiValue[i] = conv(stationLongti[i])
}
  
# judge index ----
foundindex = function(vv)
{
  if( 0 %in% vv )
    ind = which(vv==0)
  else
  {
    groups = embed(vv,2)
    mutres = groups[,1]*groups[,2]
    ind = which(mutres < 0)
  }
  return(ind)
}

staLatti = rep(NA,snum)
staLongti = rep(NA,snum)
for ( i in 1:snum )
{
  Longtigap = Longtirange - LongtiValue[i]
  staLongti[i] = foundindex(Longtigap)
  Lattigap = Lattirange - LattiValue[i]
  staLatti[i] = foundindex(Lattigap)
}

# calculating ----
stationmap = array(0,dim=c(xx,yy)) # set blank map

for ( i in 1:snum )
{
  x = 73 - staLatti[i]
  y = staLongti[i]
  stationmap[x,y] = stationmap[x,y]+1
}

# classify ----
classifyType = function(v,a1,a2,a3,a4)
{
    c=NA
      if(v==1)
        c=a1 else 
          if(v==2)
            c = a2 else
              if(v==3)
                c = a3 else
                  if(v>=4)
                    c = a4
              
        return(c)
}

a1 = "#006d2c"
a2 = "#08519c"
a3 = "#810f7c"
a4 = "#bd0026"

  #ffeda0
  #fed976
  #feb24c
  #fd8d3c
  #fc4e2a
  #e31a1c
  #bd0026
  #800026  
  
sigmap = array(NA,dim=c(xx,yy))
for(x in 1:xx)
  for(y in 1:yy)
      sigmap[x,y] = classifyType(stationmap[x,y],a1,a2,a3,a4)
  

yind=xind=sigindex=NULL
for(x in 1:xx)
  for(y in 1:yy)
  {
    yind=c(yind,71.75+0.5*y) 
    xind=c(xind,54.25-0.5*x) 
    sigindex=c(sigindex,sigmap[x,y])
  }

  
# plotting ----
# 
# sizeindex = sigindex
# colindex = sigindex
# sizeindex[which(sigindex == a1)] = 0.6
# sizeindex[which(sigindex == a2)] = 0.7
# sizeindex[which(sigindex == a3)] = 0.8
# colindex[ which(sigindex == a1) ] = "#000000"  #525252
# colindex[ which(sigindex == a2) ] = "#000000" #252525
# colindex[ which(sigindex == a3) ] = "#000000" #000000
# bgindex = colindex

#fa9fb5
#f768a1
#dd3497
#ae017e
#7a0177

par(mar=c(1,1,1,1))
mapr = raster(SPIarray[,,1]) #转为r层
extent(mapr) = c(72, 136, 18, 54)
plot(mapr,
     ylim=c(17,55),xlim=c(72,136),
     xaxs="i",yaxs="i",
     xaxt="n",yaxt="n",
     xlab="",ylab="",
     mgp=c(0,0,0))

# axis(1, xL, xlabel)
# axis(2, yL, ylabel)
lines(provincesp)
points(yind,xind,pch=22, col=sigindex, bg=sigindex, cex=0.7,) 


