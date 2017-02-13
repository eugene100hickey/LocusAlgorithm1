# LOcus Algorithm in R
# Eugene
# February 2016

#################################
require(RCurl)
require(ggplot2)

#just in case a directory needs to be specified
#some files will get written here
setwd("C:/Users/Fujitsu PC/Desktop/Oisin/Locus_data")

ObjID = "1237658203434778700"
ObjID = "1237658492816130233"
#ObjID = "1237660613977899152"
#ObjID="1237674649391661089"
#ObjID="1237663543673815599"



#SQL that downloads some info on the chosen target from SDSS. 
#ObjID from SDSS specifies the target
targetSqlQuery = paste("select ra, dec, psfmag_u, psfmag_g, psfmag_r, psfmag_i, psfmag_z
                       from star
                       where ObjID = ", ObjID)

#downlaods target data
#creates target .csv file
#dataframe target has necessary info
targetSqlQuery = gsub(pattern="\n",replacement=" ",x=targetSqlQuery)
urlBase = "http://skyserver.sdss.org/dr12/SkyserverWS/SearchTools/SqlSearch?"
X = getForm(urlBase, cmd = targetSqlQuery, format = "csv")
write(X,file="SDSStarget.csv")
target = read.table("SDSStarget.csv", header=TRUE, sep=",", dec=".", comment.char="#")

#sets some variables for convenience. Last two are the field sizes
#ra.size is automatically adjusted for each target depending on its dec
#M is the maximum colour difference
g = target$psfmag_g
r = target$psfmag_r
i = target$psfmag_i
ra = target$ra
dec = target$dec
M = 0.1
dec.size = 0.25
ra.size = dec.size / cos(dec*pi/180)
resol = 0.003
dynamic_range = 5


#SQL query that downloads data from SDSS for objects 
#potentially in the same field as the target
mySqlQuery = paste("select objID, ra, dec, psfmag_u, psfmag_g, psfmag_r, psfmag_i, psfmag_z
                   from star
                   where
                   ra between (", ra-ra.size, ") and (", ra+ra.size, ")
                   and dec between (", dec-dec.size, ") and (", dec+dec.size, ")
                   and psfmag_r between ", r-2, "and ", r+2, " 
                   and (psfmag_g - psfmag_r) between (", g-r-0.1, ") and (", g-r+0.1, ")
                   and (psfmag_r - psfmag_i) between (", r-i-0.1, ") and (", r-i+0.1, ")
                   --and psfmagerr_u < 0.05 and psfmagerr_g < 0.05 and psfmagerr_r < 0.05 and psfmagerr_i < 0.05 and psfmagerr_z < 0.05
                   and ((flags & 0x10000000) != 0)       -- detected in BINNED1
                   AND ((flags & 0x8100000c00a4) = 0)     -- not EDGE, NOPROFILE, PEAKCENTER,
                   -- NOTCHECKED, PSF_FLUX_INTERP, SATURATED, or BAD_COUNTS_ERROR
                   AND ((flags & 0x400000000000) = 0) 
                   -- not DEBLEND_NOPEAK or small PSF error
                   -- (substitute psfmagerr in other band as appropriate)
                   AND (((flags & 0x100000000000) = 0) or (flags & 0x1000) = 0)", sep=" ")

#reads in data from SDSS. Creates file on disk
#dataframe called A has all the details
mySqlQuery = gsub(pattern="\n",replacement=" ",x=mySqlQuery)
X = getForm(urlBase, cmd = mySqlQuery, format = "csv")
write(X,file="SDSSsample.csv")
A = read.table("SDSSsample.csv", header=TRUE, sep=",", dec=".", comment.char="#")
A$objID = as.character(A$objID)



#SQL query that downloads data from SDSS for objects 
#potentially crowding the target
crowdingSqlQuery = paste("select objID, ra, dec, psfmag_r from star
                         where
                         ra between (", ra-ra.size, ") and (", ra+ra.size, ")
                         and dec between (", dec-dec.size, ") and (", dec+dec.size, ")
                         and psfmag_r < ", r+7, sep=" ")
crowdingSqlQuery = gsub(pattern="\n",replacement=" ",x=crowdingSqlQuery)
crowd = getForm(urlBase, cmd = crowdingSqlQuery, format = "csv")
write(crowd,file="SDSScrowd.csv")
crowdA = read.table("SDSScrowd.csv", header=TRUE, sep=",", dec=".", comment.char="#")

#function to check for crowding
crowd = function(x, y, z, a, b, c){
  sum(abs(x - a) < resol
      & abs(y - b) < resol 
      & c-z < dynamic_range)
}

#functions that find intersection points
X1 = function(x1, x2){min(x1, x2)+ra.size/2}
Y1 = function(y1, y2){max(y1, y2)-dec.size/2}
X2 = function(x1, x2){max(x1, x2)-ra.size/2}
Y2 = function(y1, y2){min(y1, y2)+dec.size/2}

#function to calculate rating. Uses Oisín's routine
rating = function(gr, rr, ir){
  gt = g
  rt = r
  it = i
  delta.CS = (gt-rt) - (gr-rr)
  delta.CL = (rt-it) - (rr-ir)
  RS = 1 - abs(delta.CS/M)
  RL = 1 - abs(delta.CL/M)
  RS*RL
}

#looks for crowding
for(q in dim(A)[1]:1){
  if(crowd(A$ra[q], A$dec[q], A$psfmag_r[q], crowdA$ra, crowdA$dec, crowdA$psfmag_r)>1) A = A[-q,]
}


#calculate ratings for each potential reference
ratings = rating(A$psfmag_g, A$psfmag_r, A$psfmag_i)

#add ratings to the data frame
A = cbind(A, ratings)

#finds all intersection points for each pair of potential references
x1 = as.vector(sapply(A$ra, function(x){mapply(X1, x,A$ra)}))
x2 = as.vector(sapply(A$ra, function(x){mapply(X2, x,A$ra)}))
y1 = as.vector(sapply(A$dec, function(x){mapply(Y1, x,A$dec)}))
y2 = as.vector(sapply(A$dec, function(x){mapply(Y2, x,A$dec)}))

#creates dataframe with all these intersection points
x = c(x1, x1, x2, x2)
y = c(y1, y2, y1, y2)
int.pts = data.frame(x, y)
names(int.pts) = c("ra", "dec")


#function that returns the score for each intersection point
score1 = function(X, Y){
  B = A[abs(A$ra-X)<=ra.size/2+0.001 & abs(A$dec-Y)<=dec.size/2+0.001,]
  sum(B$ratings)
}

#function that returns a dataframe with all references 
#in a FOV defined by an intersection point
score2 = function(X, Y){
  A[abs(A$ra-X)<=ra.size/2+0.001 & abs(A$dec-Y)<=dec.size/2+0.001,]
}

#calculates the score for each intersection point and orders them
int.pts$score = mapply(score1, int.pts$ra, int.pts$dec)
int.pts = int.pts[order(int.pts$score, decreasing = T), ]

#prints out best pointing
#usually a bunch of ties  but just picks the first one
max.index = which(int.pts$score == max(int.pts$score))
final_pointing = as.data.frame(c(ra=int.pts[max.index[1],]$ra, dec=int.pts[max.index[1],]$dec, score = int.pts[max.index[1],]$score))
final_pointing = as.data.frame(t(final_pointing))
View(final_pointing)

#prints data frame of reference stars for best pointing
B = score2(int.pts[max.index[1],"ra"], int.pts[max.index[1], "dec"])
B = B[with(B, order(ra)), c(1:9)]
View(B)

#point is the target ra/dec, not the chosen pointing
point = as.data.frame(c(ra, dec))

#draws a picture. Green dot is the target. 
#Only chosen references and their boxes are shown
g = ggplot(data=B, aes(x=-ra, y=dec))
g = g+geom_point(aes(color=-ratings),size=5)
g = g+geom_tile(data = final_pointing, aes(width=ra.size, height=dec.size, fill=F), alpha=0, colour="black")
g = g+geom_text(aes(label=rownames(B)), size=8, nudge_x=0.02, nudge_y=0.02, colour="red")
g = g+geom_point(data=final_pointing, size=4, colour="green", shape=1)
g = g+geom_point(data=target, size=4, colour="yellow", shape=1)

g
