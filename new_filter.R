require(RCurl)
require(tidyverse)
require(magrittr)
require(knitr)
require(kableExtra)
require(data.table)
ObjID <- "1237680117417115655" # star
# ObjID <- "1237667211601248290"
# SQL that downloads some info on the chosen target from SDSS.
# ObjID from SDSS specifies the target
targetSqlQuery <- paste("SELECT top 10 ra, dec, psfmag_u, psfmag_g, psfmag_r, psfmag_i, psfmag_z FROM star WHERE ObjID = ", ObjID)

# downloads target data
# dataframe target has necessary info
targetSqlQuery <- str_squish(targetSqlQuery)
urlBase <- "http://skyserver.sdss.org/dr15/SkyserverWS/SearchTools/SqlSearch?"
X <- getForm(urlBase, cmd = targetSqlQuery, format = "csv")
target <- read.table(text = X, header = TRUE, sep = ",", dec = ".", comment.char = "#")

# sets some variables for convenience. Last two are the field sizes
# ra.size is automatically adjusted for each target depending on its dec
# M is the maximum colour difference
# resol is important for gauging crowded references
# dynamic range is to prevent saturation of either target or reference
u <- target$psfmag_u
g <- target$psfmag_g
r <- target$psfmag_r
i <- target$psfmag_i
z <- target$psfmag_z
ra <- target$ra
dec <- target$dec
M <- 0.1 # 0.1
dec.size = 10/60 # 0.167
dec.super = 12/60 # 0.20
ra.size = dec.size / cos(dec*pi/180)
ra.super = dec.super / cos(dec*pi/180)
resol <- 0.003
dynamic_range <- 2 # 5
crowd_mag_limit <- 5
'%+%' <- function(x,y) paste(x,y,sep="")

# SQL that counts objects in reference area
mySqlQuery <- str_glue(
  "SELECT objID, ra, dec, psfmag_u, psfmag_g, psfmag_r, psfmag_i, psfmag_z
FROM photoObj
WHERE (ra between ({ra - ra.size}) AND ({ra + ra.size})
OR ra BETWEEN ({360 + ra - ra.size}) AND ({360 + ra + ra.size})
OR ra BETWEEN ({-360 + ra - ra.size}) AND ({-360 + ra + ra.size}))
AND dec BETWEEN ({dec - dec.size}) AND ({dec + dec.size})
AND clean = 1
AND (calibStatus_r & 1) != 0")

# reads in data from SDSS.
# just counting how many objects in field
# doesn't care about mags
mySqlQuery <- str_squish(mySqlQuery)
X <- getForm(urlBase, cmd = mySqlQuery, format = "csv")
in_frame <- read.table(text = X, 
                       header = TRUE, 
                       sep = ",", 
                       dec = ".", 
                       comment.char = "#",
                       colClasses = c("character", rep("numeric",7)))
Object_Count <- dim(in_frame)[1]
Magnitude_Range <- in_frame %>% 
  filter(between(psfmag_r, r - dynamic_range, r + dynamic_range)) 
Magnitude_Count <- dim(Magnitude_Range)[1]
Colour_Range <- Magnitude_Range %>% 
  filter(between((psfmag_g - psfmag_r) , g - r - M, g - r + M)) %>% 
  filter(between((psfmag_r - psfmag_i), r - i - M, r - i + M))
Colour_Count <- dim(Colour_Range)[1]

# SQL query that downloads data from SDSS for objects
# potentially in the same field as the target
mySqlQuery1 <-str_glue(
  "SELECT objID, ra, dec, psfmag_u, psfmag_g, psfmag_r, psfmag_i, psfmag_z\n
FROM photoObj\n
WHERE (ra between ({ra - ra.size}) AND ({ra + ra.size})\n
OR ra BETWEEN ({360 + ra - ra.size}) AND ({360 + ra + ra.size})\n
OR ra BETWEEN ({-360 + ra - ra.size}) AND ({-360 + ra + ra.size})\n)
AND dec BETWEEN ({dec - dec.size}) AND ({dec + dec.size})\n
AND psfmag_r BETWEEN {r - dynamic_range} AND { r + dynamic_range} \n
AND (psfmag_g - psfmag_r) BETWEEN ({g - r - M}) AND ({g - r + M})\n
AND (psfmag_r - psfmag_i) BETWEEN ({r - i - M}) AND ({r - i + M})\n
AND clean = 1\n
AND (calibStatus_r & 1) != 0")

# reads in data from SDSS.
# dataframe called A has all the details
mySqlQuery <- str_squish(mySqlQuery1)

X <- getForm(urlBase, cmd = mySqlQuery, format = "csv")
A <- read.table(text = X, 
                header = TRUE, 
                sep = ",", 
                dec = ".", 
                comment.char = "#",
                colClasses = c("character", rep("numeric",7)))
Object_Count_mags <- dim(A)[1]

# wrap-around for targets near 0 RA
A$ra <- if_else(A$ra - target$ra > 180, A$ra - 360, A$ra) 
# wrap-around for targets near 360 RA
A$ra <- if_else(target$ra - A$ra > 180, A$ra + 360, A$ra)


# SQL query that downloads data from SDSS for objects
# potentially crowding the references

crowdingSqlQuery <- str_glue(
  "SELECT objID, ra, dec, psfmag_r FROM photoObj
WHERE (ra between ({ra - ra.size}) AND ({ra + ra.size})
OR ra between ({360 + ra - ra.size}) AND ({360 + ra + ra.size})
OR ra between ({-360 + ra - ra.size}) AND ({-360 + ra + ra.size}))
AND dec between ({dec - dec.size}) AND ({dec + dec.size})
AND psfmag_r < {r + crowd_mag_limit} AND clean = 1
AND (calibStatus_r & 1) != 0")
crowdingSqlQuery <- str_squish(crowdingSqlQuery)
crowd <- getForm(urlBase, cmd = crowdingSqlQuery, format = "csv")
crowdA <- read.table(text = crowd, 
                     header = TRUE, 
                     sep = ",", 
                     dec = ".", 
                     comment.char = "#",
                     colClasses = c("character", rep("numeric",3)))
# wrap-around for crowding near 0 RA
crowdA$ra <- if_else(crowdA$ra - target$ra > 180, crowdA$ra - 360, crowdA$ra) 
# wrap-around for crowding near 360 RA
crowdA$ra <- if_else(target$ra - crowdA$ra > 180, crowdA$ra + 360, crowdA$ra)


# function to check for crowding
crowd <- function(x, y, z, a, b, c) {
  sum(abs(x - a) < resol
      & abs(y - b) < resol
      & c - z < dynamic_range)
}


# function to calculate rating. Uses Oisin's routine
rating <- function(gr, rr, ir) {
  gt <- g
  rt <- r
  it <- i
  delta.CS <- (gt - rt) - (gr - rr)
  delta.CL <- (rt - it) - (rr - ir)
  RS <- 1 - abs(delta.CS / M)
  RL <- 1 - abs(delta.CL / M)
  RS * RL
}

# looks for crowding
for (q in dim(A)[1]:1) {
  if (crowd(A$ra[q], A$dec[q], A$psfmag_r[q], crowdA$ra, crowdA$dec, crowdA$psfmag_r) > 1) A <- A[-q, ]
}

object_count_crowding <- dim(A)[1]

# calculate ratings for each potential reference
ratings <- rating(A$psfmag_g, A$psfmag_r, A$psfmag_i)

# add ratings to the data frame
A <- cbind(A, ratings)

##########################
# finds all intersection points for each pair of potential references
##########################
A_coords <- A %>% arrange(ra) %>% select(objID, ra, dec)

swap_ij <- function(u, v){
  i <- v
  j <- u
}

int_pts_finder <- function(i, j) {
  if (i < j) {
    if (abs(A_coords$dec[i] - A_coords$dec[j]) < dec.size) {
      if (abs(A_coords$ra[i] - A_coords$ra[j]) < ra.size)
      {
        ifelse(
          A_coords$dec[i] > A_coords$dec[j],
          z <- data.frame(
            int_ra = c(A_coords$ra[i] + ra.size / 2,
                       A_coords$ra[j] - ra.size / 2),
            int_dec = c(A_coords$dec[j] + dec.size / 2,
                        A_coords$dec[i] - dec.size / 2),
            objID_i = c(A_coords$objID[i], A_coords$objID[i]),
            objID_j = c(A_coords$objID[j], A_coords$objID[j])
          ),
          z <- data.frame(
            int_ra = c(A_coords$ra[i] + ra.size / 2,
                       A_coords$ra[j] - ra.size / 2),
            int_dec = c(A_coords$dec[j] - dec.size / 2,
                        A_coords$dec[i] + dec.size / 2),
            objID_i = c(A_coords$objID[i], A_coords$objID[i]),
            objID_j = c(A_coords$objID[j], A_coords$objID[j])
          )
        )
        return(z)
      }
    }
  }
}       
index_matrix <- expand.grid(1:dim(A)[1], 1:dim(A)[1])
names(index_matrix) <- c("i", "j")
int.pts <- pmap(index_matrix, int_pts_finder) %>% 
  rbindlist() %>% 
  filter(between(int_ra, ra-ra.size/2, ra+ra.size/2), 
         between(int_dec, dec-dec.size/2, dec+dec.size/2)) %>%
  distinct()
#####################################



# function that returns the score for each intersection point
score1 <- function(X, Y) {
  B <- A[abs(A$ra - X) <= ra.size / 2 + 0.001 & abs(A$dec - Y) <= dec.size / 2 + 0.001, ]
  sum(B$ratings)
}

# function that returns a dataframe with all references
# in a FOV defined by an intersection point
score2 <- function(X, Y) {
  A[abs(A$ra - X) <= ra.size / 2 + 0.001 & abs(A$dec - Y) <= dec.size / 2 + 0.001, ]
}

# calculates the score for each intersection point and orders them
int.pts$score <- mapply(score1, int.pts$int_ra, int.pts$int_dec)
int.pts <- int.pts[order(int.pts$score, decreasing = T), ]

# prints out best pointing
# usually a bunch of ties  but just picks the first one
max.index <- which(int.pts$score == max(int.pts$score))
final_pointing <- data.frame(
  ra = int.pts[max.index[1], ]$int_ra,
  dec = int.pts[max.index[1], ]$int_dec,
  score = int.pts[max.index[1], ]$score
)

# makes data frame of reference stars for best pointing
B <- score2(int.pts[max.index[1], "int_ra"], int.pts[max.index[1], "int_dec"])
B <- B[with(B, order(ra)), c(1:9)]



loci_vertical <- data.frame(x1 = ifelse(A$ra < target$ra, A$ra + ra.size/2, A$ra - ra.size/2), 
                            x2 = ifelse(A$ra < target$ra, A$ra + ra.size/2, A$ra - ra.size/2),
                            y1 = A$dec - dec.size/2,
                            y2 = A$dec + dec.size/2)

loci_horizontal <- data.frame(x1 = A$ra - ra.size/2,
                              x2 = A$ra + ra.size/2,
                              y1 = ifelse(A$dec < target$dec, A$dec + dec.size/2, A$dec - dec.size/2),
                              y2 = ifelse(A$dec < target$dec, A$dec + dec.size/2, A$dec - dec.size/2))

int.pts %>% ggplot(aes(x = -int_ra, y = int_dec)) + 
  geom_point(fill = "#000000", colour = "#000000", size = 0.8) + #intersection points in black
  geom_tile(data = target, 
            aes(x=-ra, y=dec), 
            width=2*ra.size, 
            height=2*dec.size, 
            alpha=0, size = 1.5,
            colour = "#009E73") + ### candidate zone in green
  geom_point(data = A, aes(-ra, dec), colour = "red") + ### potential references in red
  geom_segment(data = loci_horizontal, aes(x = -x1, y = y1, xend = -x2, yend = y2),
            colour="blue") + ### locii in blue
  geom_segment(data = loci_vertical, aes(x = -x1, y = y1, xend = -x2, yend = y2),
               colour="blue") + ### locii in blue
  geom_tile(data = target, 
            aes(x=-ra, y=dec), 
            width=ra.size, 
            height=dec.size, 
            alpha=0,
            size = 0.6,
            colour = "blue") + ### extra locus for candidate
  geom_point(data = target, 
             aes(-ra, dec), size=4, 
             colour="#009E73", fill = "#009E73") + ### target in green
  labs(caption="Plot Showing the Target (Green Circle), \nPotential Reference Stars (Red Circles), \nAssociated Locii (Blue Squares), \nand Points of Intersection (Black Dots). \nThe Candidate Zone is Shown by the Green Square") + 
  ylab("Dec") +
  xlab("RA") +
  theme_minimal()