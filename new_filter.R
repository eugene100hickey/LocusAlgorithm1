A_coords <- A %>% arrange(ra) %>% select(ra, dec)

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
                        A_coords$dec[i] - dec.size / 2)
          ),
          z <- data.frame(
            int_ra = c(A_coords$ra[i] + ra.size / 2,
                       A_coords$ra[j] - ra.size / 2),
            int_dec = c(A_coords$dec[j] - dec.size / 2,
                        A_coords$dec[i] + dec.size / 2)
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
  rbind_list() %>% 
  filter(between(int_ra, ra-ra.size/2, ra+ra.size/2), 
         between(int_dec, dec-dec.size/2, dec+dec.size/2)) %>%
  distinct()

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

int.pts %>% ggplot(aes(x = int_ra, y = int_dec)) + 
  geom_point() + 
  geom_tile(data = data.frame(ra=ra, dec=dec), 
            aes(x=ra, y=dec), 
            width=2*ra.size, 
            height=2*dec.size, 
            alpha=0, 
            colour = "blue") + 
  geom_point(data = A, aes(ra, dec), colour = "red") + 
  geom_tile(data = A, aes(ra, dec), 
            width=ra.size, 
            height=dec.size, 
            colour="red", 
            alpha=0)