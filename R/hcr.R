tier_3 <- function(b, b40, f40) {
  b <- as.numeric(b)
  b40 <- as.numeric(b40)
  f40 <- as.numeric(f40)
  
  if(length(b) == 0 || is.na(b)) return(0)
  
  bb40 <- b / b40
  
  if(bb40 > 1) {
    return(f40)
  } else if(bb40 > 0.05 && bb40 <= 1) {
    return(f40 * (bb40 - 0.05) / (1-0.05))
  } else {
    return(0)
  }
}