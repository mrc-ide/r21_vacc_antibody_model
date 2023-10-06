get_season <- function(t = 1:365, g0, g1, g2, g3, h1, h2, h3, floor){
  
  # Fourier series
  t <- t / 365
  prediction <- data.frame(
    g0 = g0,
    g1 = g1 * cos(2 * pi * t * 1),
    g2 = g2 * cos(2 * pi * t * 2),
    g3 = g3 * cos(2 * pi * t * 3),
    h1 = h1 * sin(2 * pi * t * 1),
    h2 = h2 * sin(2 * pi * t * 2),
    h3 = h3 * sin(2 * pi * t * 3))
  prediction <- rowSums(prediction)
  prediction <- pmax(floor, prediction)
  
  return(prediction)
  
}

