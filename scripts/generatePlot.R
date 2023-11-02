# Function to generate a plot
generatePlot <- function(trait) {
  x <- 1:10
  y <- x^2
  
  data <- data.frame(x = x, y = y)
  
  ggplot(data, aes(x = x, y = y)) +
    geom_line() +
    labs(title = "Plot generated in a function")
}
