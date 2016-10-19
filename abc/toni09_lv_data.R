# Generate the sample data (from figure 1 of Toni et al. 2009)

# Times of sample data points
tt <- c(1.2, 2.4, 3.9, 5.7, 7.5, 9.6, 11.9, 14.5)

# Predators (triangles in fig. 1)
predators <-c(0.5, 2.6, 1.55, 0.03, 1.15, 1.68, 1.1, 0.94)

# Prey (circles in fig. 1)
prey <- c(1.9, 0.63, 0.2, 0.3, 1.65, 1.15, 0.25, 2.94)

sample_data <- data.frame(time = tt, x = prey, y = predators)

write.csv(sample_data, "../data/toni09_lv_data.csv", row.names = FALSE)

