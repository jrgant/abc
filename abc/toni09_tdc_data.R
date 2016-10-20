
It <- c(1, 1, 3, 7, 6, 10, 13, 13, 14, 14, 17, 10, 6, 6, 4, 3, 1, 1, 1, 1, 0)
Rt <- c(0, 0, 0, 0 , 5, 7, 8, 13, 13, 16, 16, 24, 30, 31, 33, 34, 36, 36, 36, 36, 37)

tdc_data <- data.frame(day = 1:21, I = It, R = Rt)

write.csv(tdc_data, "../data/toni09_tdc_data.csv", row.names = FALSE)
