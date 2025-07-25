# Calibration plots

library(glue)
library(ggplot2)
library(writexl)

# ---- FUNCTION FOR MAKING BINS ----

#'@description Function which creates bins for calibration plot
#'
#'@param top_prob Vector of predicted top class probabilities
#'@param n_bins Number of bins to make
#'@return List containing the following:
#'\item{bin_ids}{Factor vector containing the interval in which the predicted prob belongs}
#'\item{bin_centers}{Vector of the midpoints of the bins}
#'\item{counts}{Vector of bin sizes}
make_calibration_bins <- function(top_prob, n_bins) {
  # Compute quantile breakpoints
  probs <- seq(0, 1, length.out = n_bins + 1)
  bin_cuts <- quantile(top_prob, probs = probs, na.rm = TRUE)
  
  # Remove duplicate cut points (can happen with tied probs)
  bin_cuts <- unique(bin_cuts)
  
  # Assign each probability to a bin
  bin_ids <- cut(top_prob, breaks = bin_cuts, include.lowest = TRUE)
  
  # Check how many predictions fall in each bin
  counts <- table(bin_ids)
  print(counts)
  
  # Compute bin centers (midpoint of cut edges)
  bin_centers <- (head(bin_cuts, -1) + tail(bin_cuts, -1)) / 2
  
  return(list(bin_ids = bin_ids, bin_centers = bin_centers, counts = counts))
}

# ---- READ DATA ----
dataset <- "vehicle"
method <- "mpbart"
n_bins <- 10 # number of bins for calibration plot
rds_path <- glue("C:/Users/matth/OneDrive/Bureaublad/msc_thesis/thesis_code/calibration_plots/calibration_{method}_{dataset}.rds")
data <- readRDS(rds_path)

# ---- PREPARE DATA FOR PLOTS ----

# combine labels and predicted probs across folds
labels_all <- unlist(lapply(data, function(x) x$labels))
pred_probs_all <- do.call(rbind, lapply(data, function(x) x$pred_probs))
rownames(pred_probs_all) <- NULL

# extract top class, top class probability and indicator for correct prediction
top_class <- apply(pred_probs_all, 1, which.max) - 1
top_prob <- apply(pred_probs_all, 1, max)
is_correct <- as.integer(top_class == labels_all)

# create bins and compute center and accuracy per bin
bins <- make_calibration_bins(top_prob = top_prob, n_bins = n_bins)
bin_centers <- bins$bin_centers
avg_pred_probs <- tapply(top_prob, bins$bin_ids, mean)
bin_accuracies <- tapply(is_correct, bins$bin_ids, mean)

# ---- CREATE CALIBRATION PLOT ----

df_plot <- data.frame(
  avg_pred_probs = as.numeric(avg_pred_probs),
  empirical_accuracies = as.numeric(bin_accuracies)
)

# save to excel
file_path <- glue("C:/Users/matth/OneDrive/Bureaublad/msc_thesis/thesis_code/calibration_plots/plot_data_{dataset}_{method}.xlsx")
write_xlsx(df_plot, path = file_path)

# create plot
ggplot(df_plot, aes(x = avg_pred_probs, y = empirical_accuracies)) +
  geom_line() +
  geom_point(size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Average Predicted Top-Class Probability (per bin)",
       y = "Empirical Accuracy",
       title = "Top-Class Calibration Plot") +
  theme_minimal()


