library(readxl)
library(mclust)
library(dunn.test)
library(dbscan)
library(plotly)
library(dplyr)
library(ggplot2)
library(patchwork)

# Manually go through each individual brain by chaning ID
# This step is to semi-manually assign Vsx vs Rx identity by A-P position.
# The manual part is mostly to account for mounting imprecisions (
# therefore we couldn't have axes to be perfectly aligning to AP/DV/ML axes).

out <- data.frame()
id <- 7
file_path <- paste0("20231115_CantonS_405ArmMs_488GrhRt_555Vsx1RxGp_647OptixRb_", id, ".xls")
sn <- excel_sheets(file_path)
ol_pos <- read_xls(file_path, sheet = which(sn == "Position"), skip = 1)
colnames(ol_pos)[1:3] <- c("x", "y", "z")
ol_vr <- read_xls(file_path, sheet = which(sn == "Intensity Median Ch=2 Img=1"), skip = 1)
ol_op <- read_xls(file_path, sheet = which(sn == "Intensity Median Ch=3 Img=1"), skip = 1)
mcoptix <- Mclust(ol_op$`Intensity Median`, G = 2)
mcvr <- Mclust(ol_vr$`Intensity Median`, G = 2)

optix_pos <- which.max(mcoptix$parameters$mean)
vr_pos <- which.max(mcvr$parameters$mean)

H <- ks::Hns(x = as.matrix(ol_pos[ , c("x", "y", "z")]))
kde <- ks::kde(x = as.matrix(ol_pos[ , c("x", "y", "z")]), H = H, eval.points = as.matrix(ol_pos[ , c("x", "y", "z")]))
knn <- kNN(as.matrix(ol_pos[ , c("x", "y", "z")]), k = 5)
vr_epos <- sqrt((1 - mcvr$z[, vr_pos]) * mcoptix$z[ , optix_pos])
plot_ly(data = ol_pos[ , c("x", "y", "z")], x = ~x, y = ~y, z = ~z, color = vr_epos)

ident <- ifelse(vr_epos > 0.3, "Optix", "tbd")

## Inspect the shape of Optix domain and assign the cut between Vsx vs Rx
## This step is subjective but should be mostly robust due to the
## huge gap (the Optix domain) between Vsx and Rx
# ident[ident == "tbd"] <- ifelse(ol_pos$x[ident == "tbd"] > 25310, "Vsx", "Rx")
ident[ident == "tbd"] <- ifelse(ol_pos$x[ident == "tbd"] > 27830 & ol_pos$y[ident == "tbd"] > 28720 & ol_pos$y[ident == "tbd"] < 28845 & ol_pos$z[ident == "tbd"] < -90, "Vsx", "Rx")
# plot_ly(data = ol_pos[ , c("x", "y", "z")], x = ~x, y = ~y, z = ~z, color = mcvr$z[ , vr_pos])
plot_ly(data = ol_pos[ , c("x", "y", "z")], x = ~x, y = ~y, z = ~z, color = ident)

## After annotation, concatenate intensity information and assignment
## onto a master table
temp_out <- data.frame(
  x = ol_pos$x,
  y = ol_pos$y,
  z = ol_pos$z,
  vsx_rx = ol_vr$`Intensity Median`,
  optix = ol_op$`Intensity Median`,
  domain = ident,
  density = kde$estimate,
  knnd = rowMeans(knn$dist),
  sample = id
)
out <- rbind(out, temp_out)

# After everything is done, make domain identity categorical so we
# can assign the order of levels for aesthetic purposes.
out$domain <- factor(out$domain, levels = c("Vsx", "Optix", "Rx"))


# The brains were collected from wandering larvae, which represents a ~12 hr
# developmental window during which Rx weans down in the Rx domain (Bi/Omb
# persists but the antibody is incompatible for this experiment).
# We thus only keep animal 4, 6, and 7 that are at similar stages when
# Rx can still be confidently detected for further evaluation.
res <- out %>%
  filter(sample %in% c(4, 6, 7))

sumtbl <- res %>%
  group_by(sample, domain) %>%
  summarize(kdist = median(knnd), density = median(density))


# Re-code animal 4, 6, 7 to be larval brain 1, 2, 3 for aesthetic purposes.
res$sample <- factor(res$sample, levels = c(4, 6, 7), labels = paste("Larval brain", 1:3))


# Export quantification for future re-examination
write.csv(res, "quant.csv", row.names = FALSE, quote = TRUE)
res <- read.csv("quant.csv")

# Visualization code to reproduce the figures.
res %>%
  ggplot(aes(x = x, y = y, color = domain)) +
    geom_point() +
    facet_wrap(~sample, scales = "free") +
    scale_color_viridis_d() +
    theme_bw() +
  labs(
    title = "Domain Assignment",
    color = "Domain"
  ) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      strip.background = element_rect(fill = "white"),
      legend.position = "bottom"
    )

out %>%
  ggplot(aes(x = x, y = y, color = domain)) +
  geom_point() +
  facet_wrap(~sample, scales = "free") +
  scale_color_viridis_d() +
  theme_bw() +
  labs(
    title = "Domain Assignment",
    color = "Domain"
  ) +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    strip.background = element_rect(fill = "white"),
    legend.position = "bottom"
  )

ggsave("domain_assignment.pdf", width = 16, height = 10)

res %>%
  ggplot(aes(x = x, y = y, color = log(vsx_rx))) +
  geom_point() +
  facet_wrap(~sample, scales = "free") +
  scale_color_viridis_c() +
  theme_bw() +
  labs(
    title = "Vsx/Rx Fluorescence Intensity",
    color = "Log(Intensity)"
  ) +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    strip.background = element_rect(fill = "white")
  )

(
  res %>%
    ggplot(aes(x = x, y = y, color = log(density))) +
    geom_point() +
    facet_wrap(~sample, scales = "free") +
    scale_color_viridis_c(option = "A") +
    scale_x_reverse() +
    theme_bw() +
    labs(
      title = "Cell Density (Gaussian Kernel)",
      color = "Log(Density)"
    ) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      strip.background = element_rect(fill = "white"),
      plot.title = element_text(size = 44, face = "bold"),
      legend.key.size = unit(2, "cm"),
      legend.text = element_text(size = 20),
      legend.title = element_text(size = 24, face = "bold"),
      strip.text = element_text(size = 24),
      legend.position = "bottom"
    ) | ggplot(sumtbl, aes(x = domain, y = density, fill = domain)) +
    geom_boxplot() +
    geom_point(size = 3) +
    # geom_segment(x = 1, xend = 2, y = 3e-6, yend = 3e-6) +
    # annotate(geom = "text", x = 1.5, y = 3.1e-6, hjust = 0.5, label = "n.s.") +
    # geom_segment(x = 2, xend = 3, y = 2.6e-6, yend = 2.6e-6) +
    # annotate(geom = "text", x = 2.5, y = 2.7e-6, hjust = 0.5, label = "p = 0.018") +
    geom_segment(x = 1, xend = 3, y = 3.3e-6, yend = 3.3e-6) +
    annotate(
      geom = "text", x = 2, y = 3.5e-6, hjust = 0.5, size = 8,
      label = paste0(
        "n.s. (p = ",
        round(kruskal.test(x = sumtbl$density, g = sumtbl$domain)$p.value, digits = 3),
        ")"
        )
      ) +
    scale_fill_viridis_d() +
    scale_y_continuous(limits = c(0, 3.6e-6)) +
    labs(
      y = "Gaussian Kernel Density Estimate",
      fill = "Domain"
    ) +
    guides(fill = "none") +
    theme_bw() +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 20),
      text = element_text(size = 28)
    )
) + plot_layout(width = c(3, 1))
ggsave("density_domain.pdf", width = 20, height = 9)

res %>%
  ggplot(aes(x = x, y = y, color = log(density))) +
  geom_point() +
  facet_wrap(domain~sample, scales = "free") +
  scale_color_viridis_c(option = "A") +
  theme_bw() +
  labs(
    title = "Cell Density per domain",
    color = "Log(Density)"
  ) +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    strip.background = element_rect(fill = "white")
  )

ggsave("density_per_domain.pdf", width = 16, height = 12)

(
  res %>%
    ggplot(aes(x = x, y = y, color = knnd)) +
    geom_point() +
    facet_wrap(~sample, scales = "free") +
    scale_color_viridis_c(option = "G", direction = -1) +
    scale_x_reverse() +
    theme_bw() +
    labs(
      title = "Average Distance to Nearest Neighbors (k = 5)",
      color = "Distance"
    ) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      strip.background = element_rect(fill = "white")
    ) | ggplot(sumtbl, aes(x = domain, y = kdist, fill = domain)) +
    geom_boxplot() +
    geom_point(size = 3) +
    geom_segment(x = 1, xend = 2, y = 5.3, yend = 5.3) +
    annotate(geom = "text", x = 1.5, y = 5.4, hjust = 0.4, label = "n.s.") +
    geom_segment(x = 2, xend = 3, y = 5.8, yend = 5.8) +
    annotate(geom = "text", x = 2.5, y = 5.9, hjust = 0.5, label = "p = 0.037") +
    geom_segment(x = 1, xend = 3, y = 6.1, yend = 6.1) +
    annotate(geom = "text", x = 2, y = 6.3, hjust = 0.5, label = "p = 0.013") +
    scale_fill_viridis_d() +
    scale_y_continuous(limits = c(0, 6.5)) +
    labs(
      y = "Average Nearest Neighbor Distance\n(k = 5)",
      fill = "Domain"
    ) +
    guides(fill = "none") +
    theme_bw() +
    theme(
      axis.title.x = element_blank()
    )
) + plot_layout(width = c(3, 1))

ggsave("knn_domain.pdf", width = 16, height = 9)

res %>%
  ggplot(aes(x = x, y = y, color = log(knnd))) +
  geom_point() +
  facet_wrap(sample ~ domain, scales = "free") +
  scale_color_viridis_c(option = "G", direction = -1) +
  theme_bw() +
  labs(
    title = "Average Distance to Nearest Neighbors (k = 5)",
    color = "Distance"
  ) +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    strip.background = element_rect(fill = "white")
  )

res

## Hypothesis testing -- note that Dunn test is post-hoc while Kruskal-Wallis
## is omnibus. One shall not perform Dunn test if KW turns out insignificant.
## dunn.test() is used here because this function is implemented in a way that
## KW test is done first followed by a Dunn test, so it is empirically
## convenient, but pairwise comparison from Dunn Test should be
## ignored if KW is insignificant.
dunn.test(x = sumtbl$kdist, g = sumtbl$domain)
dunn.test(x = sumtbl$density, g = sumtbl$domain)

ggplot(sumtbl2, aes(x = domain, y = density, fill = domain)) +
  geom_boxplot() +
  geom_point(size = 3) +
  # geom_segment(x = 1, xend = 2, y = 3e-6, yend = 3e-6) +
  # annotate(geom = "text", x = 1.5, y = 3.1e-6, hjust = 0.5, label = "n.s.") +
  # geom_segment(x = 2, xend = 3, y = 2.6e-6, yend = 2.6e-6) +
  # annotate(geom = "text", x = 2.5, y = 2.7e-6, hjust = 0.5, label = "p = 0.018") +
  # geom_segment(x = 1, xend = 3, y = 3.3e-6, yend = 3.3e-6) +
  # annotate(geom = "text", x = 2, y = 3.4e-6, hjust = 0.5, label = "n.s. (p = 0.06)") +
  scale_fill_viridis_d() +
  scale_y_continuous(limits = c(0, 3.5e-6)) +
  labs(
    y = "Gaussian Kernel Density Estimate",
    fill = "Domain"
  ) +
  guides(fill = "none") +
  theme_bw() +
  theme(
    axis.title.x = element_blank()
  )
