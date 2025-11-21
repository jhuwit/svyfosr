# install.packages("hexSticker")
library(hexSticker)
library(ggplot2) # often used to make plots for the sticker
library(dplyr)

# create data for beta(s) and its two CI bands
s <- seq(0, 1, length.out = 200)
beta_hat <- sin(2 * pi * s) + 0.3 * s
ci_point <- 0.2 * exp(-3 * abs(s - 0.5))
ci_joint <- 0.4 * exp(-3 * abs(s - 0.5))

df = tibble(
  s = s,
  beta_hat = beta_hat,
  lower_point = beta_hat - ci_point,
  upper_point = beta_hat + ci_point,
  lower_joint = beta_hat - ci_joint,
  upper_joint = beta_hat + ci_joint
)

# build plot
p = ggplot(df, aes(x = s)) +
  # joint CI (wide, lighter)
  geom_ribbon(aes(ymin = lower_joint, ymax = upper_joint),
              fill = "#FAD7A0", alpha = 0.3) +
  # pointwise CI (narrow, darker)
  geom_ribbon(aes(ymin = lower_point, ymax = upper_point),
              fill = "#F5B041", alpha = 0.5) +
  # main beta(s) curve
  geom_line(aes(y = beta_hat), color = "#E59866", linewidth = 2) +
  theme_void() +
  coord_cartesian(ylim = c(-1.5, 1.5)) +
  theme(legend.position = "none") +
  theme_void() +
  theme(
    plot.margin = margin(0, 0, 0, 0),
    panel.spacing = unit(0, "lines")
  ) +
  coord_cartesian(clip = "off", expand = FALSE, ylim = c(-1.5, 1.5), xlim = c(0, 1))
p
# h_fill = "#283747",
# h_color = "#E59866"

# hex sticker


filename = "man/figures/logo.png"
res = sticker(
  subplot = p,
  package = "svyfosr",
  p_size = 18,
  p_color = "white",
  s_x    = 1.0,
  s_y      = 0.8,    # a bit higher vertically
  s_width  = 1,
  s_height = 0.9,
  h_size = 3,
  h_fill = "#283747",     # dark navy background
  h_color = "#E59866",    # orange border
  filename = filename,
  # url = "github.com/lilykoff/svyfosr",
  u_color = "white",
  white_around_sticker = FALSE,
  dpi = 300,
  spotlight = TRUE,
  l_alpha = 0.5,
  p_fontface = "bold"
)
res@layers$annotate$aes_params$size = 50
ggsave(plot = res, filename=filename)
img = magick::image_read(filename)
img %>%
  magick::image_trim() %>%
  magick::image_write(filename)
