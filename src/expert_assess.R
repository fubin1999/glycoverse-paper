library(tidyverse)
library(ggimage)
library(cowplot)

data <- tribble(
  ~expert, ~I, ~II, ~III, ~IV, ~V,
  "expert 1", 4, 4, 4, 3, 3,
  "expert 2", 2, 2, 1, 5, 3,
  "expert 3", 4, 4, 3, 2, 3,
) %>%
  pivot_longer(cols = -expert, names_to = "task", values_to = "score")

radar_plot_cartesian <- function(
  data,
  task,
  value,
  title = NULL,
  group = NULL,
  start_deg = 90,
  r_max = NULL,
  r_breaks = NULL,
  close = TRUE,
  show_points = TRUE,
  color = "black",
  fill_alpha = 0.15
) {
  task  <- rlang::enquo(task)
  value <- rlang::enquo(value)
  group <- rlang::enquo(group)

  df <- data %>%
    mutate(.task = fct_inorder(!!task),
           .val  = !!value)

  n_axes <- dplyr::n_distinct(df$.task)
  start_rad <- start_deg/180*pi

  df <- df %>%
    group_by(!!group, .add = TRUE) %>%
    arrange(.task, .by_group = TRUE) %>%
    mutate(.theta = (as.numeric(.task) - 1) / n_axes * 2*pi + start_rad,
           x = .val * sin(.theta),
           y = .val * cos(.theta)) %>%
    ungroup()

  if (is.null(r_max))  r_max <- ceiling(max(df$.val, na.rm = TRUE))
  if (is.null(r_breaks)) r_breaks <- pretty(c(0, r_max), n = 5)

  df_poly <- df %>%
    group_by(!!group, .add = TRUE) %>%
    { if (close) bind_rows(., slice(., 1)) else . } %>%
    ungroup()

  circle_df <- function(r, n = 360) {
    t <- seq(0, 2*pi, length.out = n)
    tibble(r = r, t = t, x = r*sin(t), y = r*cos(t))
  }
  grid_circles <- purrr::map_dfr(r_breaks, circle_df, .id = "gid")

  theta_all <- df %>% distinct(.task, .theta)
  axis_len  <- 1.05 * r_max
  label_len <- 1.12 * r_max

  p <- ggplot() +
    geom_path(data = grid_circles, aes(x, y, group = gid), color = "grey85") +
    geom_segment(data = theta_all,
                 aes(x = 0, y = 0,
                     xend = axis_len * sin(.theta),
                     yend = axis_len * cos(.theta)),
                 color = "grey85") +
    geom_text(data = theta_all,
              aes(x = label_len * sin(.theta),
                  y = label_len * cos(.theta),
                  label = as.character(.task)),
              color = "grey30", size = 4) +
    coord_equal(xlim = c(-1.3*r_max, 1.3*r_max),
                ylim = c(-1.3*r_max, 1.3*r_max)) +
    theme_void()

  if (rlang::quo_is_null(group)) {
    p <- p +
      geom_polygon(data = df_poly, aes(x, y), fill = color, color = color,alpha = fill_alpha, linewidth = 0.7)
    if (show_points) p <- p + geom_point(data = df, aes(x, y), size = 1.8, color = color)
  } else {
    aes_grp <- aes(x, y, group = !!group, color = !!group, fill = !!group)
    p <- p +
      geom_polygon(data = df_poly, mapping = aes_grp, alpha = 0.15, linewidth = 0.7) +
      geom_path(data = df_poly,    mapping = aes_grp, linewidth = 0.7)
    if (show_points) p <- p + geom_point(data = df, mapping = aes_grp, size = 1.8)
    p <- p + guides(fill = "none") + theme(legend.position = "right")
  }

  if (!is.null(title)) {
    p <- p +
      labs(title = title) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
      )
  }

  p
}

p_expert1 <- data %>%
  filter(expert == "expert 1") %>%
  radar_plot_cartesian(task, score, title = "expert 1", color = "#A57BAB", fill_alpha = 0.15)

p_expert2 <- data %>%
  filter(expert == "expert 2") %>%
  radar_plot_cartesian(task, score, title = "expert 2", color = "#A57BAB", fill_alpha = 0.15)

p_expert3 <- data %>%
  filter(expert == "expert 3") %>%
  radar_plot_cartesian(task, score, title = "expert 3", color = "#A57BAB", fill_alpha = 0.15)

plot_grid(p_expert1, p_expert2, p_expert3, nrow = 1)
ggsave("results/figures/expert_assess.pdf", width = 7.5, height = 2.5)
