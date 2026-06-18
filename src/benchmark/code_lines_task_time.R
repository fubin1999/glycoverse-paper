library(tidyverse)
library(ggimage)
library(rstatix)

data <- read_csv("data/benchmark/task_duration.csv") |>
  mutate(
    glycoverse = as.factor(glycoverse),
    task = factor(task, levels = 1:6)
  )

ggplot(data, (aes(code_lines, fct_rev(task)))) +
  geom_point(aes(color = glycoverse), size = 4, alpha = 0.6) +
  scale_color_manual(values = c("#7E9BBF", "#A57BAB")) +
  geom_vline(xintercept = 0) +
  labs(x = "Code Lines") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 12, color = "black"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
  )
ggpreview(width = 2.5, height = 3)
ggsave("results/figures/code_lines.pdf", width = 2.5, height = 3)


ggplot(data, (aes(time, fct_rev(task)))) +
  geom_point(aes(color = glycoverse), size = 4, alpha = 0.6) +
  scale_color_manual(values = c("#7E9BBF", "#A57BAB")) +
  geom_vline(xintercept = 0) +
  labs(x = "Duration (min)") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 12, color = "black"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
  )
ggpreview(width = 2.5, height = 3)
ggsave("results/figures/duration.pdf", width = 2.5, height = 3)

total <- data |>
  summarise(
    time = sum(time),
    code_lines = sum(code_lines),
    .by = c(glycoverse, person)
  ) |>
  mutate(task = "total") |>
  relocate(person, task, glycoverse, time, code_lines)

table <- data |>
  bind_rows(total) |>
  summarise(
    mean_time = mean(time),
    mean_code_lines = mean(code_lines),
    sd_time = sd(time),
    sd_code_lines = sd(code_lines),
    .by = c(task, glycoverse)
  )
write_csv(table, "results/data/time_code_summary.csv")
