library(tidyverse)
library(ggimage)

data <- tribble(
  ~expert, ~task, ~code_lines, ~duration,
  "Expert 1", "Task 1", 85, 240,
  "Expert 1", "Task 2", 61, 60,
  "Expert 1", "Task 3", 77, 15,
  "Expert 1", "Task 4", 52, 10,
  "Expert 1", "Task 5", 128, 60,
  "Expert 1", "Task 6", 138, 60,
  "Expert 2", "Task 1", 116, 360,
  "Expert 2", "Task 2", 69, 60,
  "Expert 2", "Task 3", 75, 15,
  "Expert 2", "Task 4", 50, 15,
  "Expert 2", "Task 5", 66, 15,
  "Expert 2", "Task 6", 90, 30,
  "Expert 3", "Task 1", 141, 180,
  "Expert 3", "Task 2", 82, 30,
  "Expert 3", "Task 3", 65, 10,
  "Expert 3", "Task 4", 45, 15,
  "Expert 3", "Task 5", 80, 30,
  "Expert 3", "Task 6", 68, 60,
  "GPT-5", "Task 1", 173, NA,
  "GPT-5", "Task 2", 83, NA,
  "GPT-5", "Task 3", 62, NA,
  "GPT-5", "Task 4", 47, NA,
  "GPT-5", "Task 5", 99, NA,
  "GPT-5", "Task 6", 81, NA,
  "Claude-4-Sonnet", "Task 1", 120, NA,
  "Claude-4-Sonnet", "Task 2", 125, NA,
  "Claude-4-Sonnet", "Task 3", 125, NA,
  "Claude-4-Sonnet", "Task 4", 117, NA,
  "Claude-4-Sonnet", "Task 5", 207, NA,
  "Claude-4-Sonnet", "Task 6", 236, NA,
  "Gemini-2.5-Pro", "Task 1", 170, NA,
  "Gemini-2.5-Pro", "Task 2", 62, NA,
  "Gemini-2.5-Pro", "Task 3", 64, NA,
  "Gemini-2.5-Pro", "Task 4", 42, NA,
  "Gemini-2.5-Pro", "Task 5", 104, NA,
  "Gemini-2.5-Pro", "Task 6", 88, NA,
) |>
  mutate(expert_type = if_else(str_starts(expert, "Expert"), "Human", "AI"))


ggplot(data, (aes(code_lines, fct_rev(task)))) +
  geom_point(aes(color = expert_type), size = 4, alpha = 0.6) +
  scale_color_manual(values = c(Human = "#7E9BBF", AI = "#A57BAB")) +
  geom_vline(xintercept = 0) +
  labs(x = "Code Lines") +
  scale_x_continuous(limits = c(0, 250)) +
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


ggplot(data, (aes(duration, fct_rev(task)))) +
  geom_point(aes(color = expert_type), size = 4, alpha = 0.6) +
  scale_color_manual(values = c(Human = "#7E9BBF", AI = "#A57BAB")) +
  geom_vline(xintercept = 0) +
  labs(x = "Duration (min)") +
  scale_x_continuous(limits = c(0, 360)) +
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
