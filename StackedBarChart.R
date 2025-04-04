# StackedBarChart.R - Multi-cohort Contribution Analysis
# Author: Assistant
# Date: 2023-10-15

# Data Preparation --------------------------------------------------------
multi_cohort_data <- bind_rows(
  rush_data %>% mutate(cohort = "RUSH Cohort"),
  mayo_data %>% mutate(cohort = "Mayo Clinic"),
  adni_data %>% mutate(cohort = "ADNI Consortium")
) %>%
  group_by(gene, cohort) %>% 
  summarise(
    score = mean(expression_level, na.rm = TRUE), # Example metric
    .groups = "drop"
  ) %>%
  group_by(gene) %>% 
  mutate(
    total = sum(score),
    proportion = score / total
  )

# Visualization Construction ----------------------------------------------
multi_cohort_data %>%
  ggplot(aes(x = reorder(gene, total), y = proportion)) +
  # Geometric Elements
  geom_col(
    aes(fill = cohort),
    width = 0.85, 
    color = "white", 
    linewidth = 0.4
  ) +
  geom_text( # Proportional labels
    aes(label = if_else(proportion > 0.15, 
                        scales::percent(proportion, 1), 
                        "")),
    position = position_stack(vjust = 0.5),
    color = "white", 
    family = "roboto",
    fontface = "bold"
  ) +
  
  # Scale System
  scale_fill_manual(
    values = c("#264653", "#2A9D8F", "#E9C46A"), # Muted professional palette
    guide = guide_legend(
      title = "Study Cohort",
      reverse = TRUE,
      keyheight = unit(0.4, "cm")
    )
  ) +
  scale_y_continuous(
    labels = scales::percent,
    expand = expansion(mult = c(0, 0.1))
  ) +
  
  # Coordinate System
  coord_flip() +
  
  # Theme System
  theme_minimal(base_family = "roboto") +
  theme(
    axis.title = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = c(0.85, 0.15),
    plot.title.position = "plot",
    plot.title = element_markdown(size = 16),
    plot.subtitle = element_text(color = "gray40")
  ) +
  
  # Annotation System
  labs(
    title = "**Cross-Cohort Contribution Analysis**",
    subtitle = "Proportional contribution of each cohort to gene-level associations",
    caption = "Script: StackedBarChart.R"
  )