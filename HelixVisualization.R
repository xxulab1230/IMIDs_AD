# HelixVisualization.R - 3D Gene Score Mapping
# Author: Assistant
# Date: 2023-10-15

# Spiral Coordinate Generator ---------------------------------------------
generate_spiral_coords <- function(data, n_turns = 5) {
  data %>%
    arrange(desc(total_score)) %>%
    mutate(
      rank = row_number(),
      angle = rank * (2*pi*n_turns)/nrow(.), # Spiral angle calculation
      radius = rank * 1.2,                   # Radius increment
      x = radius * cos(angle),
      y = radius * sin(angle)
    )
}

# Visualization Construction ----------------------------------------------
top_genes %>% 
  generate_spiral_coords(n_turns = 4) %>% # Generate 4-turn spiral
  
  ggplot(aes(x = x, y = y)) +
  # Background System
  geom_segment( # Spiral guidelines
    data = . %>% mutate(angle_group = cut(angle, 50)),
    aes(xend = x*0.95, yend = y*0.95, color = angle),
    linewidth = 0.6, 
    alpha = 0.3
  ) +
  
  # Geometric Elements
  geom_point( # Main data points
    aes(fill = total_score, size = total_score),
    shape = 21, 
    color = "white", 
    stroke = 0.8
  ) +
  geom_text_repel( # Intelligent labels
    aes(label = gene), 
    color = "gray10",
    size = 4.5,
    family = "roboto",
    segment.color = "gray80",
    box.padding = 0.5
  ) +
  
  # Scale System
  scale_fill_gradientn(
    colors = viridis::magma(10), 
    guide = guide_colorbar(
      title = "Z-Score",
      barwidth = unit(3, "cm"),
      title.position = "top"
    )
  ) +
  scale_size_continuous(range = c(3, 10)) +
  
  # Theme System
  theme_void() +
  theme(
    legend.position = "bottom",
    plot.background = element_rect(fill = "#F7F7F7", color = NA),
    plot.title = element_markdown(hjust = 0.5, size = 18)
  ) +
  
  # Annotation System
  labs(
    title = "**Genomic Landscape in 3D Helix Space**",
    caption = "Script: HelixVisualization.R | Spiral mapping represents multidimensional associations"
  )