# BubblePlot.R - Visualize Gene Association Scores
# Author: Assistant
# Date: 2023-10-15

# Load Required Packages
library(ggplot2)
library(ggtext)
library(showtext)
library(scales)

# Initialize Typography
font_add_google("Roboto Condensed", "roboto") # Modern sans-serif font
showtext_auto()

# Data Preparation -----------------------------------------------------------
# Ensure dataframe contains: gene, total_score
top_genes <- your_data %>% 
  arrange(desc(total_score)) %>% # Sort by descending score
  head(20) %>%                   # Select top 20 genes
  mutate(
    total_score = round(total_score, 1),
    category = case_when(        # Add functional categories
      gene %in% c("APOE","APP") ~ "Amyloid Pathway",
      gene %in% c("TREM2","CXCL1") ~ "Neuroinflammation",
      TRUE ~ "Other Pathways"
    )
  )

# Visualization Parameters --------------------------------------------------
CUSTOM_PALETTE <- c(
  "#FF6B6B", # Primary red
  "#4ECDC4", # Teal
  "#45B7D1"  # Cyan
) # Modern flat color scheme

# Visualization Construction ------------------------------------------------
ggplot(top_genes, aes(x = reorder(gene, total_score), 
                      y = total_score)) +
  # Geometric Elements
  geom_segment( # Baseline connectors
    aes(xend = gene, y = 0, yend = total_score),
    color = "gray80", 
    linewidth = 0.6
  ) +
  geom_point( # Main bubbles
    aes(size = total_score, fill = category),
    shape = 21, 
    color = "white", 
    stroke = 1.2, 
    alpha = 0.9
  ) +
  
  # Label System
  geom_text( # Gene labels
    aes(label = gene), 
    hjust = -0.3, 
    size = 4.5, 
    family = "roboto",
    fontface = "italic"
  ) +
  geom_text( # Score labels
    aes(label = total_score), 
    hjust = 0.3, 
    vjust = -1.5, 
    size = 4,
    color = "gray30"
  ) +
  
  # Scale System
  scale_fill_manual(
    values = CUSTOM_PALETTE,
    name = "Functional Category"
  ) + 
  scale_size_continuous(
    range = c(5, 15), 
    guide = guide_legend(
      title = "Score Intensity",
      override.aes = list(fill = "#45B7D1")
    )
  ) +
  
  # Coordinate System
  coord_flip(ylim = c(0, max(top_genes$total_score)*1.1)) +
  
  # Theme System
  theme_minimal(base_family = "roboto") +
  theme(
    axis.text.x = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = c(0.85, 0.2),
    plot.title = element_markdown(size = 22, margin = margin(b = 15)),
    plot.subtitle = element_markdown(color = "gray40", size = 12)
  ) +
  
  # Annotation System
  labs(
    title = "**Top Genetic Associations in Alzheimer's Disease**",
    subtitle = "Bubble size represents cumulative association scores across multi-omics datasets",
    caption = "Script: BubblePlot.R"
  )