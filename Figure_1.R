Qlibrary(ggplot2)
library(tidyr)
library(dplyr)

work_dir_path <- gsub('Script','Results',getwd())
resource_dir_path <- gsub('Script','Resource',getwd())

#data <- read.table(paste0(resource_dir_path,"/COVID-19_visit_calc_date.tsv"), header = TRUE, sep = "\t")
#data <- read.delim(paste0(resource_dir_path,"/COVID-19_visit_calc_date.tsv"),header = TRUE, sep = "\t")
data <- read.delim(paste0(resource_dir_path,"/Figure1B_clinical_fu_days/COVID-19_visit_calc_date_v1v4.tsv"),header = TRUE, sep = "\t")


# print(data_ids)
# quit()
data <- data %>%
  mutate(ID = ifelse(ID == "C19-C014", "C19-C014*", ID)) %>% 
  mutate(color = ifelse(Severity == "severe", "red", "darkgreen"))

data_ids <- rev(data$ID)

data_long <- data %>%
#  select(ID, v1, v2, v3, v4, color) %>%
  select(ID, Acute_phase, Recovery_phase, color) %>%
#  pivot_longer(cols = c(v1, v2, v3, v4), names_to = "variable", values_to = "days") %>%
  pivot_longer(cols = c(Acute_phase, Recovery_phase), names_to = "variable", values_to = "days") %>%
  mutate(ID = factor(ID, levels = data_ids))

# data_long <- data %>%
#   select(ID, v1, v2, v3, v4, color) %>%
#   pivot_longer(cols = c(v1, v2, v3, v4), names_to = "variable", values_to = "days")

# data_long <- data_long %>%
#   arrange(desc(factor(ID, levels = data_ids)), ID)

# print(data_long$ID)
# quit()
plot <- ggplot(data_long, aes(x = days, y = ID, color = color, group = ID)) +
  geom_line(color="black") +
  geom_point(size = 5) +
  labs(x = "Follow-up days", y = "COVID-19 patients", color = "Severity", fontsize = 20) +
  scale_x_continuous(breaks=seq(0, 35, 7))+
  theme_bw() + 
  theme(
    axis.text.x = element_text(size=15),
    axis.title.x = element_text(size=20),
    axis.text.y = element_text(size=13),
    axis.title.y = element_text(size=20),
    panel.grid.minor.x = element_blank(), 
    panel.grid.major.y = element_blank(),
    legend.position='none'
  )+
  scale_color_identity()

ggsave("covid_followup_points_v1v4_v2_250421.png", plot, width = 10, height = 15, units = "in")
