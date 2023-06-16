library(tidyverse)
library(here)
library(tidypaleo)
theme_set(theme_paleo(8))
library(rioja)
library(vegan)
library(ggpubr)
library(scales)


dataRaw <- read_csv(here::here("data/raw-data/2021-12-08_nalkane_fames_concentration.csv"))

prepare_fig2data <- function(file_path) {
  # Código de preparación de datos aquí (desde la lectura del archivo CSV hasta la creación de data_bywax4 y summary_data4)
  # Select the columns of interest
  data1 <- dataRaw[-64:-91, ]
  data <- data1[,-27:-29]
  columns_of_interest <- data[,-1:-9]

# Code for data

# Calculate the total abundance for each row in the subset of columns

  total_abundance <- columns_of_interest %>% mutate(Total = rowSums(across(everything())))

# Calculate the relative abundance percentages in the subset of columns

  porcentajes <- total_abundance %>%
    mutate_at(vars(-Total), ~ ifelse(Total != 0, (. / Total) * 100, 0)) %>%
    select(-Total) %>%
    ungroup()

  # Replace columns 9 to 20 in the 'data' tibble with the calculated percentages
  data <- data1 %>% select(-10:-26) %>% bind_cols(porcentajes)

  data_bywax <- data %>%
    dplyr::select(Belts,wax,C22:C32) %>%
    dplyr::mutate(Belts = fct_relevel(Belts,
                                      "Coastal",
                                      "Hyperarid Desert",
                                      "Prepuna",
                                      "Puna",
                                      "Steppe",
                                      "Surface Sediment")) %>%
    tidyr::pivot_longer(cols = C22:C32,
                        names_to = "waxes",
                        values_drop_na = TRUE) %>%
    dplyr::group_by(Belts = as.factor(Belts), waxes, wax) %>%
    dplyr::summarize_if(is.numeric, funs(median, mad), na.rm = TRUE) %>%
    ungroup()

  data_bywax2 <- data %>%
    dplyr::select(Belts,wax,abu_total) %>%
    mutate(Belts = fct_relevel(Belts,
                               "Coastal",
                               "Hyperarid Desert",
                               "Prepuna",
                               "Puna",
                               "Steppe",
                               "Surface Sediment")) %>%
    tidyr::gather(key = "index", value = "values", abu_total)

  data_bywax3 <- data %>%
    dplyr::select(Belts,wax,ACL) %>%
    mutate(Belts = fct_relevel(Belts,
                               "Coastal",
                               "Hyperarid Desert",
                               "Prepuna",
                               "Puna",
                               "Steppe",
                               "Surface Sediment")) %>%
    tidyr::gather(key = "index", value = "values", ACL)

  summary_data3 <- suppressWarnings(summarySE(data_bywax3,
                             measurevar = "values",
                             groupvars = c("Belts", "wax", "index"),
                             na.rm = TRUE))

  data_bywax4 <- data %>%
    dplyr::select(Belts,wax,CPI) %>%
    mutate(Belts = fct_relevel(Belts,
                               "Coastal",
                               "Hyperarid Desert",
                               "Prepuna",
                               "Puna",
                               "Steppe",
                               "Surface Sediment")) %>%
    tidyr::gather(key = "index", value = "values", CPI)

  summary_data4 = suppressWarnings(summarySE(data_bywax4, measurevar="values", groupvars=c("Belts","wax","index"),na.rm = T))

  # Escribir los dataframes a archivos de datos de R
  write.csv(data_bywax, file = "data/derived-data/data_bywax.csv",
            row.names = FALSE)
  write.csv(data_bywax2, file = "data/derived-data/data_bywax2.csv",
            row.names = FALSE)
  write.csv(summary_data3, file ="data/derived-data/summary_data3.csv",
                      row.names = FALSE)
  write.csv(summary_data4, file = "data/derived-data/summary_data4.csv",
            row.names = FALSE)

  # Cargar los archivos de datos
  data_bywax <- read.csv("data/derived-data/data_bywax.csv")
  data_bywax2 <- read.csv("data/derived-data/data_bywax2.csv")
  summary_data3 <- read.csv("data/derived-data/summary_data3.csv")
  summary_data4 <- read.csv("data/derived-data/summary_data4.csv")

  # Retornar una lista con los data frames
  data <- list(data_bywax = data_bywax,
       data_bywax2 = data_bywax2,
       summary_data3 = summary_data3,
       summary_data4 = summary_data4)

# Llamar a la función y asignar los data frames al ambiente de trabajo
list2env(data, envir = .GlobalEnv)

}

prepare_fig2data(dataRaw)

generate_plot2 <- function(data_bywax, data_bywax2, summary_data3, summary_data4) {
  # Código para generar gráficos aquí (desde la creación de base1 hasta el retorno de p3)
  wax_palette <- c("#DAA520", "#7B68EE", "#B22222", "#FFD700", "#FFA07A", "#8FBC8F", "#CD5C5C", "#708090", "#CD853F", "#C71585", "#9ACD32", "#6B8E23", "#9370DB", "#8B4513")
  # Code for Figure 2d
  base1 <- ggplot(data_bywax,
                  aes(x=waxes,
                      y=median,
                      group=interaction(waxes, wax),
                      fill=wax)) +
    geom_col(position = position_dodge2(preserve = "single"), color = "black", size = 0.2) +
    facet_wrap(~Belts, nrow =7, scale="free_x") +
    geom_errorbar(aes(x=waxes,
                      ymin=median-mad,
                      ymax=median+mad),
                  colour="black",
                  width = 0.3,
                  position = position_dodge(width = 0.9)) +
    scale_x_discrete(breaks=c("C22","C23","C24",
                              "C25","C26","C27",
                              "C28","C29","C30",
                              "C31","C32"),
                     labels = c("C22","C23","C24",
                                "C25","C26","C27",
                                "C28","C29","C30",
                                "C31","C32"),
                     expand=c(0,0.6)) +
    coord_cartesian(ylim = c(0,60))

  base1 <- base1 +
    labs(
      x = "Carbon-chain length",
      y = "Relative abundances (%)",
      fill = "Plant waxes") +
    scale_fill_manual(values = wax_palette) + # Agregar la paleta de colores aquí
    theme_tufte() +
    theme(
      strip.text.x = element_text(margin = margin(2, 0, 2, 0), size = 20),
      panel.border=element_blank(),
      legend.position = ("top"),
      plot.title = element_text(face = "bold", size = 20),
      legend.justification = c(0, 1),
      text = element_text(size = 20),
      axis.ticks = element_line(colour = "grey50", size = 0.15),
      axis.title = element_text(size = rel(0.85), face = "bold"),
      axis.text = element_text(size = rel(0.70), face = "bold"),
      axis.line = element_line(color = "black"),
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 18),
      axis.title.y = element_text(size = 22),
      axis.text.y = element_text(size = 15),
      # La légende
      legend.title = element_text(size = 20, face = "bold"),
      legend.text = element_text(size = 20, face = "bold"),
      legend.key = element_rect(fill = "transparent", colour = NA),
      legend.key.size = unit(0.8, "lines"),
      legend.background = element_rect(fill = "transparent", colour = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )

  #####
  #Code for Figure 2c
  baseA <- ggplot(data_bywax2, aes(x = Belts, y = values, color = wax, fill = wax)) +
    geom_boxplot(outlier.shape = NA, colour = "black") +
    coord_cartesian(ylim = c(0, 800)) +
    geom_jitter(width = 0.1, size = 3, stroke = 1, shape = 21, color = "black") +
    scale_fill_manual(values = wax_palette) +
    scale_color_manual(values = wax_palette) +
    scale_x_discrete(labels = c("Coastal","H. Desert","Prepuna","Puna","Steppe","S. Sediment")) +
    labs(x = "", y = expression("T. Abundances (μg/gdw)"), color = "Plant waxes") +
    theme_tufte() +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1, size = 18),
      axis.title.y = element_text(size = 23),
      axis.text.y = element_text(size = 18),
      panel.border = element_blank(),
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 20),
      legend.justification = c(0, 1),
      text = element_text(size = 20),
      axis.ticks = element_line(colour = "grey50", size = 0.15),
      axis.title = element_text(size = rel(1), face = "bold"),
      axis.text = element_text(size = rel(0.85), face = "bold"),
      axis.line = element_line(color = "black"),
      legend.title = element_text(size = rel(1), face = "bold"),
      legend.text = element_text(size = rel(0.85), face = "bold"),
      legend.key = element_rect(fill = "transparent", colour = NA),
      legend.key.size = unit(0.8, "lines"),
      legend.background = element_rect(fill = "transparent", colour = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )

  #####
  # Code for Figure 2b
  baseB <- ggplot(summary_data3, aes(x=Belts, y=values, color=wax)) +
    geom_errorbar(aes(ymin=values-sd,
                      ymax=values+sd),
                  width=.3, size=0.7, show.legend = FALSE) +
    geom_point(aes(fill=wax), shape=21, size=4, colour="black") +
    scale_fill_manual(values = wax_palette) +
    scale_color_manual(values = wax_palette)

  baseB <- baseB +
    labs(
      x = "",
      y = expression("ACL Index"),
      fill = "Plant waxes") +
    theme_tufte() +
    theme(
      strip.text.x = element_text(margin = margin(2, 0, 2, 0)),
      legend.position = ("top"),
      axis.text.x = element_blank(), # Eliminar etiquetas del eje x
      axis.ticks.x = element_blank(), # Eliminar marcas de graduación del eje x
      axis.title.y = element_text(size = 23),
      panel.border=element_blank(),
      plot.title = element_text(face = "bold", size = 20),
      legend.justification = c(0, 1),
      text = element_text(size = 20),
      axis.ticks = element_line(colour = "grey50", size = 0.15),
      axis.title = element_text(size = rel(0.85), face = "bold"),
      axis.text = element_text(size = rel(0.70), face = "bold"),
      axis.line = element_line(color = "black"),
      legend.title = element_text(size = 20, face = "bold"),
      legend.text = element_text(size = 20, face = "bold"),
      legend.key = element_rect(fill = "transparent", colour = NA),
      legend.key.size = unit(0.8, "lines"),
      legend.background = element_rect(fill = "transparent", colour = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )

  #####
  # Code for Figure 2a
  base4 <-ggplot(summary_data4, aes(x=Belts, y=values, color=wax)) +
    geom_errorbar(aes(ymin=values-sd,
                      ymax=values+sd),
                  width=.3, size=0.7) +
    geom_point(aes(fill=wax),shape=21, size=4, colour="black") +
    scale_color_grey(start=0.8, end=0.2)

  base4 <- ggplot(summary_data4, aes(x=Belts, y=values, fill=wax)) +
    geom_errorbar(aes(ymin=values-sd,
                      ymax=values+sd, color=wax),
                  width=.3, size=0.7) +
    geom_point(shape=21, size=4, stroke=1, color="black") +
    scale_fill_manual(values=wax_palette) +
    scale_color_manual(values=wax_palette)

  base4 <- base4 +
    labs(
      x = "",
      y = expression("CPI Index"),
      fill = "Plant waxes") +
    theme_tufte() +
    theme(
      axis.text.x = element_blank(),
      panel.border=element_blank(),
      axis.title.y = element_text(size = 23),
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 20),
      text = element_text(size = 20),
      axis.ticks = element_blank(),
      axis.title = element_text(size = rel(0.85), face = "bold"),
      axis.text = element_text(size = rel(0.70), face = "bold"),
      axis.line = element_line(color = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )

  # Set the size of the final figure
  options(repr.plot.width = 7, repr.plot.height = 15)
  # Replace the lines where p1, p2 and p3 are created with the following code

  p1 <- ggpubr::ggarrange(baseB, base4, baseA, ncol = 1, nrow = 3, labels = c("(b)", "(c)", "(d)"), heights = c(0.75, 0.6, 1.1))
  p2 <- ggpubr::ggarrange(base1, ncol = 1, nrow = 1, labels = c("(a)"))

  # Adjusts the relative width of the plots
  p3 <- ggpubr::ggarrange(p2, p1, ncol = 2, nrow = 1, widths = c(1.1, 1))

  # Displays the final graphic
  p3
  return(p3)
}

generate_plot2(data_bywax, data_bywax2, summary_data3, summary_data4)
