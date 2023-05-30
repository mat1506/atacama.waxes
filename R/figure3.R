# Load the necessary packages at the beginning of the script

library(tidyverse)
library(here)
library(tidypaleo)
theme_set(theme_paleo(8))
library(rioja)
library(vegan)
library(ggpubr)
library(scales)

# Define your function to prepare the data in Figure 3

prepare_macro3data <- function(file_path) {
  #n-alkane stratiplot
  #select the columns of interest
  macro0 <- file_path  %>%                                    # Apply filter & is.na
    dplyr::select(median,8:23)

  macro <- mutate(macro0,
                  Ex = Baccharis_aff_tola+Brassicaceae_aff_mathewsia_nivea+Pappostipa_frigida+
                    Cristaria+Fabiana_sp+Chenopodiaceae+Malvaceae+Adesmia_sp,
                  Lc = Phacelia_cuminingii+Phacelia_pinnatifida+Junellia_bryoides+
                    Cistanthe_sp+Ephedra_breana+Gilia_sp+Haplopappus_sp+Cryptantha_sp,
                  Total= Phacelia_cuminingii+Phacelia_pinnatifida+Junellia_bryoides+
                    Cistanthe_sp+Ephedra_breana+Gilia_sp+Haplopappus_sp+Cryptantha_sp+
                    Baccharis_aff_tola+Brassicaceae_aff_mathewsia_nivea+Pappostipa_frigida+
                    Cristaria+Fabiana_sp+Chenopodiaceae+Malvaceae+Adesmia_sp)
  ###### Variable selects

  macro3 <- macro  %>%
    dplyr::select(median,2:17) %>%
    rename(age=median,
           `Phacelia cuminingii`=Phacelia_cuminingii,
           `Phacelia pinnatifida`=Phacelia_pinnatifida,
           `Junellia bryoides`=Junellia_bryoides,
           `Cistanthe sp.`=Cistanthe_sp,
           `Ephedra breana`=Ephedra_breana,
           `Gilia sp.`=Gilia_sp,
           `Haplopappus sp.`=Haplopappus_sp,
           `Cryptantha sp.`=Cryptantha_sp,
           `Baccharis aff. tola`=Baccharis_aff_tola,
           `Brassicaceae aff. mathewsia nivea`=Brassicaceae_aff_mathewsia_nivea,
           `Pappostipa frígida`=Pappostipa_frigida,
           `Cristaria`=Cristaria,
           `Fabiana sp.`=Fabiana_sp,
           `Chenopodiaceae`=Chenopodiaceae,
           `Malvaceae`=Malvaceae,
           `Adesmia sp.`=Adesmia_sp) %>%
    gather(key = macro, value = abn,2:17)  %>%
    dplyr::mutate(macro = fct_relevel(macro,"Baccharis aff. tola",
                                      "Brassicaceae aff. mathewsia nivea",
                                      "Pappostipa frígida",
                                      "Cristaria",
                                      "Fabiana sp.",
                                      "Chenopodiaceae",
                                      "Malvaceae",
                                      "Adesmia sp.",
                                      "Phacelia cuminingii",
                                      "Phacelia pinnatifida",
                                      "Junellia bryoides",
                                      "Cistanthe sp.",
                                      "Ephedra breana",
                                      "Gilia sp.",
                                      "Haplopappus sp.",
                                      "Cryptantha sp."))


  # Escribir los dataframes a archivos de datos de R
  readr::write_csv(macro3, file = "data/derived-data/macro3.csv")
  readr::write_csv(macro, file = "data/derived-data/macro.csv")
  # Cargar los archivos de datos
  macro3 <- readr::read_csv("data/derived-data/macro3.csv")
  macro <- readr::read_csv("data/derived-data/macro.csv")

  # Retornar una lista con los data frames
  data <- list(macro3 = macro3,macro = macro)
  # Llamar a la función y asignar los data frames al ambiente de trabajo
  list2env(data, envir = .GlobalEnv)

}

# Read the data
dataplot3 <- read_csv(here::here("data/raw-data/2021-12-07_nalkane_concentration_middens.csv"))

# Call the function with the data frame read

resultado <- prepare_macro3data(dataplot3)

print(resultado)


# Function to prepare zone data

prepare_zone_data <- function() {
  tibble(ymin = c(18000,13000,11700,3200,1400),
         ymax = c(14800,8600,11900,2100,700),
         xmin = -Inf, xmax = Inf)
}

# Function to generate the graph of macrofossil in middens

generate_alkane_plot <- function(macro3, zone_data) {
  macroplot2 <- ggplot(macro3, aes(x = abn, y = age)) +
    geom_col_segsh() +
    geom_lineh(linewidth = 0.6) +
    geom_point(shape=20, size=4, colour="black") +
    geom_rect(mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax),
              data = zone_data,
              alpha = 0.3,
              fill = "blue",inherit.aes = FALSE ) +
    scale_y_reverse(breaks = c(0,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,11000,12000,13000,14000,15000,16000,17000)) +
    facet_abundanceh(vars(macro), scales = "fixed") +
    scale_x_continuous(breaks = c(1, 3, 5)) +
    labs(x = "Relative abundance index (0-5)", y = "Age a cal BP") +
    theme(
      legend.justification = c(0, 1),
      text = element_text(size = 18),
      axis.ticks = element_line(colour = "grey50", linewidth = 0.15),
      axis.title = element_text(size = rel(0.85), face = "bold"),
      axis.text = element_text(size = rel(0.70), face = "bold"),
      axis.text.x.bottom = element_text(angle=0, hjust=1),
      axis.line = element_line(color = "black"),
      axis.title.y = element_text(size = 22),
      axis.title.x = element_text(size = 22),
      # La leyenda
      legend.title = element_text(size = rel(0.85), face = "bold"),
      legend.text = element_text(size = rel(0.70), face = "bold"),
      legend.key = element_rect(fill = "transparent", colour = NA),
      legend.key.size = unit(0.8, "lines"),
      legend.background = element_rect(fill = "transparent", colour = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
}

# Function to generate the coniss graph

generate_coniss_plot <- function(coniss2) {
  coniss_plot <- ggplot() +
    layer_dendrogram(coniss2, aes(y = age)) +
    scale_y_reverse() +
    facet_geochem_gridh(vars("CONISS")) +
    labs(x = NULL) +
    theme(
      legend.justification = c(0, 1),
      text = element_text(size = 18),
      axis.ticks = element_line(colour = "grey50", linewidth = 0.15),
      axis.title = element_text(size = rel(0.85), face = "bold"),
      axis.text = element_text(size = rel(0.70), face = "bold"),
      axis.text.x.bottom = element_text(angle=0, hjust=1),
      axis.line = element_line(color = "black"),
      # La leyenda
      legend.title = element_text(size = rel(0.85), face = "bold"),
      legend.text = element_text(size = rel(0.70), face = "bold"),
      legend.key = element_rect(fill = "transparent", colour = NA),
      legend.key.size = unit(0.8, "lines"),
      legend.background = element_rect(fill = "transparent", colour = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
}

# Function to generate the resumen de los macros graph

generate_isotope_plot <- function(macro1) {
  macroplot1 <- ggplot(macro1, aes(x = abn, y = age)) +
    geom_lineh() +
    geom_point() +
    scale_y_reverse() +
    scale_x_continuous(n.breaks = 3) +
    facet_geochem_gridh(vars(macro)) +
    labs(x = NULL,y = NULL) +
    theme(
      legend.justification = c(0, 1),
      text = element_text(size = 18),
      axis.ticks = element_line(colour = "grey50", size = 0.15),
      axis.title = element_text(size = rel(0.85), face = "bold"),
      axis.text = element_text(size = rel(0.70), face = "bold"),
      axis.text.x.bottom = element_text(angle=0, hjust=1),
      axis.line = element_line(color = "black"),
      # La leyenda
      legend.title = element_text(size = rel(0.85), face = "bold"),
      legend.text = element_text(size = rel(0.70), face = "bold"),
      legend.key = element_rect(fill = "transparent", colour = NA),
      legend.key.size = unit(0.8, "lines"),
      legend.background = element_rect(fill = "transparent", colour = NA),
      panel.background = element_rect(fill = "transparent", colour = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
}

generate_plot3 <- function(macro3) {
  # Preparar los datos de zona
  zone_data <- prepare_zone_data()

  # Preparar los datos de coniss
  coniss2 <- macro3 %>%
    na.omit %>%
    nested_data(qualifiers = c(age), key = macro, value = abn, trans = scale) %>%
    nested_chclust_coniss(n_groups=10)

  # Preparar los resumen de los macros
  macro1 <- macro  %>%
    dplyr::select(median,Ex,Lc,Total) %>%
    rename(age=median) %>%
    gather(key = macro, value = abn,Ex,Lc,Total)  %>%
    dplyr::mutate(macro = fct_relevel(macro))

  # Generar gráficas individuales
  alk3plot <- generate_alkane_plot(macro3, zone_data)
  coniss1_plot <- generate_coniss_plot(coniss2)
  isoplot1 <- generate_isotope_plot(macro1)

  # Combinar gráficas
  library(patchwork)
  p3 <- wrap_plots(
    alk3plot +
      theme(strip.background = element_blank(), strip.text.y = element_blank()),
    isoplot1 +
      theme(axis.text.y.left = element_blank(),
            axis.ticks.y.left = element_blank(),
            panel.background = element_rect(fill = "transparent", colour = NA),
            plot.background = element_rect(fill = "transparent", colour = NA)) +
      coniss1_plot +
      theme(axis.text.y.left = element_blank(),
            axis.ticks.y.left = element_blank(),
            ,
            panel.background = element_rect(fill = "transparent", colour = NA),
            plot.background = element_rect(fill = "transparent", colour = NA)) +
      labs(y = NULL),
    nrow = 1,
    widths = c(3, 1)
  )

  return(p3)
}

generate_plot3(macro3)

