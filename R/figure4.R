# Load the necessary packages at the beginning of the script

library(tidyverse)
library(here)
library(tidypaleo)
theme_set(theme_paleo(8))
library(rioja)
library(vegan)
library(ggpubr)
library(scales)

# Define your function to prepare the data in Figure 4

prepare_alk4data <- function(file_path) {
  #n-alkane stratiplot
  #select the columns of interest
  alk0 <- file_path  %>%
    dplyr::select(median,24:43)

  alk <- mutate(alk0,
                Wt = C21+C22+C23+C24+C25+C26+C27+C28+C29+C30+C31+C32+C33+C34+C35, #	Ficken et al., 2000, Sikes et al., 2009
                ACL=(23*C23+25*C25+27*C27+29*C29+31*C31+33*C33)/(C23+C25+C27+C29+C31+C33),
                CPI=0.5*((C23+C25+C27+C29+C31+C33)/(C24+C26+C28+C30+C32))+(C23+C25+C27+C29+C31+C33)/(C24+C26+C28+C30+C32))

  ###### Variable selects

  alk4 <- alk  %>%
    dplyr::select(median,C21,C23,C25,C27,C29,C31,C33,Wt,ACL,CPI) %>%
    dplyr::rename(age = median) %>%
    tidyr::gather(key = alk, value = abn,2:11) %>%
    dplyr::mutate(alk = forcats::fct_relevel(alk,"C21","C23","C25","C27","C29",
                                             "C31","C33","Wt","ACL","CPI"))


  # Escribir los dataframes a archivos de datos de R
  readr::write_csv(alk4, file = "data/derived-data/alk4.csv")
  readr::write_csv(alk, file = "data/derived-data/alk.csv")
  # Cargar los archivos de datos
  alk4 <- readr::read_csv("data/derived-data/alk4.csv")
  alk <- readr::read_csv("data/derived-data/alk.csv")

  # Retornar una lista con los data frames
  data <- list(alk4 = alk4,alk = alk)
  # Llamar a la función y asignar los data frames al ambiente de trabajo
  list2env(data, envir = .GlobalEnv)

}

# Read the data
dataplot4 <- read_csv(here::here("data/raw-data/2021-12-07_nalkane_concentration_middens.csv"))

# Call the function with the data frame read

resultado <- prepare_alk4data(dataplot4)

print(resultado)


# Function to prepare zone data

prepare_zone_data <- function() {
  tibble(ymin = c(18000,13000,11700,3200,1400),
         ymax = c(14800,8600,11900,2100,700),
         xmin = -Inf, xmax = Inf)
}

# Function to generate the graph of alkanes in middens

generate_alkane_plot <- function(alk4, zone_data) {
  alk4 <- alk4 %>%
    mutate(alk = factor(alk, levels = c("C21","C23","C25","C27","C29",
                                        "C31","C33","Wt","ACL","CPI")))

  alk4plot <- ggplot(alk4, aes(x = abn, y = age)) +
    geom_col_segsh(size = 0.6) +
    geom_lineh(linewidth = 0.6) +
    geom_point(shape=20, size=4, colour="black") +
    geom_rect(mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax),
              data = zone_data,
              alpha = 0.2,
              fill = "blue",
              inherit.aes = FALSE ) +
    scale_y_reverse(breaks = c(0,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,11000,12000,13000,14000,15000,16000,17000)) +
    scale_x_continuous(n.breaks = 4,
                       breaks = function(x) pretty(x, n = 2),
                       labels = label_number()) +
    facet_geochem_gridh(vars(alk), scales = "free", space = "fixed") +
    labs(x = "n-Alkane abundance (μg/gdw)", y = "Age a cal BP") +
    theme(
      legend.justification = c(0, 1),
      text = element_text(size = 18),
      axis.ticks = element_line(colour = "grey50", linewidth = 0.15),
      axis.title = element_text(size = rel(0.85), face = "bold"),
      axis.text = element_text(size = rel(0.70), face = "bold"),
      axis.text.x.bottom = element_text(angle=0, hjust=1, size = 12),
      axis.line = element_line(color = "black"),
      axis.title.y = element_text(size = 23),
      axis.text.y = element_text(size = 15),
      axis.title.x = element_text(size = 23),
      # La légende
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

generate_coniss_plot <- function(coniss4) {
  ggplot() +
    layer_dendrogram(coniss4, aes(y = age)) +
    geom_lineh(linewidth = 0.6) +
    scale_y_reverse() +
    scale_x_continuous(n.breaks = 4) +
    facet_geochem_gridh(vars("CONISS")) +
    labs(x = NULL) +
    theme(
      legend.justification = c(0, 1),
      text = element_text(size = 18),
      axis.ticks = element_line(colour = "grey50", linewidth = 0.15),
      axis.title = element_text(size = rel(0.85), face = "bold"),
      axis.text = element_text(size = rel(0.70), face = "bold"),
      axis.text.x.bottom = element_text(angle = 0, hjust = 1, size = rel(1)),
      axis.line = element_line(color = "black"),
      # La leyenda
      legend.title = element_text(size = rel(0.85), face = "bold"),
      legend.text = element_text(size = rel(0.70), face = "bold"),
      legend.key = element_rect(fill = "transparent", colour = NA),
      legend.key.size = unit(0.8, "lines"),
      legend.background = element_rect(fill = "transparent", colour = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
    )
}

# Function to generate the isotopic graph

generate_isotope_plot <- function(iso) {
  ggplot(na.omit(iso), aes(x = abn, y = age)) +
    geom_lineh(linewidth = 0.6) +
    geom_point(shape=20, size=4, colour="black") +
    scale_y_reverse() +
    scale_x_continuous(n.breaks = 2,
                       breaks = function(x) pretty(x, n = 2),
                       labels = label_number()) +
    facet_geochem_gridh(vars(alk)) +
    labs(x = NULL,y = NULL) +
    theme(
      legend.justification = c(0, 1),
      text = element_text(size = 18),
      axis.ticks = element_line(colour = "grey50", linewidth = 0.15),
      axis.title = element_text(size = rel(0.85), face = "bold"),
      axis.text = element_text(size = rel(0.70), face = "bold"),
      axis.text.x.bottom = element_text(angle=0, hjust=1, size = rel(1)),
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


generate_plot4 <- function(alk4) {
  # Preparar los datos de zona
  zone_data <- prepare_zone_data()

  # Preparar los datos de alkane
  alkconis <- alk  %>%
    dplyr::select(median,C21,C23,C25,C27,C29,C31,C33,CN,d13C,d15N,Wt,ACL,CPI) %>%
    dplyr::rename(age = median, "C/N" =CN) %>%
    tidyr::gather(key = alk, value = abn,2:14) %>%
    dplyr::mutate(alk = fct_relevel(alk,"C21","C23","C25","C27","C29",
                                    "C31","C33","Wt","ACL","CPI","C/N","d13C","d15N"))

  # Preparar los datos de coniss
  coniss4 <- alkconis %>%
    tidypaleo::nested_data(qualifiers = c(age), key = alk, value = abn, trans = scale) %>%
    tidypaleo::nested_chclust_coniss(n_groups=10)

  # Preparar los datos isotópicos
  iso <- alk  %>%
    dplyr::select(median,CN,d13C,d15N) %>%
    rename(age=median,"C/N" =CN) %>%
    gather(key = alk, value = abn,"C/N",d13C,d15N)  %>%
    dplyr::mutate(alk = fct_relevel(alk,"C/N","d13C","d15N"))

  # Generar gráficas individuales
  alk4plot <- generate_alkane_plot(alk4, zone_data)
  coniss1_plot <- generate_coniss_plot(coniss4)
  isoplot1 <- generate_isotope_plot(iso)

  # Combinar gráficas
  library(patchwork)
  p3 <- wrap_plots(
    alk4plot +
      theme(strip.background = element_blank(), strip.text.y = element_blank()),
    isoplot1 +
      theme(axis.text.y.left = element_blank(), axis.ticks.y.left = element_blank()) +
      coniss1_plot +
      theme(axis.text.y.left = element_blank(), axis.ticks.y.left = element_blank()) +
      labs(y = NULL),
    nrow = 1,
    widths = c(3, 1)
  )

  return(p3)
}

generate_plot4(alk4)

