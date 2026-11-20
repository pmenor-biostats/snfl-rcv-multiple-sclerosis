# theme_pmenor.R ----------------------------------------------------------
# pMenor · Estadística clínica en R
# Tema de ggplot2, paleta y exportación para carruseles, web y GitHub.
#
# Uso:
#   source("theme_pmenor.R")
#   usar_pmenor()                       # tema y colores por defecto
#   p <- ggplot(datos, aes(x, y)) + geom_point()
#   guardar_pmenor(p, "grafico.svg")    # 888 × 600 px, listo para Figma
#
# Requisitos: ggplot2 (>= 3.4), svglite, systemfonts y las fuentes
# IBM Plex Sans e IBM Plex Mono instaladas en el ordenador (Google Fonts).
# ------------------------------------------------------------------------

library(ggplot2)

# Colores de marca --------------------------------------------------------
pmenor_col <- c(
  burdeos = "#7A1E2E",  # solo donde está la idea
  tinta   = "#1C1A22",
  gris    = "#5B5660",
  papel   = "#F7F4F0",
  linea   = "#D9D3CC"
)

# Paleta para grupos: tinta y burdeos primero; con más de 3 grupos,
# mejor facetas que colores.
pmenor_paleta <- unname(pmenor_col[c("tinta", "burdeos", "gris")])

# Conversión de px a unidades de ggplot -----------------------------------
# guardar_pmenor() exporta a 72 ppp, así que 1 pt = 1 px en el SVG.
# El grosor de línea de ggplot va en mm y R lo dibuja a 0,75 pt por unidad
# de lwd, así que: px / (.pt × 0,75)
px_linea <- function(px) px / (ggplot2::.pt * 0.75)

# Grosor en px pensado para el carrusel (base 26) que se escala solo
# si el tema activo tiene otro tamaño de texto (informes, README).
grosor_pmenor <- function(px) {
  base <- theme_get()$text$size
  if (is.null(base)) base <- 26
  px_linea(px * base / 26)
}

# Números con coma decimal, punto de millares y signo menos tipográfico
numeros_pmenor <- function(...) {
  scales::label_number(decimal.mark = ",", big.mark = ".", style_negative = "minus", ...)
}

# Fuentes ------------------------------------------------------------------
registrar_fuentes_pmenor <- function() {
  if (requireNamespace("systemfonts", quietly = TRUE)) {
    disponibles <- systemfonts::system_fonts()$family
    faltan <- setdiff(c("IBM Plex Sans", "IBM Plex Mono"), disponibles)
    if (length(faltan) > 0) {
      warning("Faltan fuentes: ", paste(faltan, collapse = ", "),
              ". Instálalas desde Google Fonts para que el gráfico salga con la tipografía de marca.",
              call. = FALSE)
    }
    try(systemfonts::register_variant(
      name = "IBM Plex Sans Medium", family = "IBM Plex Sans", weight = "medium"
    ), silent = TRUE)
  }
  invisible(TRUE)
}

# Tema ---------------------------------------------------------------------
#' @param base_size  Tamaño del texto en px del lienzo final. 26 para
#'                   carruseles (1080 px de ancho); 11-12 para informes.
#' @param fondo      "transparente" (por defecto, para colocar sobre la
#'                   diapositiva), "papel" o "tinta".
#' @param rejilla    "y", "x", "xy" o "ninguna".
theme_pmenor <- function(base_size = 26,
                         fondo = c("transparente", "papel", "tinta"),
                         rejilla = c("y", "x", "xy", "ninguna")) {
  fondo   <- match.arg(fondo)
  rejilla <- match.arg(rejilla)

  oscuro      <- fondo == "tinta"
  col_texto   <- if (oscuro) pmenor_col[["papel"]] else pmenor_col[["tinta"]]
  col_sec     <- if (oscuro) pmenor_col[["linea"]] else pmenor_col[["gris"]]
  col_rejilla <- if (oscuro) "#3A3742" else pmenor_col[["linea"]]
  col_fondo   <- switch(fondo,
                        transparente = NA,
                        papel = pmenor_col[["papel"]],
                        tinta = pmenor_col[["tinta"]])

  linea_rejilla <- element_line(colour = col_rejilla, linewidth = px_linea(base_size / 13))

  theme_minimal(base_size = base_size, base_family = "IBM Plex Mono") %+replace%
    theme(
      text  = element_text(family = "IBM Plex Mono", colour = col_texto, size = base_size),

      # Títulos (en los carruseles el titular va en Figma; esto es para informes y GitHub)
      plot.title    = element_text(family = "IBM Plex Sans Medium", colour = col_texto,
                                   size = rel(1.4), hjust = 0, margin = margin(b = base_size * 0.5)),
      plot.subtitle = element_text(family = "IBM Plex Sans", colour = col_sec,
                                   size = rel(1), hjust = 0, margin = margin(b = base_size * 0.8)),
      plot.caption  = element_text(colour = col_sec, size = rel(0.8), hjust = 0,
                                   margin = margin(t = base_size * 0.8)),
      plot.title.position   = "plot",
      plot.caption.position = "plot",

      # Ejes: cifras en Mono, en gris
      axis.text   = element_text(colour = col_sec, size = rel(0.9)),
      axis.text.x = element_text(colour = col_sec, size = rel(1), margin = margin(t = base_size * 0.35)),
      axis.text.y = element_text(colour = col_sec, size = rel(1), hjust = 1, margin = margin(r = base_size * 0.35)),
      axis.title  = element_text(colour = col_sec, size = rel(0.9)),
      axis.title.x = element_text(colour = col_sec, size = rel(1), margin = margin(t = base_size * 0.6)),
      axis.title.y = element_text(colour = col_sec, size = rel(1), angle = 90, margin = margin(r = base_size * 0.6)),
      axis.line.x = element_line(colour = col_sec, linewidth = px_linea(base_size / 13)),
      axis.ticks  = element_blank(),

      # Rejilla discreta
      panel.grid.major   = linea_rejilla,
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = if (rejilla %in% c("x", "xy")) linea_rejilla else element_blank(),
      panel.grid.major.y = if (rejilla %in% c("y", "xy")) linea_rejilla else element_blank(),

      # Fondo
      plot.background  = element_rect(fill = col_fondo, colour = NA),
      panel.background = element_blank(),

      # Leyenda arriba a la izquierda, sin título
      legend.position      = "top",
      legend.justification = "left",
      legend.title         = element_blank(),
      legend.text          = element_text(colour = col_texto, size = rel(0.9)),
      legend.key.size      = unit(base_size * 0.9, "pt"),
      legend.margin        = margin(0, 0, base_size * 0.3, 0),

      # Facetas
      strip.text = element_text(family = "IBM Plex Sans Medium", colour = col_texto,
                                size = rel(1), hjust = 0, margin = margin(b = base_size * 0.4)),

      plot.margin = margin(base_size * 0.4, base_size * 0.6, base_size * 0.4, base_size * 0.4)
    )
}

# Escalas ------------------------------------------------------------------
scale_colour_pmenor <- function(...) {
  scale_colour_manual(values = pmenor_paleta, ...)
}
scale_fill_pmenor <- function(...) {
  scale_fill_manual(values = pmenor_paleta, ...)
}
scale_color_pmenor <- scale_colour_pmenor

# Activar por defecto ------------------------------------------------------
#' Fija theme_pmenor() y los colores de las geometrías para la sesión.
usar_pmenor <- function(base_size = 26, ...) {
  registrar_fuentes_pmenor()
  theme_set(theme_pmenor(base_size = base_size, ...))
  tinta   <- pmenor_col[["tinta"]]
  burdeos <- pmenor_col[["burdeos"]]
  gris    <- pmenor_col[["gris"]]
  esc <- base_size / 26  # todo escala con el tamaño de texto
  update_geom_defaults("point",   list(colour = tinta, size = 4 * esc))
  update_geom_defaults("line",    list(colour = tinta, linewidth = px_linea(4 * esc)))
  update_geom_defaults("path",    list(colour = tinta, linewidth = px_linea(4 * esc)))
  update_geom_defaults("col",     list(fill = tinta))
  update_geom_defaults("bar",     list(fill = tinta))
  update_geom_defaults("boxplot", list(colour = tinta, fill = NA))
  update_geom_defaults("smooth",  list(colour = burdeos, fill = pmenor_col[["linea"]]))
  update_geom_defaults("hline",   list(colour = gris, linewidth = px_linea(3 * esc)))
  update_geom_defaults("vline",   list(colour = gris, linewidth = px_linea(3 * esc)))
  update_geom_defaults("text",    list(family = "IBM Plex Mono", colour = gris, size = 0.85 * base_size / ggplot2::.pt))
  invisible(TRUE)
}

# Exportar -----------------------------------------------------------------
#' svglite quita el relleno de rectángulos y líneas con una regla CSS que
#' Figma ignora al importar (y entonces los pinta en negro). Esta función
#' escribe esos valores en cada elemento, quita los rectángulos de fondo
#' y pone el tamaño en px.
limpiar_svg_figma <- function(archivo, ancho, alto) {
  svg <- readLines(archivo, encoding = "UTF-8", warn = FALSE)
  svg <- svg[!grepl("^<\\?xml", svg)]
  # Rectángulos de fondo (lienzo y panel): fuera
  svg <- svg[!grepl("^<rect (width='100%'|x='[-0-9.e]+' y='0.00' width='[0-9.]+' height='[0-9.]+' style='stroke-width: [0-9.]+; stroke: none;' ?/>)", svg)]
  # Todo elemento con estilo pero sin relleno explícito: sin relleno
  sin_fill <- grepl("style='", svg) & !grepl("fill:", svg) & !grepl("^<text", svg)
  svg[sin_fill] <- sub("style='", "style='fill: none; ", svg[sin_fill])
  # Extremos de línea que daba el CSS
  es_trazo <- grepl("^<(line|polyline|path|polygon)", svg) & !grepl("stroke-linecap", svg)
  svg[es_trazo] <- sub("style='", "style='stroke-linecap: round; stroke-linejoin: round; ", svg[es_trazo])
  # Tamaño en px en lugar de pt
  svg <- sub("width='[0-9.]+pt' height='[0-9.]+pt'",
             sprintf("width='%s' height='%s'", ancho, alto), svg)
  writeLines(svg, archivo, useBytes = TRUE)
  invisible(archivo)
}

#' Guarda en SVG a 72 ppp: 1 px del SVG = 1 px de la diapositiva.
#' Por defecto, 888 × 600 px: el ancho útil del carrusel (1080 - 2 × 96).
#' Con extensión .png guarda un PNG al doble de resolución (web, GitHub).
guardar_pmenor <- function(plot = last_plot(), archivo, ancho = 888, alto = 600) {
  ext <- tolower(tools::file_ext(archivo))
  if (ext == "svg") {
    svglite::svglite(archivo, width = ancho / 72, height = alto / 72, bg = "transparent")
    print(plot)
    grDevices::dev.off()
    limpiar_svg_figma(archivo, ancho, alto)
  } else if (ext == "png") {
    ggsave(archivo, plot, width = ancho / 72, height = alto / 72, dpi = 144, bg = "transparent")
  } else {
    stop("Usa .svg (carruseles) o .png (web, GitHub).", call. = FALSE)
  }
  invisible(archivo)
}
