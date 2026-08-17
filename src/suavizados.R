# =====================================================================
#  App Shiny: métodos de suavizamiento para series de tiempo
#  Media móvil unilateral | Suavización exponencial simple | Holt-Winters
#
#  Curso: Análisis de Series de Tiempo — UdeA
#
#  Cómo ejecutarla (desde PowerShell, parado en la raíz del repo):
#     Rscript.exe -e "shiny::runApp('src/suavizados.R', launch.browser=TRUE)"
#
#  Paquetes requeridos: shiny, forecast, ggplot2
#     install.packages(c("shiny","forecast","ggplot2"))
# =====================================================================

library(shiny)
library(forecast)
library(ggplot2)

# ---------------------------------------------------------------------
#  UTILIDADES
# ---------------------------------------------------------------------

# Convierte un objeto ts en data.frame para graficar con ggplot
ts_df <- function(x, serie) {
  data.frame(
    t     = as.numeric(stats::time(x)),
    v     = as.numeric(x),
    serie = serie,
    stringsAsFactors = FALSE
  )
}

# Interpreta la primera columna del CSV como eje temporal.
# Acepta tiempo numérico (1960, 1960.25) o fechas tipo "1960-03-31", "31/03/1960", etc.
parse_tiempo <- function(col) {
  if (is.numeric(col)) {
    return(list(t = col, tipo = "numérico"))
  }
  txt <- as.character(col)

  # Números guardados como texto: "1960", "1960.25"
  num <- suppressWarnings(as.numeric(txt))
  if (sum(!is.na(num)) > length(num) / 2) {
    return(list(t = num, tipo = "numérico"))
  }

  # Fechas: intentar múltiples formatos (ymd, dmy, mdy) para soportar regiones
  # Formato ISO: YYYY-MM-DD
  d <- tryCatch(suppressWarnings(as.Date(txt)),
                error = function(e) as.Date(rep(NA_character_, length(txt))))
  if (sum(!is.na(d)) > length(d) / 2) {
    anio <- as.numeric(format(d, "%Y"))
    dia  <- as.numeric(format(d, "%j"))
    return(list(t = anio + (dia - 1) / 365.25, tipo = "fecha"))
  }

  # Formato regional latino/europeo: DD/MM/YYYY o DD-MM-YYYY
  d <- tryCatch(suppressWarnings(as.Date(txt, format = "%d/%m/%Y")),
                error = function(e) as.Date(rep(NA_character_, length(txt))))
  if (sum(!is.na(d)) > length(d) / 2) {
    anio <- as.numeric(format(d, "%Y"))
    dia  <- as.numeric(format(d, "%j"))
    return(list(t = anio + (dia - 1) / 365.25, tipo = "fecha"))
  }

  # Formato alternativo: DD-MM-YYYY
  d <- tryCatch(suppressWarnings(as.Date(txt, format = "%d-%m-%Y")),
                error = function(e) as.Date(rep(NA_character_, length(txt))))
  if (sum(!is.na(d)) > length(d) / 2) {
    anio <- as.numeric(format(d, "%Y"))
    dia  <- as.numeric(format(d, "%j"))
    return(list(t = anio + (dia - 1) / 365.25, tipo = "fecha"))
  }

  # Formato MM/DD/YYYY (US)
  d <- tryCatch(suppressWarnings(as.Date(txt, format = "%m/%d/%Y")),
                error = function(e) as.Date(rep(NA_character_, length(txt))))
  if (sum(!is.na(d)) > length(d) / 2) {
    anio <- as.numeric(format(d, "%Y"))
    dia  <- as.numeric(format(d, "%j"))
    return(list(t = anio + (dia - 1) / 365.25, tipo = "fecha"))
  }

  # Último recurso: tratarla como un simple índice 1, 2, 3, ...
  list(t = seq_along(col), tipo = "índice")
}

# Deduce la frecuencia a partir del salto mínimo entre tiempos consecutivos
detectar_freq <- function(t) {
  u <- sort(unique(t[!is.na(t)]))
  if (length(u) < 2) return(1)
  d <- round(diff(u), 4)
  d <- d[d > 0]
  if (!length(d)) return(1)
  f <- round(1 / min(d))
  if (!is.finite(f) || f < 1 || f > 366) 1 else f
}

# Arma el objeto ts respetando el ciclo en que arranca la serie
construir_ts <- function(t, v, freq) {
  t0  <- min(t, na.rm = TRUE)
  anio <- floor(t0 + 1e-8)
  ciclo <- round((t0 - anio) * freq) + 1
  if (ciclo < 1 || ciclo > freq) ciclo <- 1
  stats::ts(v, start = c(anio, ciclo), frequency = freq)
}

# Métricas de error dentro de muestra
metricas <- function(y, ajustado) {
  yy <- as.numeric(y)
  ff <- as.numeric(ajustado)
  e  <- yy - ff
  ok <- is.finite(e)
  if (!any(ok)) return(c(RMSE = NA, MAE = NA, MAPE = NA))
  rmse <- sqrt(mean(e[ok]^2))
  mae  <- mean(abs(e[ok]))
  okp  <- ok & is.finite(yy) & yy != 0
  mape <- if (any(okp)) mean(abs(e[okp] / yy[okp])) * 100 else NA
  c(RMSE = rmse, MAE = mae, MAPE = mape)
}

# ---------------------------------------------------------------------
#  MÉTODOS DE SUAVIZAMIENTO
#  Todos devuelven la misma estructura para que la app los trate igual.
# ---------------------------------------------------------------------

# --- 1. Media móvil unilateral (trailing) ---
#   Suavizado:  S_t = (y_t + y_{t-1} + ... + y_{t-k+1}) / k
#   Pronóstico: un paso adelante -> ŷ_{t+1} = S_t ; a h pasos se mantiene plano.
ajustar_ma <- function(y, k, h) {
  n  <- length(y)
  fq <- stats::frequency(y)
  if (k > n) stop("El orden k de la media móvil supera el número de observaciones.")

  suav <- stats::filter(y, rep(1 / k, k), method = "convolution", sides = 1)
  suav <- as.numeric(suav)

  # Alineación de pronóstico: lo ajustado en t es el suavizado de t-1
  ajustado <- stats::ts(c(NA, suav[-n]), start = stats::start(y), frequency = fq)

  punto <- suav[n]
  media <- stats::ts(rep(punto, h),
                     start = stats::tsp(y)[2] + 1 / fq, frequency = fq)

  resid <- as.numeric(y) - as.numeric(ajustado)
  s     <- stats::sd(resid, na.rm = TRUE)

  list(
    nombre   = paste0("Media móvil unilateral (k = ", k, ")"),
    ajustado = ajustado,
    media    = media,
    lower    = stats::ts(rep(punto - 1.96 * s, h), start = stats::start(media), frequency = fq),
    upper    = stats::ts(rep(punto + 1.96 * s, h), start = stats::start(media), frequency = fq),
    resid    = stats::ts(resid, start = stats::start(y), frequency = fq),
    resumen  = paste0(
      "MEDIA MÓVIL UNILATERAL\n",
      "----------------------------------------\n",
      "Orden k              : ", k, "\n",
      "Observaciones        : ", n, "\n",
      "Último suavizado S_n : ", round(punto, 4), "\n",
      "Desv. est. residuos  : ", round(s, 4), "\n\n",
      "Nota: el pronóstico es plano (se repite S_n). El intervalo es\n",
      "aproximado: ±1.96·s, sin crecimiento con el horizonte."
    ),
    nota = NULL
  )
}

# --- 2. Suavización exponencial simple (SES) ---
ajustar_ses <- function(y, alpha, optimizar, h) {
  fc <- if (optimizar) forecast::ses(y, h = h) else forecast::ses(y, h = h, alpha = alpha)
  list(
    nombre   = if (optimizar) "SES (α optimizado)" else paste0("SES (α = ", alpha, ")"),
    ajustado = stats::fitted(fc),
    media    = fc$mean,
    lower    = fc$lower[, "95%"],
    upper    = fc$upper[, "95%"],
    resid    = fc$residuals,
    resumen  = paste(utils::capture.output(summary(fc)), collapse = "\n"),
    nota     = NULL
  )
}

# --- 3. Holt-Winters ---
#   Con frecuencia > 1 usa el modelo estacional completo (nivel+tendencia+estación).
#   Con frecuencia = 1 no hay estacionalidad posible: cae a Holt (nivel+tendencia).
ajustar_hw <- function(y, h, estacional, optimizar, alpha, beta, gamma) {
  fq   <- stats::frequency(y)
  nota <- NULL

  if (fq <= 1 || length(y) < 2 * fq) {
    nota <- paste0(
      "Holt-Winters estacional necesita frecuencia > 1 y al menos dos ciclos ",
      "completos. Esta serie no cumple, así que se ajustó Holt (nivel + tendencia, ",
      "sin componente estacional)."
    )
    fc <- if (optimizar) forecast::holt(y, h = h)
          else forecast::holt(y, h = h, alpha = alpha, beta = beta)
    nombre <- "Holt (sin estacionalidad)"
  } else {
    ajuste <- function(tipo) {
      if (optimizar) forecast::hw(y, h = h, seasonal = tipo)
      else forecast::hw(y, h = h, seasonal = tipo, alpha = alpha, beta = beta, gamma = gamma)
    }
    fc <- tryCatch(ajuste(estacional), error = function(e) NULL)
    if (is.null(fc) && estacional == "multiplicative") {
      nota <- paste0("El modelo multiplicativo falló (suele pasar con valores ",
                     "cero o negativos). Se usó el aditivo.")
      estacional <- "additive"
      fc <- ajuste("additive")
    }
    if (is.null(fc)) stop("No fue posible ajustar Holt-Winters con estos datos.")
    nombre <- paste0("Holt-Winters ",
                     if (estacional == "additive") "aditivo" else "multiplicativo",
                     " (m = ", fq, ")")
  }

  list(
    nombre   = nombre,
    ajustado = stats::fitted(fc),
    media    = fc$mean,
    lower    = fc$lower[, "95%"],
    upper    = fc$upper[, "95%"],
    resid    = fc$residuals,
    resumen  = paste(utils::capture.output(summary(fc)), collapse = "\n"),
    nota     = nota
  )
}

# ---------------------------------------------------------------------
#  INTERFAZ DE USUARIO
# ---------------------------------------------------------------------

COLORES <- c("Observado"  = "#37474f",
             "Ajustado"   = "#1e88e5",
             "Pronóstico" = "#e53935")

ui <- fluidPage(
  titlePanel("Métodos de suavizamiento: media móvil, SES y Holt-Winters"),

  sidebarLayout(
    sidebarPanel(
      width = 4,

      h4("1. Datos"),
      radioButtons("fuente", NULL,
                   choices = c("Cargar un CSV"    = "csv",
                               "Datos de ejemplo" = "demo"),
                   selected = "csv"),

      conditionalPanel(
        "input.fuente == 'csv'",
        fileInput("file", "Archivo CSV",
                  accept = c("text/csv", "text/comma-separated-values", ".csv"),
                  buttonLabel = "Buscar...", placeholder = "Ningún archivo"),
        helpText("Dos columnas: 1ª tiempo (numérico o fecha), 2ª observaciones.")
      ),
      conditionalPanel(
        "input.fuente == 'demo'",
        selectInput("demo_set", "Serie de ejemplo",
                    choices = c("AirPassengers (mensual)"        = "air",
                                "Johnson & Johnson (trimestral)" = "jj"))
      ),

      checkboxInput("auto_freq", "Detectar la frecuencia automáticamente", TRUE),
      conditionalPanel(
        "!input.auto_freq",
        numericInput("freq_manual", "Frecuencia (períodos por año)",
                     value = 4, min = 1, max = 366, step = 1),
        helpText("1 = anual · 4 = trimestral · 12 = mensual")
      ),

      hr(),
      h4("2. Método"),
      radioButtons("metodo", NULL,
                   choices = c("Media móvil unilateral" = "ma",
                               "Exponencial simple (SES)" = "ses",
                               "Holt-Winters" = "hw"),
                   selected = "ses"),

      conditionalPanel(
        "input.metodo == 'ma'",
        sliderInput("k", "Orden de la media móvil (k):",
                    min = 2, max = 24, value = 4, step = 1)
      ),

      conditionalPanel(
        "input.metodo == 'ses'",
        checkboxInput("ses_opt", "Optimizar α automáticamente", FALSE),
        conditionalPanel(
          "!input.ses_opt",
          sliderInput("alpha", "Nivel — α:", min = 0.01, max = 0.99,
                      value = 0.2, step = 0.01)
        )
      ),

      conditionalPanel(
        "input.metodo == 'hw'",
        radioButtons("hw_tipo", "Componente estacional:",
                     choices = c("Aditivo" = "additive",
                                 "Multiplicativo" = "multiplicative"),
                     selected = "additive", inline = TRUE),
        checkboxInput("hw_opt", "Optimizar α, β y γ automáticamente", TRUE),
        conditionalPanel(
          "!input.hw_opt",
          sliderInput("hw_alpha", "Nivel — α:",       0.01, 0.99, 0.20, 0.01),
          sliderInput("hw_beta",  "Tendencia — β:",   0.01, 0.99, 0.10, 0.01),
          sliderInput("hw_gamma", "Estación — γ:",    0.01, 0.99, 0.10, 0.01)
        )
      ),

      hr(),
      h4("3. Pronóstico"),
      sliderInput("h", "Horizonte (períodos adelante):",
                  min = 1, max = 12, value = 8, step = 1),

      hr(),
      div(style = "font-size:12px; color:#546e7a;", uiOutput("infoSerie"))
    ),

    mainPanel(
      width = 8,
      uiOutput("aviso"),
      tabsetPanel(
        type = "tabs",

        tabPanel(
          "Descomposición",
          br(),
          radioButtons("dec_tipo", "Tipo de descomposición:",
                       choices = c("Aditiva"       = "additive",
                                   "Multiplicativa" = "multiplicative"),
                       selected = "multiplicative", inline = TRUE),
          helpText(
            "Aditiva: y = tendencia + estacionalidad + residuo. Supone que la ",
            "amplitud estacional es constante. Multiplicativa: y = tendencia × ",
            "estacionalidad × residuo, para cuando la estacionalidad crece con el ",
            "nivel de la serie. Compara ambas mirando el panel de residuo: si queda ",
            "un patrón visible en forma de abanico, el tipo elegido no es el adecuado."
          ),
          plotOutput("plotDescomp", height = "500px"),
          br(),
          h4("Índices estacionales estimados"),
          tableOutput("tablaEstacional")
        ),

        tabPanel(
          "Pronóstico",
          br(),
          uiOutput("banner"),
          plotOutput("plotMain", height = "420px"),
          br(),
          h4("Ajuste dentro de muestra"),
          tableOutput("tablaMetricas"),
          h4("Resumen del modelo"),
          verbatimTextOutput("resumen")
        ),

        tabPanel(
          "Residuales",
          br(),
          helpText("Si el método capturó la estructura de la serie, los residuales ",
                   "deben verse sin patrón, centrados en cero y con la ACF dentro ",
                   "de las bandas."),
          plotOutput("plotResid", height = "230px"),
          plotOutput("plotAcf",   height = "230px"),
          plotOutput("plotHist",  height = "230px")
        ),

        tabPanel(
          "Comparación",
          br(),
          helpText("Los tres métodos ajustados a la misma serie, con los parámetros ",
                   "actuales de la barra lateral. Menor error no siempre significa ",
                   "mejor: un k pequeño memoriza la serie sin explicarla."),
          plotOutput("plotComp", height = "420px"),
          h4("Errores dentro de muestra"),
          tableOutput("tablaComp")
        )
      )
    )
  )
)

# ---------------------------------------------------------------------
#  SERVIDOR
# ---------------------------------------------------------------------

server <- function(input, output, session) {

  # --- Serie de tiempo (reactiva) ---
  serie <- reactive({
    if (input$fuente == "demo") {
      if (input$demo_set == "jj") {
        if (!requireNamespace("astsa", quietly = TRUE)) {
          validate(need(FALSE, "El paquete 'astsa' no está instalado. Usa AirPassengers o instala astsa."))
        }
        y <- astsa::jj
      } else {
        y <- datasets::AirPassengers
      }
      if (!input$auto_freq) {
        y <- construir_ts(as.numeric(stats::time(y)), as.numeric(y), input$freq_manual)
      }
      return(y)
    }

    req(input$file)
    # Fallback: intentar read.csv (separador ,), luego read.csv2 (separador ;)
    df <- tryCatch(
      utils::read.csv(input$file$datapath, stringsAsFactors = FALSE),
      error = function(e) {
        utils::read.csv2(input$file$datapath, stringsAsFactors = FALSE)
      }
    )
    validate(need(ncol(df) >= 2,
                  "El archivo debe tener al menos dos columnas (tiempo y observaciones)."))

    tt <- parse_tiempo(df[[1]])
    # Normalizar decimales: convertir comas a puntos antes de pasar a numérico
    col2_txt <- as.character(df[[2]])
    col2_norm <- gsub(",", ".", col2_txt)
    vv <- suppressWarnings(as.numeric(col2_norm))
    validate(need(sum(!is.na(vv)) > 2,
                  "La segunda columna no parece numérica. Revisa separadores decimales."))

    freq <- if (input$auto_freq) detectar_freq(tt$t) else input$freq_manual
    construir_ts(tt$t, vv, freq)
  })

  # --- Ajuste según el método elegido ---
  ajustar <- function(y, metodo) {
    switch(
      metodo,
      ma  = ajustar_ma(y, input$k, input$h),
      ses = ajustar_ses(y, input$alpha, input$ses_opt, input$h),
      hw  = ajustar_hw(y, input$h, input$hw_tipo, input$hw_opt,
                       input$hw_alpha, input$hw_beta, input$hw_gamma)
    )
  }

  modelo <- reactive({
    y <- serie()
    ajustar(y, input$metodo)
  })

  # --- Información de la serie cargada ---
  output$infoSerie <- renderUI({
    y <- serie()
    tp <- stats::tsp(y)
    HTML(paste0(
      "<b>Serie cargada</b><br>",
      "Observaciones: ", length(y), "<br>",
      "Frecuencia: ", stats::frequency(y),
      " (", switch(as.character(stats::frequency(y)),
                   "1" = "anual", "4" = "trimestral", "12" = "mensual",
                   "personalizada"), ")<br>",
      "Rango: ", round(tp[1], 3), " – ", round(tp[2], 3)
    ))
  })

  # --- Aviso cuando el método se degrada solo ---
  output$aviso <- renderUI({
    m <- modelo()
    if (is.null(m$nota)) return(NULL)
    div(style = "background:#fff8e1; border-left:5px solid #ffa000; padding:10px 14px;
                 border-radius:6px; margin-bottom:14px; font-size:14px;",
        strong("Aviso: "), m$nota)
  })

  # --- Descomposición clásica de la serie ---
  descomposicion <- reactive({
    y  <- serie()
    fq <- stats::frequency(y)

    validate(need(
      fq > 1,
      paste0("La descomposición necesita una serie estacional (frecuencia > 1). ",
             "Esta serie tiene frecuencia 1, así que no hay componente estacional ",
             "que separar. Si la frecuencia se detectó mal, corrígela en la barra lateral.")
    ))
    validate(need(
      length(y) >= 2 * fq,
      paste0("Hacen falta al menos dos ciclos completos (", 2 * fq,
             " observaciones) para estimar la componente estacional.")
    ))
    if (input$dec_tipo == "multiplicative") {
      validate(need(
        all(as.numeric(y) > 0, na.rm = TRUE),
        "La descomposición multiplicativa exige valores estrictamente positivos. Usa la aditiva."
      ))
    }
    stats::decompose(y, type = input$dec_tipo)
  })

  # Etiquetas legibles para los períodos del ciclo
  etiquetas_ciclo <- function(fq) {
    if (fq == 4)  return(paste0("T", 1:4))
    if (fq == 12) return(c("Ene", "Feb", "Mar", "Abr", "May", "Jun",
                           "Jul", "Ago", "Sep", "Oct", "Nov", "Dic"))
    as.character(seq_len(fq))
  }

  output$plotDescomp <- renderPlot({
    d <- descomposicion()

    comps <- list("Serie observada" = d$x,
                  "Tendencia"       = d$trend,
                  "Estacionalidad"  = d$seasonal,
                  "Residuo"         = d$random)

    df <- do.call(rbind, Map(function(x, nm) ts_df(x, nm), comps, names(comps)))
    df$serie <- factor(df$serie, levels = names(comps))

    ggplot(df, aes(t, v)) +
      geom_line(colour = "#1e88e5", linewidth = 0.7, na.rm = TRUE) +
      facet_grid(serie ~ ., scales = "free_y", switch = "y") +
      labs(
        title = paste0("Descomposición ",
                       if (input$dec_tipo == "additive") "aditiva" else "multiplicativa",
                       " (m = ", stats::frequency(d$x), ")"),
        subtitle = if (input$dec_tipo == "additive")
          "y(t) = tendencia + estacionalidad + residuo"
        else
          "y(t) = tendencia × estacionalidad × residuo",
        x = "Tiempo", y = NULL
      ) +
      theme_minimal(base_size = 13) +
      theme(strip.placement = "outside",
            strip.text.y.left = element_text(angle = 90, face = "bold"),
            panel.spacing = grid::unit(0.8, "lines"))
  })

  output$tablaEstacional <- renderTable({
    d  <- descomposicion()
    fq <- stats::frequency(d$x)
    data.frame(
      Período = etiquetas_ciclo(fq),
      Índice  = round(as.numeric(d$figure), 4),
      check.names = FALSE
    )
  }, striped = TRUE, width = "420px")

  # --- Banner con el pronóstico del próximo período ---
  output$banner <- renderUI({
    m  <- modelo()
    tn <- stats::time(m$media)[1]
    v  <- round(as.numeric(m$media)[1], 4)
    lo <- round(as.numeric(m$lower)[1], 4)
    hi <- round(as.numeric(m$upper)[1], 4)

    div(
      style = "background: linear-gradient(135deg,#e3f2fd 0%,#bbdefb 100%);
               border-left:6px solid #1e88e5; padding:20px; border-radius:8px;
               margin-bottom:20px; box-shadow:0 4px 6px rgba(0,0,0,.1);",
      div(style = "display:flex; justify-content:space-between; align-items:center;",
        div(
          h4(style = "margin:0; color:#0d47a1; font-weight:bold;",
             paste0("Pronóstico próximo período (t = ", round(tn, 3), ")")),
          h1(style = "margin:8px 0 0 0; color:#1565c0; font-size:42px; font-weight:bold;", v),
          p(style = "margin:4px 0 0 0; color:#546e7a; font-size:13px;", m$nombre)
        ),
        div(style = "text-align:right; color:#37474f; background:rgba(255,255,255,.7);
                     padding:10px 15px; border-radius:6px;",
          p(style = "margin:0; font-size:13px; font-weight:bold;", "INTERVALO 95%"),
          p(style = "margin:3px 0 0 0; font-size:17px; font-weight:bold; color:#2e7d32;",
            paste0("[ ", lo, "  —  ", hi, " ]"))
        )
      )
    )
  })

  # --- Gráfica principal ---
  output$plotMain <- renderPlot({
    m <- modelo()
    y <- serie()

    d_obs <- ts_df(y, "Observado")
    d_aju <- ts_df(m$ajustado, "Ajustado")

    # Se engancha el último dato observado al pronóstico para que la línea no salte
    d_fc  <- rbind(
      data.frame(t = d_obs$t[nrow(d_obs)], v = d_obs$v[nrow(d_obs)], serie = "Pronóstico"),
      ts_df(m$media, "Pronóstico")
    )
    d_int <- data.frame(t  = as.numeric(stats::time(m$media)),
                        lo = as.numeric(m$lower),
                        hi = as.numeric(m$upper))

    ggplot() +
      geom_ribbon(data = d_int, aes(x = t, ymin = lo, ymax = hi),
                  fill = "#e53935", alpha = 0.15) +
      geom_line(data = d_obs, aes(t, v, colour = serie), linewidth = 0.7) +
      geom_line(data = d_aju, aes(t, v, colour = serie), linewidth = 0.8, na.rm = TRUE) +
      geom_line(data = d_fc,  aes(t, v, colour = serie), linewidth = 1.0) +
      scale_colour_manual(values = COLORES, name = NULL) +
      labs(title = m$nombre,
           subtitle = paste0("Pronóstico a ", input$h, " período(s), banda del 95%"),
           x = "Tiempo", y = "Observaciones") +
      theme_minimal(base_size = 13) +
      theme(legend.position = "bottom")
  })

  # --- Métricas del método activo ---
  output$tablaMetricas <- renderTable({
    m  <- modelo()
    mt <- metricas(serie(), m$ajustado)
    data.frame(Métrica = c("RMSE", "MAE", "MAPE (%)"),
               Valor   = round(as.numeric(mt), 4))
  }, striped = TRUE, width = "320px")

  output$resumen <- renderPrint({
    cat(modelo()$resumen)
  })

  # --- Diagnóstico de residuales ---
  output$plotResid <- renderPlot({
    m <- modelo()
    d <- ts_df(m$resid, "res")
    ggplot(d, aes(t, v)) +
      geom_hline(yintercept = 0, colour = "#e53935", linetype = "dashed") +
      geom_line(colour = "#1e88e5", na.rm = TRUE) +
      labs(title = "Residuales en el tiempo", x = "Tiempo", y = "Residual") +
      theme_minimal(base_size = 13)
  })

  output$plotAcf <- renderPlot({
    r <- modelo()$resid
    forecast::ggAcf(r) +
      labs(title = "ACF de los residuales") +
      theme_minimal(base_size = 13)
  })

  output$plotHist <- renderPlot({
    r <- as.numeric(modelo()$resid)
    r <- r[is.finite(r)]
    ggplot(data.frame(r = r), aes(r)) +
      geom_histogram(bins = 20, fill = "#1e88e5", colour = "white") +
      labs(title = "Distribución de los residuales", x = "Residual", y = "Frecuencia") +
      theme_minimal(base_size = 13)
  })

  # --- Comparación de los tres métodos ---
  comparacion <- reactive({
    y <- serie()
    res <- list()
    for (met in c("ma", "ses", "hw")) {
      res[[met]] <- tryCatch(ajustar(y, met), error = function(e) NULL)
    }
    res[!vapply(res, is.null, logical(1))]
  })

  output$plotComp <- renderPlot({
    y   <- serie()
    cmp <- comparacion()
    validate(need(length(cmp) > 0, "No se pudo ajustar ningún método a esta serie."))

    d_aju <- do.call(rbind, lapply(cmp, function(m) ts_df(m$ajustado, m$nombre)))
    d_obs <- ts_df(y, "Observado")

    ggplot() +
      geom_line(data = d_obs, aes(t, v), colour = "#37474f",
                linewidth = 0.9, alpha = 0.85) +
      geom_line(data = d_aju, aes(t, v, colour = serie),
                linewidth = 0.8, na.rm = TRUE) +
      scale_colour_brewer(palette = "Set1", name = NULL) +
      labs(title = "Valores ajustados de los tres métodos",
           subtitle = "La línea gris oscura es la serie observada",
           x = "Tiempo", y = "Observaciones") +
      theme_minimal(base_size = 13) +
      theme(legend.position = "bottom", legend.direction = "vertical")
  })

  output$tablaComp <- renderTable({
    y   <- serie()
    cmp <- comparacion()
    validate(need(length(cmp) > 0, "No se pudo ajustar ningún método."))

    do.call(rbind, lapply(cmp, function(m) {
      mt <- metricas(y, m$ajustado)
      data.frame(Método = m$nombre,
                 RMSE   = round(mt[["RMSE"]], 4),
                 MAE    = round(mt[["MAE"]],  4),
                 `MAPE (%)` = round(mt[["MAPE"]], 4),
                 check.names = FALSE)
    }))
  }, striped = TRUE, width = "100%")
}

shinyApp(ui = ui, server = server)
