# =====================================================================
#  App Shiny: modelación SARIMA de la serie co2 (Mauna Loa, 1959-1997)
#
#  Curso: Series de Tiempo I, Universidad de Antioquia, 2026-2
#  Autor: Juan Camilo Melo Ortíz
#
#  Sigue la metodología de los notebooks del curso (Amazon, Cmort y
#  Air Passengers): descriptivo, suavizamiento, varianza, estacionariedad,
#  identificación, estimación, diagnóstico residual, pronóstico y
#  capacidad de pronóstico, comparando varios modelos a la vez.
#
#  Cómo ejecutarla (desde PowerShell, parado en la raíz del repo):
#     Rscript.exe -e "shiny::runApp('aplicacion-co2/app', launch.browser=TRUE)"
#
#  Paquetes:
#     install.packages(c("shiny", "bslib", "forecast", "tseries", "nortest",
#                        "car", "ggplot2", "plotly", "DT", "gridExtra"))
# =====================================================================

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(forecast)
  library(tseries)
  library(nortest)
  library(car)
  library(ggplot2)
  library(plotly)
  library(DT)
})

options(warn = -1)

# ---------------------------------------------------------------------
#  CONSTANTES Y UTILIDADES
# ---------------------------------------------------------------------

MESES <- c("Ene", "Feb", "Mar", "Abr", "May", "Jun", "Jul", "Ago", "Sep", "Oct", "Nov", "Dic")

# Modelos predefinidos: los cuatro candidatos del notebook
PREDEFINIDOS <- list(
  M1 = list(order = c(0, 1, 1), seasonal = c(0, 1, 1)),
  M2 = list(order = c(1, 1, 1), seasonal = c(0, 1, 1)),
  M3 = list(order = c(1, 1, 1), seasonal = c(1, 1, 2)),
  M4 = list(order = c(1, 1, 1), seasonal = c(2, 1, 2))
)

ETIQUETAS_BASE <- c(
  M1     = "M1: SARIMA(0,1,1)(0,1,1)",
  M2     = "M2: SARIMA(1,1,1)(0,1,1)",
  M3     = "M3: SARIMA(1,1,1)(1,1,2)",
  M4     = "M4: SARIMA(1,1,1)(2,1,2)",
  CUSTOM = "SARIMA personalizado",
  AUTO   = "auto.arima",
  HW     = "Holt-Winters",
  ETS    = "ETS automático",
  SNAIVE = "Ingenuo estacional"
)

COLORES <- c("#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e",
             "#8c564b", "#e377c2", "#17becf", "#7f7f7f")

tema <- theme_minimal(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom")

nombre_periodo <- function(freq) {
  if (freq == 12) MESES else if (freq == 4) paste0("T", 1:4) else as.character(seq_len(freq))
}

etiqueta_tiempo <- function(y) {
  f <- frequency(y)
  paste(nombre_periodo(f)[cycle(y)], floor(round(time(y), 6)))
}

# Pruebas ADF y KPSS en una fila
fila_estac <- function(x, nombre) {
  x <- na.omit(x)
  if (length(x) < 20) return(NULL)
  adf  <- suppressWarnings(adf.test(x))
  kpss <- suppressWarnings(kpss.test(x, null = "Level"))
  data.frame(
    Serie = nombre,
    `ADF estadístico` = round(unname(adf$statistic), 3),
    `ADF p-valor` = round(adf$p.value, 3),
    `KPSS estadístico` = round(unname(kpss$statistic), 3),
    `KPSS p-valor` = round(kpss$p.value, 3),
    `Conclusión` = ifelse(adf$p.value < 0.05 & kpss$p.value > 0.05, "Estacionaria",
                 ifelse(adf$p.value >= 0.05 & kpss$p.value <= 0.05, "No estacionaria",
                        "Pruebas en conflicto")),
    check.names = FALSE
  )
}

# ACF y PACF en ggplot con bandas del 95 %
graf_acf_pacf <- function(x, lags, titulo = "", freq = 12) {
  x <- as.numeric(na.omit(x))
  ci <- qnorm(0.975) / sqrt(length(x))
  a <- acf(x, lag.max = lags, plot = FALSE)$acf[-1]
  p <- pacf(x, lag.max = lags, plot = FALSE)$acf[, 1, 1]
  df <- rbind(data.frame(Panel = "ACF", Rezago = 1:lags, Valor = a),
              data.frame(Panel = "PACF", Rezago = 1:lags, Valor = p))
  g <- ggplot(df, aes(Rezago, Valor)) +
    geom_segment(aes(xend = Rezago, yend = 0), color = "steelblue", linewidth = 1) +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = c(-ci, ci), linetype = "dashed", color = "red") +
    facet_wrap(~ Panel, ncol = 1) +
    labs(title = titulo, x = "Rezago", y = "") + tema
  if (freq > 1) g <- g + geom_vline(xintercept = seq(freq, lags, freq), linetype = "dotted", color = "grey50")
  g
}

# Media móvil: centrada (2xk si k es par) o unilateral
media_movil <- function(y, k, centrada = TRUE) {
  if (centrada) {
    w <- if (k %% 2 == 0) c(0.5, rep(1, k - 1), 0.5) / k else rep(1 / k, k)
    stats::filter(y, w, sides = 2)
  } else {
    stats::filter(y, rep(1 / k, k), sides = 1)
  }
}

# ---------------------------------------------------------------------
#  AJUSTE Y PRONÓSTICO DE LOS MÉTODOS
# ---------------------------------------------------------------------

ajustar_arima <- function(y, order, seasonal, lambda) {
  f <- frequency(y)
  if (f == 1) seasonal <- c(0, 0, 0)
  tryCatch(
    Arima(y, order = order, seasonal = list(order = seasonal, period = f),
          method = "CSS-ML", lambda = lambda, biasadj = !is.null(lambda)),
    error = function(e) tryCatch(
      Arima(y, order = order, seasonal = list(order = seasonal, period = f),
            method = "ML", lambda = lambda, biasadj = !is.null(lambda)),
      error = function(e2) Arima(y, order = order, seasonal = list(order = seasonal, period = f),
                                 method = "CSS", lambda = lambda, biasadj = !is.null(lambda))))
}

# Ajusta un método. Devuelve list(id, tipo, fit, nombre)
ajustar_metodo <- function(id, y, lambda = NULL, custom = NULL, hw_tipo = "additive") {
  f <- frequency(y)
  if (id %in% names(PREDEFINIDOS)) {
    e <- PREDEFINIDOS[[id]]
    return(list(id = id, tipo = "arima", fit = ajustar_arima(y, e$order, e$seasonal, lambda)))
  }
  if (id == "CUSTOM") {
    return(list(id = id, tipo = "arima", fit = ajustar_arima(y, custom$order, custom$seasonal, lambda)))
  }
  if (id == "AUTO") {
    return(list(id = id, tipo = "arima",
                fit = auto.arima(y, lambda = lambda, biasadj = !is.null(lambda))))
  }
  if (id == "HW") {
    yt <- if (is.null(lambda)) y else BoxCox(y, lambda)
    tipo_s <- if (f > 1) hw_tipo else NULL
    fit <- if (f > 1) {
      tryCatch(HoltWinters(yt, seasonal = tipo_s),
               error = function(e) HoltWinters(yt, seasonal = "additive"))
    } else HoltWinters(yt, gamma = FALSE)
    return(list(id = id, tipo = "hw", fit = fit, lambda = lambda, y = y))
  }
  if (id == "ETS") {
    return(list(id = id, tipo = "ets", fit = ets(y, lambda = lambda, biasadj = !is.null(lambda))))
  }
  if (id == "SNAIVE") {
    return(list(id = id, tipo = "snaive", y = y, lambda = lambda,
                fit = snaive(y, h = f, lambda = lambda)))
  }
  stop("Método desconocido")
}

pronosticar <- function(obj, h, level = 95) {
  if (obj$tipo %in% c("arima", "ets")) {
    return(forecast(obj$fit, h = h, level = level))
  }
  if (obj$tipo == "hw") {
    fc <- forecast(obj$fit, h = h, level = level)
    if (!is.null(obj$lambda)) {
      fc$mean  <- InvBoxCox(fc$mean, obj$lambda)
      fc$lower <- InvBoxCox(fc$lower, obj$lambda)
      fc$upper <- InvBoxCox(fc$upper, obj$lambda)
      fc$x     <- obj$y
    }
    return(fc)
  }
  if (obj$tipo == "snaive") {
    return(snaive(obj$y, h = h, level = level, lambda = obj$lambda))
  }
}

valores_ajustados <- function(obj, y) {
  if (obj$tipo == "arima") return(fitted(obj$fit))
  if (obj$tipo == "ets") return(fitted(obj$fit))
  if (obj$tipo == "hw") {
    fv <- obj$fit$fitted[, "xhat"]
    if (!is.null(obj$lambda)) fv <- InvBoxCox(fv, obj$lambda)
    out <- y * NA
    out[match(round(time(fv), 6), round(time(y), 6))] <- as.numeric(fv)
    return(out)
  }
  if (obj$tipo == "snaive") return(fitted(obj$fit))
}

# Residuales para diagnóstico (en la escala en que se estimó el modelo)
residuales_diag <- function(obj, y) {
  if (obj$tipo == "arima") {
    r <- residuals(obj$fit)
    a <- obj$fit$arma            # p, q, P, Q, S, d, D
    corte <- a[6] + a[7] * a[5]
    if (corte > 0) r <- r[-seq_len(corte)]
    return(ts(as.numeric(r), end = end(y), frequency = frequency(y)))
  }
  if (obj$tipo == "ets") return(residuals(obj$fit, type = "innovation"))
  if (obj$tipo == "hw") {
    fv <- obj$fit$fitted[, "xhat"]
    yt <- if (is.null(obj$lambda)) y else BoxCox(y, obj$lambda)
    return(window(yt, start = start(fv)) - fv)
  }
  if (obj$tipo == "snaive") return(na.omit(residuals(obj$fit)))
}

# Raíces de los polinomios AR y MA completos (parte regular por parte estacional)
raices_arima <- function(fit) {
  quitar_ceros <- function(v) { while (length(v) > 0 && abs(tail(v, 1)) < 1e-10) v <- head(v, -1); v }
  phi <- quitar_ceros(fit$model$phi)
  theta <- quitar_ceros(fit$model$theta)
  list(ar = if (length(phi)) polyroot(c(1, -phi)) else complex(0),
       ma = if (length(theta)) polyroot(c(1, theta)) else complex(0))
}

n_parametros <- function(obj) {
  if (obj$tipo == "arima") return(length(coef(obj$fit)))
  if (obj$tipo == "ets") return(length(obj$fit$par))
  if (obj$tipo == "hw") return(3)
  0
}

nombre_modelo <- function(obj) {
  if (obj$tipo == "arima") {
    a <- obj$fit$arma
    s <- if (a[5] > 1 && (a[3] + a[4] + a[7]) > 0) sprintf("(%d,%d,%d)[%d]", a[3], a[7], a[4], a[5]) else ""
    return(sprintf("SARIMA(%d,%d,%d)%s", a[1], a[6], a[2], s))
  }
  if (obj$tipo == "ets") return(obj$fit$method)
  if (obj$tipo == "hw") return("Holt-Winters")
  "Ingenuo estacional"
}

# Tabla de pruebas sobre residuales
pruebas_residuales <- function(r, k, freq, lags_lb = c(12, 24, 36)) {
  r <- as.numeric(na.omit(r))
  lags_lb <- lags_lb[lags_lb > k & lags_lb < length(r) / 2]
  lb <- lapply(lags_lb, function(h) {
    b <- Box.test(r, lag = h, type = "Ljung-Box", fitdf = k)
    data.frame(Prueba = paste0("Ljung-Box (rezago ", h, ")"), `Hipótesis_nula` = "No autocorrelación",
               `Estadístico` = round(unname(b$statistic), 3), p_valor = b$p.value)
  })
  grupos <- factor(ceiling(seq_along(r) / max(freq, 12)))
  lev <- leveneTest(r ~ grupos)
  sw <- if (length(r) <= 5000) shapiro.test(r) else NULL
  ks <- lillie.test(r)
  jb <- jarque.bera.test(r)
  ml <- Box.test(r^2, lag = min(12, floor(length(r) / 4)), type = "Ljung-Box")
  tt <- t.test(r, mu = 0)
  out <- rbind(
    data.frame(Prueba = "t de media cero", `Hipótesis_nula` = "Media = 0",
               `Estadístico` = round(unname(tt$statistic), 3), p_valor = tt$p.value),
    do.call(rbind, lb),
    data.frame(Prueba = "McLeod-Li (residuales²)", `Hipótesis_nula` = "Sin efecto ARCH",
               `Estadístico` = round(unname(ml$statistic), 3), p_valor = ml$p.value),
    data.frame(Prueba = "Levene (por bloques)", `Hipótesis_nula` = "Varianza constante",
               `Estadístico` = round(lev$`F value`[1], 3), p_valor = lev$`Pr(>F)`[1]),
    if (!is.null(sw)) data.frame(Prueba = "Shapiro-Wilk", `Hipótesis_nula` = "Normalidad",
                                 `Estadístico` = round(unname(sw$statistic), 4), p_valor = sw$p.value),
    data.frame(Prueba = "Lilliefors (K-S)", `Hipótesis_nula` = "Normalidad",
               `Estadístico` = round(unname(ks$statistic), 4), p_valor = ks$p.value),
    data.frame(Prueba = "Jarque-Bera", `Hipótesis_nula` = "Normalidad",
               `Estadístico` = round(unname(jb$statistic), 3), p_valor = jb$p.value)
  )
  out$Decision <- ifelse(out$p_valor < 0.05, "Se rechaza H0", "No se rechaza H0")
  out$p_valor <- signif(out$p_valor, 4)
  names(out) <- c("Prueba", "Hipótesis nula", "Estadístico", "p-valor", "Decisión al 5 %")
  out
}

metricas <- function(real, pron, lo = NULL, hi = NULL, escala = NULL) {
  e <- real - pron
  out <- data.frame(ME = mean(e), RMSE = sqrt(mean(e^2)), MAE = mean(abs(e)),
                    MAPE = mean(abs(e / real)) * 100)
  if (!is.null(escala)) out$MASE <- mean(abs(e)) / escala
  if (!is.null(lo)) out$`Cobertura (%)` <- mean(real >= lo & real <= hi) * 100
  out
}

# ---------------------------------------------------------------------
#  INTERFAZ
# ---------------------------------------------------------------------

# ggplotly con leyenda horizontal limpia (sin las entradas duplicadas de las bandas)
plotly_h <- function(g, tooltip = "all") {
  p <- ggplotly(g, tooltip = tooltip)
  vistos <- character(0)
  for (i in seq_along(p$x$data)) {
    tr <- p$x$data[[i]]
    nm <- tr$name
    es_banda <- !is.null(tr$fill) && tr$fill %in% c("toself", "tonexty")
    if (is.null(nm) || !nzchar(nm)) next
    nm <- sub("^\\((.*),[0-9]+\\)$", "\\1", nm)
    p$x$data[[i]]$name <- nm
    p$x$data[[i]]$legendgroup <- nm
    if (es_banda || nm %in% vistos) {
      p$x$data[[i]]$showlegend <- FALSE
    } else {
      p$x$data[[i]]$showlegend <- TRUE
      vistos <- c(vistos, nm)
    }
  }
  layout(p, legend = list(orientation = "h", x = 0, y = -0.25), margin = list(b = 90))
}

caja <- function(titulo, ...) card(full_screen = TRUE, card_header(titulo), card_body(...))

nota <- function(...) div(class = "text-muted small", style = "margin-top: 6px;", ...)

ui <- page_navbar(
  title = "CO₂ Mauna Loa",
  theme = bs_theme(version = 5, bootswatch = "flatly", primary = "#1b5e20"),
  fillable = FALSE,
  sidebar = sidebar(
    width = 300,
    h6("Datos"),
    radioButtons("fuente", NULL, c("Serie co2 de R" = "co2", "Subir un CSV" = "csv")),
    conditionalPanel(
      "input.fuente == 'csv'",
      fileInput("archivo", "Archivo CSV", accept = ".csv"),
      nota("Se usa la última columna numérica del archivo."),
      numericInput("csv_freq", "Frecuencia", 12, min = 1, max = 52),
      numericInput("csv_anio", "Año inicial", 2000),
      numericInput("csv_per", "Periodo inicial", 1, min = 1)
    ),
    hr(),
    h6("Transformación de Box-Cox"),
    selectInput("bc", NULL, c("Ninguna" = "no", "Logaritmo (λ = 0)" = "log",
                              "λ de Guerrero" = "guerrero", "λ manual" = "manual")),
    conditionalPanel("input.bc == 'manual'",
                     sliderInput("lambda_man", "λ", min = -2, max = 2, value = 0.5, step = 0.05)),
    nota("Se aplica a todos los modelos. Los pronósticos se devuelven a la escala original con corrección de sesgo."),
    hr(),
    h6("Modelos a comparar"),
    checkboxGroupInput("metodos", NULL, choices = setNames(names(ETIQUETAS_BASE), ETIQUETAS_BASE),
                       selected = c("M1", "M2", "M3", "M4", "HW", "ETS", "SNAIVE")),
    conditionalPanel(
      "input.metodos.indexOf('CUSTOM') > -1",
      h6("SARIMA personalizado (p,d,q)(P,D,Q)"),
      fluidRow(column(4, numericInput("cp", "p", 2, 0, 5)), column(4, numericInput("cd", "d", 1, 0, 2)),
               column(4, numericInput("cq", "q", 1, 0, 5))),
      fluidRow(column(4, numericInput("cP", "P", 0, 0, 3)), column(4, numericInput("cD", "D", 1, 0, 1)),
               column(4, numericInput("cQ", "Q", 1, 0, 3)))
    ),
    radioButtons("hw_tipo", "Estacionalidad de Holt-Winters",
                 c("Aditiva" = "additive", "Multiplicativa" = "multiplicative"), inline = TRUE)
  ),

  # ---------------- 1. SERIE ----------------
  nav_panel(
    "1. Serie",
    layout_columns(
      col_widths = c(3, 3, 3, 3),
      value_box("Observaciones", textOutput("vb_n"), theme = "primary"),
      value_box("Periodo", textOutput("vb_periodo"), theme = "secondary"),
      value_box("Rango", textOutput("vb_rango"), theme = "info"),
      value_box("Crecimiento anual medio", textOutput("vb_crec"), theme = "success")
    ),
    caja("Serie de tiempo", plotlyOutput("p_serie", height = "380px"),
         nota("La serie co2 son las concentraciones mensuales de CO₂ (ppm) en Mauna Loa, Hawái, entre 1959 y 1997.")),
    layout_columns(
      col_widths = c(6, 6),
      caja("Desviación de cada periodo respecto a la tendencia (2×S)", plotOutput("p_box_mes", height = "340px")),
      caja("Perfil estacional por año", plotlyOutput("p_perfil", height = "340px"))
    ),
    caja("Tabla Buys-Ballot (filas: años, columnas: periodos)", DTOutput("t_bb"))
  ),

  # ---------------- 2. SUAVIZAMIENTO ----------------
  nav_panel(
    "2. Suavizamiento",
    navset_card_tab(
      nav_panel(
        "Medias móviles",
        layout_columns(
          col_widths = c(3, 9),
          div(
            radioButtons("ma_tipo", "Tipo", c("Centrada" = "c", "Unilateral (hacia atrás)" = "u")),
            checkboxGroupInput("ma_k", "Órdenes a mostrar", choices = c(3, 5, 7, 12, 24),
                               selected = c(3, 12), inline = TRUE),
            sliderInput("ma_extra", "Orden adicional", 2, 60, 9),
            checkboxInput("ma_extra_on", "Incluir orden adicional", FALSE),
            sliderInput("ma_zoom", "Ventana de años a mostrar", 1959, 1998, c(1959, 1998), sep = ""),
            nota("Con S par, la media centrada de orden S es la 2×S: pesos 1/(2S) en los extremos y 1/S en el centro. ",
                 "Es la única de estas ventanas que borra por completo la estacionalidad.")
          ),
          div(plotlyOutput("p_ma", height = "430px"))
        ),
        layout_columns(
          col_widths = c(6, 6),
          caja("Crecimiento por periodo de la tendencia 2×S (diferencia a S pasos)", plotOutput("p_crec", height = "300px")),
          caja("Promedio del crecimiento anual por década", tableOutput("t_crec"))
        )
      ),
      nav_panel(
        "Descomposición",
        radioButtons("desc_tipo", NULL, c("Aditiva" = "additive", "Multiplicativa" = "multiplicative"), inline = TRUE),
        plotOutput("p_desc", height = "560px"),
        layout_columns(col_widths = c(6, 6),
                       caja("Índices estacionales", tableOutput("t_indices")),
                       caja("Fuerza de los componentes", tableOutput("t_fuerza"),
                            nota("Fuerza = 1 − Var(aleatoria) / Var(componente + aleatoria), como en Hyndman y Athanasopoulos. Cerca de 1 es un componente fuerte.")))
      ),
      nav_panel(
        "Holt-Winters",
        layout_columns(
          col_widths = c(3, 9),
          div(
            checkboxInput("hw_opt", "Optimizar α, β, γ", TRUE),
            conditionalPanel("!input.hw_opt",
                             sliderInput("hw_a", "α (nivel)", 0.01, 1, 0.5, 0.01),
                             sliderInput("hw_b", "β (tendencia)", 0, 1, 0.01, 0.01),
                             sliderInput("hw_g", "γ (estacionalidad)", 0, 1, 0.3, 0.01)),
            tableOutput("t_hw_par")
          ),
          div(plotlyOutput("p_hw", height = "420px"))
        )
      )
    )
  ),

  # ---------------- 3. ESTACIONARIEDAD ----------------
  nav_panel(
    "3. Estacionariedad",
    layout_columns(
      col_widths = c(6, 6),
      caja("Homogeneidad de varianzas (Levene)",
           sliderInput("lev_k", "Número de ventanas", 5, 60, 39),
           verbatimTextOutput("o_levene"),
           plotOutput("p_lev", height = "260px")),
      caja("Box-Cox",
           tableOutput("t_bc"),
           plotOutput("p_bc", height = "260px"),
           nota("Si la log-verosimilitud es casi plana, la transformación no cambia nada importante y conviene no transformar."))
    ),
    caja("Diferenciación",
         layout_columns(
           col_widths = c(3, 9),
           div(
             radioButtons("dif_d", "Diferencias ordinarias d", 0:2, selected = 1, inline = TRUE),
             radioButtons("dif_D", "Diferencias estacionales D", 0:1, selected = 1, inline = TRUE),
             sliderInput("dif_lags", "Rezagos ACF/PACF", 12, 120, 48, step = 6),
             tableOutput("t_ndiffs")
           ),
           div(plotOutput("p_dif", height = "260px"))
         ),
         DTOutput("t_estac"),
         nota("ADF: H0 = raíz unitaria (no estacionaria). KPSS: H0 = estacionaria. Las dos de acuerdo dan una conclusión firme."),
         plotOutput("p_dif_acf", height = "440px"))
  ),

  # ---------------- 4. IDENTIFICACIÓN ----------------
  nav_panel(
    "4. Identificación",
    layout_columns(
      col_widths = c(3, 9),
      card(card_header("Búsqueda por criterios de información"),
           fluidRow(column(6, numericInput("g_p", "p máx", 2, 0, 3)), column(6, numericInput("g_q", "q máx", 2, 0, 3))),
           fluidRow(column(6, numericInput("g_P", "P máx", 2, 0, 2)), column(6, numericInput("g_Q", "Q máx", 2, 0, 2))),
           fluidRow(column(6, numericInput("g_d", "d", 1, 0, 2)), column(6, numericInput("g_D", "D", 1, 0, 1))),
           radioButtons("g_crit", "Ordenar por", c("BIC", "AIC", "AICc"), inline = TRUE),
           actionButton("g_run", "Buscar modelos", class = "btn-primary"),
           nota("Los criterios se reportan como en el curso: −2·logL/T + penalización/T. La columna AIC (R) es el AIC sin dividir por T.")),
      div(
        caja("Modelos ordenados", DTOutput("t_grid")),
        caja("auto.arima()", actionButton("auto_run", "Ejecutar auto.arima"), verbatimTextOutput("o_auto"))
      )
    )
  ),

  # ---------------- 5. MODELOS ----------------
  nav_panel(
    "5. Estimación",
    caja("Criterios de información y ajuste de los modelos seleccionados", DTOutput("t_ic"),
         nota("El AIC de los ETS no es comparable con el de los SARIMA (verosimilitudes sobre datos distintos). Holt-Winters y el ingenuo no tienen verosimilitud.")),
    layout_columns(
      col_widths = c(7, 5),
      caja("Coeficientes estimados y significancia", DTOutput("t_coef")),
      caja("Raíces inversas (deben quedar dentro del círculo unitario)",
           plotOutput("p_raices", height = "420px"),
           nota("Raíces inversas de los polinomios AR y MA completos (incluida la parte estacional). Una raíz pegada al círculo indica un modelo casi no estacionario o casi no invertible; raíces AR y MA superpuestas indican que se cancelan."))
    )
  ),

  # ---------------- 6. DIAGNÓSTICO ----------------
  nav_panel(
    "6. Diagnóstico",
    layout_columns(
      col_widths = c(3, 9),
      div(
        uiOutput("ui_diag_mod"),
        sliderInput("diag_lags", "Rezagos", 12, 72, 36, step = 6),
        nota("Para los SARIMA se omiten los primeros d + D·S residuales, que corresponden a la inicialización de las diferencias.")
      ),
      plotOutput("p_diag", height = "720px")
    ),
    caja("Pruebas sobre los residuales del modelo seleccionado", DTOutput("t_diag")),
    caja("Resumen de p-valores para todos los modelos", DTOutput("t_diag_todos"),
         nota("En rojo los p-valores menores que 0.05."))
  ),

  # ---------------- 7. PRONÓSTICO ----------------
  nav_panel(
    "7. Pronóstico",
    layout_columns(
      col_widths = c(3, 9),
      div(
        sliderInput("fc_h", "Horizonte (periodos)", 1, 60, 24),
        sliderInput("fc_nivel", "Nivel de confianza (%)", 80, 99, 95),
        sliderInput("fc_hist", "Años de historia a mostrar", 1, 39, 8),
        checkboxInput("fc_bandas", "Mostrar bandas de todos los modelos", FALSE),
        uiOutput("ui_fc_mod"),
        downloadButton("fc_desc", "Descargar pronósticos (CSV)", class = "btn-primary")
      ),
      plotlyOutput("p_fc", height = "520px")
    ),
    caja("Tabla de pronósticos del modelo seleccionado", DTOutput("t_fc"))
  ),

  # ---------------- 8. CAPACIDAD DE PRONÓSTICO ----------------
  nav_panel(
    "8. Capacidad predictiva",
    navset_card_tab(
      nav_panel(
        "Dentro de muestra",
        DTOutput("t_in"),
        plotlyOutput("p_in", height = "420px"),
        nota("El ajuste dentro de muestra siempre favorece a los modelos con más parámetros; la comparación justa es fuera de muestra.")
      ),
      nav_panel(
        "Partición entrenamiento/prueba",
        layout_columns(
          col_widths = c(3, 9),
          div(sliderInput("pt_train", "Porcentaje de entrenamiento", 50, 95, 80, step = 1),
              textOutput("pt_info"),
              nota("Se reestima cada modelo con el entrenamiento y se pronostica todo el periodo de prueba de una sola vez.")),
          plotlyOutput("p_pt", height = "440px")
        ),
        DTOutput("t_pt")
      ),
      nav_panel(
        "Validación cruzada",
        layout_columns(
          col_widths = c(3, 9),
          div(sliderInput("cv_k", "Número de pliegues", 3, 20, 10),
              sliderInput("cv_min", "Entrenamiento inicial (%)", 30, 80, 50),
              radioButtons("cv_h", "Horizonte en cada pliegue",
                           c("Todo el bloque (como en el curso)" = "bloque", "Fijo" = "fijo")),
              conditionalPanel("input.cv_h == 'fijo'", sliderInput("cv_hfijo", "h", 1, 36, 12)),
              actionButton("cv_run", "Ejecutar validación cruzada", class = "btn-primary"),
              nota("Ventana expansiva: en cada pliegue el entrenamiento crece un bloque y se pronostica el bloque siguiente.")),
          div(plotOutput("p_cv", height = "380px"))
        ),
        layout_columns(col_widths = c(7, 5),
                       caja("Resumen por modelo", DTOutput("t_cv")),
                       caja("Distribución del RMSE", plotOutput("p_cv_box", height = "320px")))
      ),
      nav_panel(
        "Un paso adelante",
        layout_columns(
          col_widths = c(3, 9),
          div(sliderInput("up_train", "Porcentaje de entrenamiento", 50, 95, 80),
              nota("Se estiman los parámetros con el entrenamiento y, sin reestimar, se pronostica un paso adelante cada mes del periodo de prueba (así se usaría el modelo en la práctica).")),
          plotlyOutput("p_up", height = "400px")
        ),
        DTOutput("t_up")
      )
    )
  ),

  # ---------------- 9. METODOLOGÍA ----------------
  nav_panel(
    "Metodología",
    card(card_body(
      h4("Ruta del análisis"),
      tags$ol(
        tags$li(strong("Serie: "), "gráfico, tabla Buys-Ballot y perfiles estacionales para ver tendencia y estacionalidad."),
        tags$li(strong("Suavizamiento: "), "medias móviles (la 2×12 extrae la tendencia), descomposición clásica y Holt-Winters."),
        tags$li(strong("Varianza y estacionariedad: "), "Levene por ventanas y Box-Cox para la varianza; ADF y KPSS, ndiffs y nsdiffs, y la ACF/PACF para escoger d y D."),
        tags$li(strong("Identificación: "), "lectura de la ACF/PACF de la serie diferenciada, búsqueda por AIC/BIC y auto.arima()."),
        tags$li(strong("Estimación: "), "coeficientes con su significancia, criterios de información y raíces de los polinomios."),
        tags$li(strong("Diagnóstico residual: "), "media cero, Ljung-Box, McLeod-Li, Levene, Shapiro-Wilk, Lilliefors y Jarque-Bera."),
        tags$li(strong("Pronóstico: "), "trayectorias con intervalos de predicción para varios modelos a la vez."),
        tags$li(strong("Capacidad de pronóstico: "), "dentro de muestra, partición entrenamiento/prueba, validación cruzada con ventana expansiva y pronóstico a un paso.")
      ),
      h4("Resultado con la serie co2"),
      p("La serie necesita una diferencia ordinaria y una estacional, (1−B)(1−B¹²), y no necesita Box-Cox. ",
        "La ACF y la PACF de la serie diferenciada y el BIC apuntan al modelo de aerolíneas ",
        strong("SARIMA(0,1,1)(0,1,1)₁₂"), " (M1), con θ ≈ −0.35 y Θ ≈ −0.85, ambos muy significativos y con residuales de ruido blanco normal. ",
        "Los modelos más grandes (M2, M3 de auto.arima y M4 de mínimo AIC) agregan parámetros no significativos o polinomios estacionales que se contrarrestan, ",
        "y en la validación cruzada no pronostican mejor que M1. Por parsimonia se escoge M1."),
      p("El análisis completo, paso a paso, está en ", code("aplicacion-co2/notebooks/co2_sarima.ipynb"), "."),
      hr(),
      p(class = "text-muted", "Series de Tiempo I · Universidad de Antioquia · 2026-2 · Juan Camilo Melo Ortíz")
    ))
  )
)

# ---------------------------------------------------------------------
#  SERVIDOR
# ---------------------------------------------------------------------

server <- function(input, output, session) {

  # ---------- Datos ----------
  serie <- reactive({
    if (input$fuente == "co2") return(datasets::co2)
    req(input$archivo)
    d <- tryCatch(read.csv(input$archivo$datapath), error = function(e) NULL)
    validate(need(!is.null(d), "No se pudo leer el archivo."))
    num <- which(sapply(d, is.numeric))
    validate(need(length(num) > 0, "El CSV no tiene columnas numéricas."))
    v <- d[[tail(num, 1)]]
    validate(need(sum(!is.na(v)) >= 2 * input$csv_freq + 10, "La serie es muy corta para esta frecuencia."))
    ts(na.interp(ts(v, frequency = input$csv_freq)), start = c(input$csv_anio, input$csv_per),
       frequency = input$csv_freq)
  })

  freq <- reactive(frequency(serie()))

  lambda <- reactive({
    switch(input$bc,
           no = NULL, log = 0,
           guerrero = BoxCox.lambda(serie(), method = "guerrero"),
           manual = input$lambda_man)
  })

  observeEvent(serie(), {
    rango <- range(floor(time(serie())))
    updateSliderInput(session, "ma_zoom", min = rango[1], max = rango[2] + 1, value = c(rango[1], rango[2] + 1))
    updateSliderInput(session, "fc_hist", max = diff(rango) + 1, value = min(8, diff(rango) + 1))
    updateSliderInput(session, "lev_k", value = min(60, max(5, floor(length(serie()) / max(freq(), 4)))))
  })

  custom_spec <- reactive(list(order = c(input$cp, input$cd, input$cq),
                               seasonal = c(input$cP, input$cD, input$cQ)))

  etiquetas <- reactive({
    e <- ETIQUETAS_BASE
    e["CUSTOM"] <- sprintf("SARIMA(%d,%d,%d)(%d,%d,%d)", input$cp, input$cd, input$cq, input$cP, input$cD, input$cQ)
    e
  })

  # Ajuste de todos los métodos seleccionados sobre la serie completa
  ajustes <- reactive({
    req(length(input$metodos) > 0)
    y <- serie()
    withProgress(message = "Ajustando modelos", value = 0, {
      out <- list()
      for (m in input$metodos) {
        incProgress(1 / length(input$metodos), detail = etiquetas()[m])
        out[[m]] <- tryCatch(ajustar_metodo(m, y, lambda(), custom_spec(), input$hw_tipo),
                             error = function(e) NULL)
      }
      out[!sapply(out, is.null)]
    })
  })

  # Re-ajuste sobre una submuestra (para partición y validación cruzada)
  pronostico_submuestra <- function(m, y_tr, h, nivel = 95) {
    obj <- ajustar_metodo(m, y_tr, lambda(), custom_spec(), input$hw_tipo)
    pronosticar(obj, h, nivel)
  }

  colores_mod <- reactive(setNames(COLORES[seq_along(ETIQUETAS_BASE)], etiquetas()))

  # ---------- 1. Serie ----------
  output$vb_n <- renderText(length(serie()))
  output$vb_periodo <- renderText(paste(etiqueta_tiempo(serie())[1], "a", tail(etiqueta_tiempo(serie()), 1)))
  output$vb_rango <- renderText(paste(round(min(serie()), 2), "a", round(max(serie()), 2)))
  output$vb_crec <- renderText({
    y <- serie(); f <- freq()
    if (f == 1) return(round(mean(diff(y)), 3))
    paste(round(mean(diff(y, lag = f)), 3), "por año")
  })

  output$p_serie <- renderPlotly({
    y <- serie()
    df <- data.frame(Tiempo = as.numeric(time(y)), Valor = as.numeric(y), Etiqueta = etiqueta_tiempo(y))
    g <- ggplot(df, aes(Tiempo, Valor, text = paste(Etiqueta, "<br>", round(Valor, 2)), group = 1)) +
      geom_line(color = "darkgreen", linewidth = 0.6) +
      labs(x = "Tiempo", y = "Valor") + tema
    ggplotly(g, tooltip = "text")
  })

  output$p_box_mes <- renderPlot({
    y <- serie(); f <- freq()
    validate(need(f > 1, "La serie no tiene estacionalidad (frecuencia 1)."))
    tend <- media_movil(y, f, TRUE)
    df <- na.omit(data.frame(Periodo = factor(nombre_periodo(f)[cycle(y)], levels = nombre_periodo(f)),
                             Desv = as.numeric(y - tend)))
    ggplot(df, aes(Periodo, Desv)) + geom_boxplot(fill = "skyblue", color = "darkblue") +
      geom_hline(yintercept = 0, linetype = "dashed") + labs(x = "", y = "Desviación") + tema
  })

  output$p_perfil <- renderPlotly({
    y <- serie(); f <- freq()
    validate(need(f > 1, "La serie no tiene estacionalidad (frecuencia 1)."))
    df <- data.frame(Periodo = as.numeric(cycle(y)), Anio = floor(round(time(y), 6)), Valor = as.numeric(y))
    g <- ggplot(df, aes(Periodo, Valor, color = Anio, group = Anio,
                        text = paste(Anio, nombre_periodo(f)[Periodo], round(Valor, 2)))) +
      geom_line() + scale_color_viridis_c() +
      scale_x_continuous(breaks = 1:f, labels = nombre_periodo(f)) +
      labs(x = "", y = "Valor", color = "Año") + tema
    ggplotly(g, tooltip = "text")
  })

  output$t_bb <- renderDT({
    y <- serie(); f <- freq()
    validate(need(f > 1, "La tabla Buys-Ballot necesita frecuencia mayor que 1."))
    anios <- floor(round(time(y), 6))
    m <- tapply(as.numeric(y), list(anios, factor(nombre_periodo(f)[cycle(y)], levels = nombre_periodo(f))), mean)
    m <- cbind(m, Promedio = rowMeans(m, na.rm = TRUE))
    m <- rbind(m, Promedio = colMeans(m, na.rm = TRUE))
    datatable(round(m, 2), options = list(pageLength = 10, scrollX = TRUE, dom = "tip"))
  })

  # ---------- 2. Suavizamiento ----------
  output$p_ma <- renderPlotly({
    y <- serie()
    ks <- as.numeric(input$ma_k)
    if (input$ma_extra_on) ks <- c(ks, input$ma_extra)
    ks <- unique(ks)
    validate(need(length(ks) > 0, "Escoja al menos un orden."))
    t <- as.numeric(time(y))
    df <- do.call(rbind, lapply(ks, function(k) {
      nom <- if (input$ma_tipo == "c" && k %% 2 == 0) paste0("MA 2x", k) else paste0("MA(", k, ")")
      data.frame(Tiempo = t, Valor = as.numeric(media_movil(y, k, input$ma_tipo == "c")), Media = nom)
    }))
    zoom <- input$ma_zoom
    df <- subset(na.omit(df), Tiempo >= zoom[1] & Tiempo < zoom[2])
    dfo <- subset(data.frame(Tiempo = t, Valor = as.numeric(y)), Tiempo >= zoom[1] & Tiempo < zoom[2])
    g <- ggplot() +
      geom_line(data = dfo, aes(Tiempo, Valor), color = "grey60") +
      geom_line(data = df, aes(Tiempo, Valor, color = Media), linewidth = 0.9) +
      labs(title = "Serie original (gris) y medias móviles", x = "Tiempo", y = "Valor", color = "") + tema
    plotly_h(g)
  })

  output$p_crec <- renderPlot({
    y <- serie(); f <- freq()
    validate(need(f > 1, "Necesita frecuencia mayor que 1."))
    cr <- diff(media_movil(y, f, TRUE), lag = f)
    df <- na.omit(data.frame(Tiempo = as.numeric(time(cr)), Crec = as.numeric(cr)))
    ggplot(df, aes(Tiempo, Crec)) + geom_line(color = "firebrick") +
      geom_smooth(method = "loess", formula = y ~ x, se = FALSE, color = "black", linetype = "dashed") +
      labs(x = "Tiempo", y = "Crecimiento en S periodos") + tema
  })

  output$t_crec <- renderTable({
    y <- serie(); f <- freq()
    validate(need(f > 1, ""))
    cr <- diff(media_movil(y, f, TRUE), lag = f)
    df <- na.omit(data.frame(Decada = paste0(floor(time(cr) / 10) * 10, "s"), Crec = as.numeric(cr)))
    a <- aggregate(Crec ~ Decada, df, mean)
    names(a) <- c("Década", "Crecimiento anual medio")
    a
  }, digits = 3)

  descomp <- reactive({
    y <- serie()
    validate(need(freq() > 1 && length(y) >= 2 * freq(), "La descomposición necesita estacionalidad y al menos dos ciclos."))
    validate(need(input$desc_tipo == "additive" || all(y > 0), "La multiplicativa necesita valores positivos."))
    decompose(y, type = input$desc_tipo)
  })

  output$p_desc <- renderPlot({
    d <- descomp(); t <- as.numeric(time(d$x)); n <- length(t)
    df <- data.frame(Tiempo = rep(t, 4),
                     Valor = c(d$x, d$trend, d$seasonal, d$random),
                     Comp = factor(rep(c("Observada", "Tendencia", "Estacionalidad", "Aleatoria"), each = n),
                                   levels = c("Observada", "Tendencia", "Estacionalidad", "Aleatoria")))
    ggplot(na.omit(df), aes(Tiempo, Valor)) + geom_line(color = "steelblue4") +
      facet_grid(Comp ~ ., scales = "free_y") + labs(x = "Tiempo", y = "") + tema
  })

  output$t_indices <- renderTable({
    d <- descomp()
    data.frame(Periodo = nombre_periodo(freq()), `Índice` = d$figure)
  }, digits = 4)

  output$t_fuerza <- renderTable({
    d <- descomp()
    if (input$desc_tipo == "additive") {
      ft <- max(0, 1 - var(d$random, na.rm = TRUE) / var(d$trend + d$random, na.rm = TRUE))
      fs <- max(0, 1 - var(d$random, na.rm = TRUE) / var(d$seasonal + d$random, na.rm = TRUE))
    } else {
      ft <- max(0, 1 - var(log(d$random), na.rm = TRUE) / var(log(d$trend * d$random), na.rm = TRUE))
      fs <- max(0, 1 - var(log(d$random), na.rm = TRUE) / var(log(d$seasonal * d$random), na.rm = TRUE))
    }
    data.frame(Componente = c("Tendencia", "Estacionalidad"), Fuerza = c(ft, fs))
  }, digits = 4)

  hw_suav <- reactive({
    y <- serie(); f <- freq()
    if (input$hw_opt) {
      if (f > 1) HoltWinters(y, seasonal = input$hw_tipo) else HoltWinters(y, gamma = FALSE)
    } else {
      if (f > 1) HoltWinters(y, alpha = input$hw_a, beta = input$hw_b, gamma = input$hw_g, seasonal = input$hw_tipo)
      else HoltWinters(y, alpha = input$hw_a, beta = input$hw_b, gamma = FALSE)
    }
  })

  output$t_hw_par <- renderTable({
    h <- hw_suav()
    fv <- h$fitted[, "xhat"]; y <- window(serie(), start = start(fv))
    data.frame(`Parámetro` = c("α", "β", "γ", "RMSE", "MAPE (%)"),
               Valor = c(h$alpha, h$beta, ifelse(is.logical(h$gamma), NA, h$gamma),
                         sqrt(mean((y - fv)^2)), mean(abs((y - fv) / y)) * 100))
  }, digits = 4)

  output$p_hw <- renderPlotly({
    h <- hw_suav(); y <- serie()
    fv <- h$fitted[, "xhat"]
    fc <- forecast(h, h = 2 * freq(), level = 95)
    df <- rbind(data.frame(Tiempo = as.numeric(time(y)), Valor = as.numeric(y), Serie = "Observada"),
                data.frame(Tiempo = as.numeric(time(fv)), Valor = as.numeric(fv), Serie = "Suavizada (Holt-Winters)"),
                data.frame(Tiempo = as.numeric(time(fc$mean)), Valor = as.numeric(fc$mean), Serie = "Pronóstico"))
    dfb <- data.frame(Tiempo = as.numeric(time(fc$mean)), Li = as.numeric(fc$lower), Ls = as.numeric(fc$upper))
    g <- ggplot() +
      geom_ribbon(data = dfb, aes(Tiempo, ymin = Li, ymax = Ls), fill = "skyblue", alpha = 0.5) +
      geom_line(data = df, aes(Tiempo, Valor, color = Serie)) +
      scale_color_manual(values = c("Observada" = "grey40", "Suavizada (Holt-Winters)" = "firebrick", "Pronóstico" = "blue")) +
      labs(x = "Tiempo", y = "Valor", color = "") + tema
    plotly_h(g)
  })

  # ---------- 3. Varianza y estacionariedad ----------
  output$o_levene <- renderPrint({
    y <- serie()
    g <- cut(seq_along(y), breaks = input$lev_k)
    pl <- leveneTest(as.numeric(y) ~ g)
    print(pl)
    cat(if (pl$`Pr(>F)`[1] < 0.05) "\nSe rechaza H0: la varianza cambia entre ventanas (considerar Box-Cox)."
        else "\nNo se rechaza H0: la varianza es homogénea entre ventanas.")
  })

  output$p_lev <- renderPlot({
    y <- serie()
    g <- cut(seq_along(y), breaks = input$lev_k, labels = FALSE)
    ggplot(data.frame(V = factor(g), Valor = as.numeric(y)), aes(V, Valor)) +
      geom_boxplot(fill = "skyblue", color = "darkblue") +
      labs(x = "Ventana", y = "Valor") + tema +
      theme(axis.text.x = element_text(size = 7))
  })

  output$t_bc <- renderTable({
    y <- serie()
    validate(need(all(y > 0), "Box-Cox necesita valores positivos."))
    pt <- tryCatch(car::powerTransform(as.numeric(y)), error = function(e) NULL)
    data.frame(`Método` = c("Guerrero (forecast)", "Máxima verosimilitud (forecast)", "Máxima verosimilitud (car)"),
               `λ` = c(BoxCox.lambda(y, "guerrero"), BoxCox.lambda(y, "loglik"),
                     if (is.null(pt)) NA else unname(pt$lambda)))
  }, digits = 4)

  output$p_bc <- renderPlot({
    y <- serie()
    validate(need(all(y > 0), ""))
    b <- MASS::boxcox(as.numeric(y) ~ 1, lambda = seq(-6, 4, 0.05), plotit = FALSE)
    df <- data.frame(l = b$x, ll = b$y)
    ggplot(df, aes(l, ll)) + geom_line(color = "darkblue", linewidth = 1) +
      geom_hline(yintercept = max(df$ll) - qchisq(0.95, 1) / 2, linetype = "dashed", color = "red") +
      labs(x = expression(lambda), y = "Log-verosimilitud", title = "Perfil de Box-Cox (línea roja: IC 95 %)") + tema
  })

  serie_dif <- reactive({
    y <- serie(); f <- freq()
    if (!is.null(lambda())) y <- BoxCox(y, lambda())
    d <- as.numeric(input$dif_d); D <- as.numeric(input$dif_D)
    if (D > 0 && f > 1) y <- diff(y, lag = f, differences = D)
    if (d > 0) y <- diff(y, differences = d)
    y
  })

  nombre_dif <- reactive({
    d <- as.numeric(input$dif_d); D <- as.numeric(input$dif_D); f <- freq()
    s <- c(if (d == 1) "(1-B)" else if (d == 2) "(1-B)²", if (D == 1 && f > 1) paste0("(1-B^", f, ")"))
    if (length(s) == 0) "Serie sin diferenciar" else paste(s, collapse = "")
  })

  output$t_ndiffs <- renderTable({
    y <- serie(); f <- freq()
    if (!is.null(lambda())) y <- BoxCox(y, lambda())
    data.frame(`Función` = c("ndiffs (KPSS)", "ndiffs (ADF)", if (f > 1) "nsdiffs"),
               Sugerencia = c(ndiffs(y, test = "kpss"), ndiffs(y, test = "adf"), if (f > 1) nsdiffs(y)))
  })

  output$p_dif <- renderPlot({
    yd <- serie_dif()
    ggplot(data.frame(Tiempo = as.numeric(time(yd)), Valor = as.numeric(yd)), aes(Tiempo, Valor)) +
      geom_line(color = "steelblue4") + geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      labs(title = nombre_dif(), x = "Tiempo", y = "") + tema
  })

  output$t_estac <- renderDT({
    y <- serie(); f <- freq()
    if (!is.null(lambda())) y <- BoxCox(y, lambda())
    tab <- rbind(fila_estac(y, "Original"),
                 fila_estac(diff(y), "(1-B)"),
                 if (f > 1) fila_estac(diff(y, lag = f), paste0("(1-B^", f, ")")),
                 if (f > 1) fila_estac(diff(diff(y), lag = f), paste0("(1-B)(1-B^", f, ")")),
                 fila_estac(serie_dif(), paste("Seleccionada:", nombre_dif())))
    datatable(tab, rownames = FALSE, options = list(dom = "t", ordering = FALSE)) %>%
      formatStyle("Conclusión", color = styleEqual(c("Estacionaria", "No estacionaria"), c("darkgreen", "firebrick")),
                  fontWeight = "bold")
  })

  output$p_dif_acf <- renderPlot({
    graf_acf_pacf(serie_dif(), input$dif_lags, paste("ACF y PACF de", nombre_dif()), freq())
  })

  # ---------- 4. Identificación ----------
  grid_res <- eventReactive(input$g_run, {
    y <- serie(); f <- freq()
    comb <- expand.grid(p = 0:input$g_p, q = 0:input$g_q,
                        P = if (f > 1) 0:input$g_P else 0, Q = if (f > 1) 0:input$g_Q else 0)
    res <- vector("list", nrow(comb))
    withProgress(message = "Ajustando modelos", value = 0, {
      for (k in seq_len(nrow(comb))) {
        incProgress(1 / nrow(comb))
        o <- c(comb$p[k], input$g_d, comb$q[k]); s <- c(comb$P[k], input$g_D, comb$Q[k])
        m <- tryCatch(Arima(y, order = o, seasonal = list(order = s, period = f), method = "CSS-ML",
                            lambda = lambda()), error = function(e) NULL)
        if (is.null(m)) next
        T <- length(residuals(m)); np <- length(coef(m))
        rz <- raices_arima(m)
        ok <- all(Mod(c(rz$ar, rz$ma)) > 1.001)
        res[[k]] <- data.frame(Modelo = sprintf("(%d,%d,%d)(%d,%d,%d)", o[1], o[2], o[3], s[1], s[2], s[3]),
                               p = o[1], q = o[3], P = s[1], Q = s[3], `Parámetros` = np,
                               AIC = -2 * m$loglik / T + 2 * np / T,
                               BIC = -2 * m$loglik / T + np * log(T) / T,
                               AICc = m$aicc, `AIC (R)` = m$aic,
                               `Raíces fuera del círculo` = ifelse(ok, "Sí", "No"), check.names = FALSE)
      }
    })
    do.call(rbind, res)
  })

  output$t_grid <- renderDT({
    validate(need(input$g_run > 0, "Pulse 'Buscar modelos' para empezar (puede tardar unos segundos)."))
    g <- grid_res()
    g <- g[order(g[[input$g_crit]]), ]
    datatable(g, rownames = FALSE, options = list(pageLength = 10)) %>%
      formatRound(c("AIC", "BIC"), 5) %>% formatRound(c("AICc", "AIC (R)"), 2)
  })

  auto_res <- eventReactive(input$auto_run, auto.arima(serie(), lambda = lambda(), trace = FALSE))
  output$o_auto <- renderPrint({
    validate(need(input$auto_run > 0, "Pulse el botón para ejecutar auto.arima()."))
    auto_res()
  })

  # ---------- 5. Estimación ----------
  output$t_ic <- renderDT({
    a <- ajustes(); y <- serie()
    tab <- do.call(rbind, lapply(names(a), function(m) {
      o <- a[[m]]
      ll <- if (o$tipo %in% c("arima", "ets")) o$fit$loglik else NA
      aic <- if (o$tipo %in% c("arima", "ets")) o$fit$aic else NA
      bic <- if (o$tipo %in% c("arima", "ets")) o$fit$bic else NA
      aicc <- if (o$tipo %in% c("arima", "ets")) o$fit$aicc else NA
      r <- residuales_diag(o, y)
      data.frame(Modelo = etiquetas()[m], `Especificación` = nombre_modelo(o), `Parámetros` = n_parametros(o),
                 `RMSE dentro de muestra` = sqrt(mean(r^2, na.rm = TRUE)), LogLik = ll, AIC = aic, AICc = aicc, BIC = bic,
                 check.names = FALSE)
    }))
    datatable(tab, rownames = FALSE, options = list(dom = "t", ordering = TRUE)) %>%
      formatRound(c("RMSE dentro de muestra"), 4) %>% formatRound(c("LogLik", "AIC", "AICc", "BIC"), 2)
  })

  output$t_coef <- renderDT({
    a <- ajustes()
    tab <- do.call(rbind, lapply(names(a), function(m) {
      o <- a[[m]]
      if (o$tipo == "arima") {
        cf <- coef(o$fit)
        if (length(cf) == 0) return(NULL)
        se <- sqrt(diag(o$fit$var.coef))
        z <- cf / se
        data.frame(Modelo = etiquetas()[m], `Parámetro` = names(cf), `Estimación` = cf, `Error estándar` = se,
                   z = z, `p-valor` = 2 * pnorm(-abs(z)), check.names = FALSE)
      } else if (o$tipo == "ets") {
        p <- o$fit$par[names(o$fit$par) %in% c("alpha", "beta", "gamma", "phi")]
        data.frame(Modelo = etiquetas()[m], `Parámetro` = names(p), `Estimación` = p, `Error estándar` = NA,
                   z = NA, `p-valor` = NA, check.names = FALSE)
      } else if (o$tipo == "hw") {
        data.frame(Modelo = etiquetas()[m], `Parámetro` = c("alpha", "beta", "gamma"),
                   `Estimación` = c(o$fit$alpha, o$fit$beta, ifelse(is.logical(o$fit$gamma), NA, o$fit$gamma)),
                   `Error estándar` = NA, z = NA, `p-valor` = NA, check.names = FALSE)
      } else NULL
    }))
    validate(need(!is.null(tab), "No hay modelos con parámetros."))
    tab$Significativo <- ifelse(is.na(tab$`p-valor`), "", ifelse(tab$`p-valor` < 0.05, "Sí", "No"))
    datatable(tab, rownames = FALSE, class = "compact", options = list(pageLength = 15, dom = "tp", scrollX = TRUE)) %>%
      formatRound(c("Estimación", "Error estándar", "z"), 4) %>% formatSignif("p-valor", 4) %>%
      formatStyle("Significativo", color = styleEqual(c("Sí", "No"), c("darkgreen", "firebrick")), fontWeight = "bold")
  })

  output$p_raices <- renderPlot({
    a <- ajustes()
    a <- a[sapply(a, function(o) o$tipo == "arima")]
    validate(need(length(a) > 0, "Seleccione al menos un SARIMA."))
    df <- do.call(rbind, lapply(names(a), function(m) {
      o <- a[[m]]$fit
      rz <- raices_arima(o)
      ar <- rz$ar
      ma <- rz$ma
      rbind(if (length(ar)) data.frame(Modelo = etiquetas()[m], Tipo = "AR", Re = Re(1 / ar), Im = Im(1 / ar)),
            if (length(ma)) data.frame(Modelo = etiquetas()[m], Tipo = "MA", Re = Re(1 / ma), Im = Im(1 / ma)))
    }))
    validate(need(!is.null(df), "Los modelos no tienen parte AR ni MA."))
    circ <- data.frame(x = cos(seq(0, 2 * pi, length.out = 200)), y = sin(seq(0, 2 * pi, length.out = 200)))
    ggplot() + geom_path(data = circ, aes(x, y), color = "grey50") +
      geom_hline(yintercept = 0, color = "grey80") + geom_vline(xintercept = 0, color = "grey80") +
      geom_point(data = df, aes(Re, Im, color = Tipo, shape = Tipo), size = 2.5) +
      scale_color_manual(values = c(AR = "firebrick", MA = "steelblue")) +
      facet_wrap(~ Modelo) + coord_equal() + labs(x = "Real", y = "Imaginaria", color = "", shape = "") + tema
  })

  # ---------- 6. Diagnóstico ----------
  output$ui_diag_mod <- renderUI({
    a <- ajustes()
    selectInput("diag_mod", "Modelo", choices = setNames(names(a), etiquetas()[names(a)]),
                selected = if ("M1" %in% names(a)) "M1" else names(a)[1])
  })

  output$p_diag <- renderPlot({
    req(input$diag_mod)
    a <- ajustes(); o <- a[[input$diag_mod]]; req(o)
    r <- residuales_diag(o, serie())
    rr <- as.numeric(na.omit(r)); lags <- input$diag_lags
    re <- (rr - mean(rr)) / sd(rr); ci <- qnorm(0.975) / sqrt(length(rr))
    tt <- as.numeric(time(na.omit(r)))
    g1 <- ggplot(data.frame(t = tt, r = re), aes(t, r)) + geom_line(color = "steelblue4") +
      geom_hline(yintercept = c(-3, 0, 3), linetype = c("dashed", "solid", "dashed"), color = c("red", "black", "red")) +
      labs(title = "Residuales estandarizados", x = "", y = "") + tema
    g2 <- ggplot(data.frame(r = re), aes(r)) +
      geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "skyblue", color = "darkblue") +
      geom_density(color = "red", linewidth = 1) +
      stat_function(fun = dnorm, color = "darkgreen", linetype = "dashed", linewidth = 1) +
      labs(title = "Histograma (rojo: densidad, verde: normal)", x = "", y = "") + tema
    av <- acf(rr, lag.max = lags, plot = FALSE)$acf[-1]
    pv <- pacf(rr, lag.max = lags, plot = FALSE)$acf[, 1, 1]
    g3 <- ggplot(data.frame(l = 1:lags, v = av), aes(l, v)) + geom_segment(aes(xend = l, yend = 0), color = "steelblue", linewidth = 1) +
      geom_hline(yintercept = c(-ci, ci), linetype = "dashed", color = "red") + geom_hline(yintercept = 0) +
      labs(title = "ACF", x = "Rezago", y = "") + tema
    g4 <- ggplot(data.frame(l = 1:lags, v = pv), aes(l, v)) + geom_segment(aes(xend = l, yend = 0), color = "steelblue", linewidth = 1) +
      geom_hline(yintercept = c(-ci, ci), linetype = "dashed", color = "red") + geom_hline(yintercept = 0) +
      labs(title = "PACF", x = "Rezago", y = "") + tema
    k <- n_parametros(o) * (o$tipo == "arima")
    lb <- sapply(1:lags, function(h) if (h <= k) NA else Box.test(rr, h, "Ljung-Box", fitdf = k)$p.value)
    g5 <- ggplot(na.omit(data.frame(l = 1:lags, p = lb)), aes(l, p)) +
      geom_point(aes(color = p < 0.05), size = 2) +
      scale_color_manual(values = c("FALSE" = "darkblue", "TRUE" = "red"), guide = "none") +
      geom_hline(yintercept = 0.05, linetype = "dashed", color = "blue") + ylim(0, 1) +
      labs(title = "p-valores Ljung-Box", x = "Rezago", y = "") + tema
    qq <- qqnorm(re, plot.it = FALSE)
    g6 <- ggplot(data.frame(x = qq$x, y = qq$y), aes(x, y)) + geom_point(color = "steelblue4", alpha = 0.7) +
      geom_abline(color = "red") + labs(title = "Q-Q normal", x = "Teóricos", y = "Muestrales") + tema
    gridExtra::grid.arrange(g1, g2, g3, g4, g5, g6, ncol = 2,
                            top = paste("Diagnóstico residual:", etiquetas()[input$diag_mod]))
  })

  output$t_diag <- renderDT({
    req(input$diag_mod)
    o <- ajustes()[[input$diag_mod]]; req(o)
    k <- n_parametros(o) * (o$tipo == "arima")
    tab <- pruebas_residuales(residuales_diag(o, serie()), k, freq())
    datatable(tab, rownames = FALSE, options = list(dom = "t", ordering = FALSE)) %>%
      formatStyle("Decisión al 5 %", color = styleEqual(c("Se rechaza H0", "No se rechaza H0"), c("firebrick", "darkgreen")),
                  fontWeight = "bold")
  })

  output$t_diag_todos <- renderDT({
    a <- ajustes(); y <- serie()
    tab <- do.call(rbind, lapply(names(a), function(m) {
      o <- a[[m]]; k <- n_parametros(o) * (o$tipo == "arima")
      p <- tryCatch(pruebas_residuales(residuales_diag(o, y), k, freq()), error = function(e) NULL)
      if (is.null(p)) return(NULL)
      v <- setNames(as.list(p$`p-valor`), p$Prueba)
      data.frame(Modelo = etiquetas()[m], v, check.names = FALSE)
    }))
    cols <- setdiff(names(tab), "Modelo")
    datatable(tab, rownames = FALSE, options = list(dom = "t", scrollX = TRUE)) %>%
      formatRound(cols, 3) %>%
      formatStyle(cols, color = styleInterval(0.05, c("firebrick", "black")))
  })

  # ---------- 7. Pronóstico ----------
  output$ui_fc_mod <- renderUI({
    a <- ajustes()
    selectInput("fc_mod", "Modelo para la tabla", choices = setNames(names(a), etiquetas()[names(a)]),
                selected = if ("M1" %in% names(a)) "M1" else names(a)[1])
  })

  pronosticos <- reactive({
    a <- ajustes()
    lapply(a, function(o) tryCatch(pronosticar(o, input$fc_h, input$fc_nivel), error = function(e) NULL))
  })

  output$p_fc <- renderPlotly({
    fcs <- pronosticos(); y <- serie()
    fcs <- fcs[!sapply(fcs, is.null)]
    ini <- max(time(y)) - input$fc_hist
    dfo <- subset(data.frame(Tiempo = as.numeric(time(y)), Valor = as.numeric(y)), Tiempo > ini)
    df <- do.call(rbind, lapply(names(fcs), function(m) {
      f <- fcs[[m]]
      data.frame(Tiempo = as.numeric(time(f$mean)), Pron = as.numeric(f$mean),
                 Li = as.numeric(f$lower[, 1]), Ls = as.numeric(f$upper[, 1]), Modelo = etiquetas()[m])
    }))
    sel <- etiquetas()[input$fc_mod]
    bandas <- if (input$fc_bandas) df else subset(df, Modelo == sel)
    g <- ggplot() +
      geom_ribbon(data = bandas, aes(Tiempo, ymin = Li, ymax = Ls, fill = Modelo), alpha = 0.2) +
      geom_line(data = dfo, aes(Tiempo, Valor), color = "black") +
      geom_line(data = df, aes(Tiempo, Pron, color = Modelo), linewidth = 0.8) +
      scale_color_manual(values = colores_mod()) + scale_fill_manual(values = colores_mod(), guide = "none") +
      labs(title = paste0("Pronósticos con bandas del ", input$fc_nivel, " %"), x = "Tiempo", y = "Valor", color = "") + tema
    plotly_h(g)
  })

  tabla_fc <- reactive({
    req(input$fc_mod)
    f <- pronosticos()[[input$fc_mod]]; req(f)
    data.frame(Periodo = etiqueta_tiempo(f$mean), h = seq_along(f$mean), `Pronóstico` = as.numeric(f$mean),
               `Límite inferior` = as.numeric(f$lower[, 1]), `Límite superior` = as.numeric(f$upper[, 1]),
               check.names = FALSE)
  })

  output$t_fc <- renderDT({
    datatable(tabla_fc(), rownames = FALSE, options = list(pageLength = 12, dom = "tp")) %>%
      formatRound(3:5, 3)
  })

  output$fc_desc <- downloadHandler(
    filename = function() paste0("pronosticos_", input$fc_mod, ".csv"),
    content = function(file) write.csv(tabla_fc(), file, row.names = FALSE)
  )

  # ---------- 8. Capacidad de pronóstico ----------
  escala_mase <- reactive({
    y <- serie(); f <- freq()
    mean(abs(diff(y, lag = max(f, 1))))
  })

  output$t_in <- renderDT({
    a <- ajustes(); y <- serie()
    tab <- do.call(rbind, lapply(names(a), function(m) {
      fv <- valores_ajustados(a[[m]], y)
      ok <- !is.na(fv)
      cbind(Modelo = etiquetas()[m], metricas(as.numeric(y[ok]), as.numeric(fv[ok]), escala = escala_mase()))
    }))
    datatable(tab, rownames = FALSE, options = list(dom = "t", order = list(list(2, "asc")))) %>%
      formatRound(2:6, 4)
  })

  output$p_in <- renderPlotly({
    a <- ajustes(); y <- serie()
    ini <- max(time(y)) - 6
    dfo <- subset(data.frame(Tiempo = as.numeric(time(y)), Valor = as.numeric(y)), Tiempo > ini)
    df <- do.call(rbind, lapply(names(a), function(m) {
      fv <- valores_ajustados(a[[m]], y)
      data.frame(Tiempo = as.numeric(time(y)), Valor = as.numeric(fv), Modelo = etiquetas()[m])
    }))
    df <- subset(na.omit(df), Tiempo > ini)
    g <- ggplot() + geom_line(data = dfo, aes(Tiempo, Valor), color = "black", linewidth = 1) +
      geom_line(data = df, aes(Tiempo, Valor, color = Modelo), linetype = "dashed") +
      scale_color_manual(values = colores_mod()) +
      labs(title = "Trayectoria ajustada (últimos 6 años; negro = serie)", x = "Tiempo", y = "Valor", color = "") + tema
    plotly_h(g)
  })

  particion <- reactive({
    y <- serie(); n <- length(y)
    n_tr <- round(input$pt_train / 100 * n)
    list(train = ts(y[1:n_tr], start = start(y), frequency = freq()),
         test = ts(y[(n_tr + 1):n], start = time(y)[n_tr + 1], frequency = freq()), n_tr = n_tr)
  })

  output$pt_info <- renderText({
    p <- particion()
    paste0("Entrenamiento: ", p$n_tr, " obs. (hasta ", tail(etiqueta_tiempo(p$train), 1), "). Prueba: ",
           length(p$test), " obs.")
  })

  fc_particion <- reactive({
    p <- particion(); h <- length(p$test)
    withProgress(message = "Pronosticando el periodo de prueba", value = 0, {
      out <- list()
      for (m in input$metodos) {
        incProgress(1 / length(input$metodos), detail = etiquetas()[m])
        out[[m]] <- tryCatch(pronostico_submuestra(m, p$train, h, 95), error = function(e) NULL)
      }
      out[!sapply(out, is.null)]
    })
  })

  output$t_pt <- renderDT({
    fcs <- fc_particion(); p <- particion(); real <- as.numeric(p$test)
    esc <- mean(abs(diff(p$train, lag = max(freq(), 1))))
    tab <- do.call(rbind, lapply(names(fcs), function(m) {
      f <- fcs[[m]]
      cbind(Modelo = etiquetas()[m], metricas(real, as.numeric(f$mean), as.numeric(f$lower), as.numeric(f$upper), esc))
    }))
    datatable(tab, rownames = FALSE, options = list(dom = "t", order = list(list(2, "asc")))) %>%
      formatRound(2:6, 4) %>% formatRound(7, 1)
  })

  output$p_pt <- renderPlotly({
    fcs <- fc_particion(); p <- particion(); y <- serie()
    ini <- as.numeric(time(p$test))[1] - 5
    dfo <- subset(data.frame(Tiempo = as.numeric(time(y)), Valor = as.numeric(y)), Tiempo > ini)
    df <- do.call(rbind, lapply(names(fcs), function(m)
      data.frame(Tiempo = as.numeric(time(fcs[[m]]$mean)), Valor = as.numeric(fcs[[m]]$mean), Modelo = etiquetas()[m])))
    g <- ggplot() + geom_line(data = dfo, aes(Tiempo, Valor), color = "black", linewidth = 1) +
      geom_line(data = df, aes(Tiempo, Valor, color = Modelo)) +
      geom_vline(xintercept = as.numeric(time(p$test))[1], linetype = "dashed") +
      scale_color_manual(values = colores_mod()) +
      labs(title = "Periodo de prueba: real (negro) vs pronósticos", x = "Tiempo", y = "Valor", color = "") + tema
    plotly_h(g)
  })

  cv_res <- eventReactive(input$cv_run, {
    y <- serie(); n <- length(y); k <- input$cv_k
    min_tr <- floor(n * input$cv_min / 100)
    step <- floor((n - min_tr) / k)
    validate(need(step >= 1, "Muy pocos datos para ese número de pliegues."))
    res <- list()
    withProgress(message = "Validación cruzada", value = 0, {
      for (i in 1:k) {
        tr_end <- min_tr + (i - 1) * step
        te_end <- if (i == k) n else min_tr + i * step
        if (input$cv_h == "fijo") te_end <- min(n, tr_end + input$cv_hfijo)
        y_tr <- ts(y[1:tr_end], start = start(y), frequency = freq())
        real <- as.numeric(y[(tr_end + 1):te_end])
        for (m in input$metodos) {
          incProgress(1 / (k * length(input$metodos)), detail = paste("Pliegue", i, "-", etiquetas()[m]))
          pr <- tryCatch(as.numeric(pronostico_submuestra(m, y_tr, length(real))$mean),
                         error = function(e) rep(NA, length(real)))
          res[[length(res) + 1]] <- data.frame(Pliegue = i, Modelo = etiquetas()[m],
                                               MAPE = mean(abs((real - pr) / real)) * 100,
                                               RMSE = sqrt(mean((real - pr)^2)),
                                               MAE = mean(abs(real - pr)))
        }
      }
    })
    do.call(rbind, res)
  })

  output$p_cv <- renderPlot({
    validate(need(input$cv_run > 0, "Pulse 'Ejecutar validación cruzada' (puede tardar alrededor de un minuto)."))
    d <- cv_res()
    ggplot(d, aes(factor(Pliegue), RMSE, color = Modelo, group = Modelo)) + geom_line() + geom_point(size = 2) +
      scale_color_manual(values = colores_mod()) + scale_y_log10() +
      labs(title = "RMSE por pliegue (escala logarítmica)", x = "Pliegue", y = "RMSE", color = "") + tema
  })

  output$t_cv <- renderDT({
    req(input$cv_run > 0)
    d <- cv_res()
    gan <- sapply(split(d, d$Pliegue), function(x) x$Modelo[which.min(x$RMSE)])
    tab <- do.call(rbind, lapply(split(d, d$Modelo), function(x)
      data.frame(Modelo = x$Modelo[1], `MAPE medio` = mean(x$MAPE, na.rm = TRUE), `MAPE DE` = sd(x$MAPE, na.rm = TRUE),
                 `RMSE medio` = mean(x$RMSE, na.rm = TRUE), `RMSE DE` = sd(x$RMSE, na.rm = TRUE),
                 `Pliegues ganados` = sum(gan == x$Modelo[1]), check.names = FALSE)))
    datatable(tab, rownames = FALSE, options = list(dom = "t", order = list(list(3, "asc")))) %>%
      formatRound(2:5, 4)
  })

  output$p_cv_box <- renderPlot({
    req(input$cv_run > 0)
    d <- cv_res()
    ggplot(d, aes(reorder(Modelo, RMSE, median), RMSE, fill = Modelo)) + geom_boxplot(alpha = 0.6) +
      scale_fill_manual(values = colores_mod(), guide = "none") + scale_y_log10() + coord_flip() +
      labs(x = "", y = "RMSE (escala log)") + tema
  })

  un_paso <- reactive({
    y <- serie(); n <- length(y); f <- freq()
    n_tr <- round(input$up_train / 100 * n)
    y_tr <- ts(y[1:n_tr], start = start(y), frequency = f)
    idx <- (n_tr + 1):n
    withProgress(message = "Pronóstico a un paso", value = 0, {
      out <- list()
      for (m in input$metodos) {
        incProgress(1 / length(input$metodos), detail = etiquetas()[m])
        pr <- tryCatch({
          o <- ajustar_metodo(m, y_tr, lambda(), custom_spec(), input$hw_tipo)
          if (o$tipo == "arima") {
            as.numeric(fitted(Arima(y, model = o$fit)))[idx]
          } else if (o$tipo == "ets") {
            as.numeric(fitted(ets(y, model = o$fit, use.initial.values = TRUE)))[idx]
          } else if (o$tipo == "hw") {
            h <- o$fit
            sapply(idx, function(t) {
              yy <- ts(y[1:(t - 1)], start = start(y), frequency = f)
              if (!is.null(lambda())) yy <- BoxCox(yy, lambda())
              hh <- if (f > 1) HoltWinters(yy, alpha = h$alpha, beta = h$beta, gamma = h$gamma, seasonal = input$hw_tipo)
                    else HoltWinters(yy, alpha = h$alpha, beta = h$beta, gamma = FALSE)
              p1 <- as.numeric(predict(hh, 1))
              if (!is.null(lambda())) InvBoxCox(p1, lambda()) else p1
            })
          } else {
            as.numeric(y[idx - f])
          }
        }, error = function(e) NULL)
        if (!is.null(pr)) out[[m]] <- data.frame(Tiempo = as.numeric(time(y))[idx], Real = as.numeric(y[idx]),
                                                 Pron = pr, Modelo = etiquetas()[m])
      }
      out
    })
  })

  output$t_up <- renderDT({
    u <- un_paso()
    tab <- do.call(rbind, lapply(u, function(d) cbind(Modelo = d$Modelo[1], metricas(d$Real, d$Pron, escala = escala_mase()))))
    datatable(tab, rownames = FALSE, options = list(dom = "t", order = list(list(2, "asc")))) %>% formatRound(2:6, 4)
  })

  output$p_up <- renderPlotly({
    u <- un_paso()
    df <- do.call(rbind, u)
    df$Error <- df$Real - df$Pron
    g <- ggplot(df, aes(Tiempo, Error, color = Modelo)) + geom_line(alpha = 0.8) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      scale_color_manual(values = colores_mod()) +
      labs(title = "Errores de pronóstico a un paso en el periodo de prueba", x = "Tiempo", y = "Real − pronóstico", color = "") + tema
    plotly_h(g)
  })
}

shinyApp(ui, server)
