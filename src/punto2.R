# =====================================================================
# Taller 1 - Punto 2 : filtro de media movil sobre distintas series
# Series de Tiempo I - UdeA
# =====================================================================

# ---- (a) AR(2): x_t = -0.9 x_{t-2} + w_t ----------------------------
set.seed(123)

# 150 y no 100: los primeros 50 son burn-in. La recursion arranca en
# x = 0, un valor artificial; se descartan para que la serie ya este
# en su regimen estacionario.
w <- rnorm(150, mean = 0, sd = 1)
x <- stats::filter(w, filter = c(0, -0.9), method = "recursive")[-(1:50)]
x <- ts(x)

# filtro de media movil unilateral de orden 4
v <- stats::filter(x, filter = rep(1/4, 4), sides = 1)

plot(x, type = "l", ylab = expression(x[t]),
     main = "(a)  AR(2) con phi2 = -0.9 y su media movil de orden 4")
lines(v, lty = 2, lwd = 2, col = "blue")
legend("topright", c("x_t", "v_t (MA 4)"), lty = c(1, 2),
       col = c("black", "blue"), bty = "n")

sd(x); sd(v, na.rm = TRUE)   # esperado: 1.5012  y  0.3531
dev.cur()
list.files(pattern = "Rplots")
ls()
# ---- (b) coseno puro de periodo 4 -----------------------------------
t  <- 1:100
x2 <- ts(cos(2 * pi * t / 4))
v2 <- stats::filter(x2, filter = rep(1/4, 4), sides = 1)

plot(x2, type = "l", ylim = c(-1.1, 1.1), ylab = expression(x[t]),
     main = "(b)  cos(2*pi*t/4) y su media movil de orden 4")
lines(v2, lty = 2, lwd = 2, col = "blue")
abline(h = 0, col = "gray70")
legend("topright", c("x_t", "v_t (MA 4)"), lty = c(1, 2),
       col = c("black", "blue"), bty = "n")

x2[1:8]
range(v2, na.rm = TRUE)

# ---- (c) coseno de periodo 4 + ruido N(0,1) -------------------------
set.seed(123)
t  <- 1:100
w3 <- rnorm(100, mean = 0, sd = 1)
x3 <- ts(cos(2 * pi * t / 4) + w3)
v3 <- stats::filter(x3, filter = rep(1/4, 4), sides = 1)

plot(x3, type = "l", ylab = expression(x[t]),
     main = "(c)  cos(2*pi*t/4) + w_t  y su media movil de orden 4")
lines(v3, lty = 2, lwd = 2, col = "blue")
abline(h = 0, col = "gray70")
legend("topright", c("x_t", "v_t (MA 4)"), lty = c(1, 2),
       col = c("black", "blue"), bty = "n")

sd(x3); sd(v3, na.rm = TRUE)
acf(v3[!is.na(v3)], lag.max = 6, main = "ACF de v_t en (c)")
# =====================================================================
# COMENTARIOS
# =====================================================================
#
# --- (a) AR(2): x_t = -0.9 x_{t-2} + w_t -----------------------------
# x_t oscila de forma pseudoperiodica con periodo aproximado 4. El
# polinomio caracteristico 1 + 0.9 z^2 tiene raices complejas de modulo
# 1/sqrt(0.9) = 1.054, apenas por fuera del circulo unitario: raices
# complejas -> oscilacion, modulo cercano a 1 -> muy persistente.
# La ACF cumple rho(h) = -0.9 rho(h-2): cero en rezagos impares,
# rho(2) = -0.9, rho(4) = 0.81.
#
# La media movil casi elimina la oscilacion: sd(x) = 1.50 vs
# sd(v) = 0.35, una reduccion de 4 a 6 veces segun la realizacion.
# La reduccion NO viene de promediar ruido: si x fuera ruido blanco de
# la misma varianza, Var(v) seria gamma(0)/4 = 1.32, y el valor real es
# 0.132, diez veces menor. La causa es que gamma(2) = -0.9*gamma(0) es
# fuertemente negativo y los terminos se cancelan:
#   Var(v) = (1/16)[4*gamma(0) + 4*gamma(2)] = (1/16)[21.05 - 18.95]
# Queda un residuo pequeno porque el ciclo del AR(2) no es una sinusoide
# pura: su espectro esta concentrado alrededor de pi/2, no en un punto.
#
# --- (b) x_t = cos(2*pi*t/4) -----------------------------------------
# Senal deterministica de periodo exactamente 4. Muestreada en los
# enteros solo toma cuatro valores: 0, -1, 0, 1.
# El filtro la anula por completo: v_t = 0 exacto (lo que se observa,
# del orden de 1e-15, es error de redondeo de punto flotante). Cada
# ventana de 4 valores consecutivos suma 0 + (-1) + 0 + 1 = 0.
# Es el caso extremo de (a): alli habia atenuacion, aqui aniquilacion
# total, porque el periodo de la senal coincide exactamente con la
# longitud de la ventana.
#
# --- (c) x_t = cos(2*pi*t/4) + w_t -----------------------------------
# La descomposicion es lineal, asi que el filtro actua por separado
# sobre cada componente:
#   v_t = MA4(coseno) + MA4(ruido) = 0 + MA4(w_t)
# La parte periodica desaparece igual que en (b); lo que queda en v_t
# es ruido filtrado. sd(x) = 1.14 (teorico sqrt(1.5) = 1.22) frente a
# sd(v) = 0.44 (teorico sqrt(1/4) = 0.50).
# Aqui si aplica el resultado clasico: promediar 4 valores
# independientes divide la varianza entre 4.
# El filtro ademas INTRODUCE autocorrelacion donde no la habia:
# rho_v(h) = (4-h)/4, es decir 0.75, 0.50, 0.25 y cero desde h = 4,
# porque dos ventanas separadas h posiciones comparten 4-h datos.
#
# --- (d) Comparacion -------------------------------------------------
# Los tres casos son la misma explicacion vista a distintas escalas: el
# filtro de media movil de orden 4 es un filtro PASA-BAJAS cuya funcion
# de respuesta en frecuencia es
#   |H(w)| = |sin(2w)| / (4*|sin(w/2)|),   con H(0) = 1
# y que tiene CEROS exactos en w = pi/2 (periodo 4) y w = pi (periodo 2):
#
#   periodo   16      8      6      4      3      2
#   |H(w)|   0.906  0.653  0.433  0.000  0.250  0.000
#
# (a) El AR(2) concentra su energia alrededor de w = pi/2, justo sobre
#     un cero del filtro -> atenuacion fuerte pero no total, porque el
#     espectro esta repartido en una banda, no en un solo punto.
# (b) Toda la energia esta EXACTAMENTE en w = pi/2 -> se elimina todo.
# (c) La energia esta repartida: una parte exactamente en w = pi/2 (el
#     coseno, se elimina) y el resto uniforme en todas las frecuencias
#     (el ruido blanco, se atenua parcialmente). Sobrevive solo la parte
#     de baja frecuencia del ruido.
#
# Conclusion: la media movil no "suaviza" de forma generica. Elimina
# selectivamente las frecuencias donde su ganancia se anula, y como
# H(0) = 1, deja intacto el nivel y casi intacta la tendencia. Por eso
# se usa para estimar tendencia: mata la estacionalidad de periodo igual
# a la ventana y conserva el movimiento lento.