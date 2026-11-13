# Aplicación: concentración de CO₂ en Mauna Loa con modelos SARIMA

Trabajo de aplicación del curso **Series de Tiempo I** (Estadística, Universidad de Antioquia, 2026-2).

**Autor:** Juan Camilo Melo Ortíz
**Profesor:** Duván Cataño

## La serie

`co2` viene incluida en R (paquete `datasets`): concentraciones atmosféricas mensuales
de CO₂, en partes por millón (ppm), medidas en el observatorio de Mauna Loa (Hawái)
entre enero de 1959 y diciembre de 1997. Son 468 observaciones con una tendencia
creciente y un ciclo anual muy marcado (pico en mayo, valle en septiembre-octubre).

## Ruta del análisis

Sigue la metodología de los notebooks del curso (Amazon, Cmort y Air Passengers):

1. Análisis descriptivo: tabla Buys-Ballot, perfiles estacionales y boxplots por mes.
2. Suavizamiento con medias móviles (3, 5, 2×12 y 2×24) y descomposición clásica.
3. Homogeneidad de varianzas (Levene) y Box-Cox.
4. Pruebas de estacionariedad (ADF y KPSS) y diferenciación.
5. Identificación con ACF/PACF, búsqueda por AIC/BIC y `auto.arima()`.
6. Estimación de cuatro modelos candidatos, significancia y raíces.
7. Diagnóstico residual: Ljung-Box, Shapiro-Wilk, Lilliefors, Jarque-Bera, Levene y McLeod-Li.
8. Pronóstico a 24 meses.
9. Capacidad de pronóstico: dentro de muestra, partición 80/20, validación cruzada con
   ventana expansiva (10 pliegues) y pronóstico a un paso, comparando contra
   Holt-Winters, ETS y el ingenuo estacional.

## Estructura

```
aplicacion-co2/
├── notebooks/
│   └── co2_sarima.ipynb   análisis completo (kernel de R), ya ejecutado
├── figuras/               figuras principales que guarda el notebook
└── app/
    └── app.R              aplicación Shiny para comparar modelos
```

## Resultados principales

| | Resultado |
|---|---|
| Varianza | Homogénea (Levene p ≈ 1), no se transforma |
| Serie original | No estacionaria (ADF p = 0.23, KPSS p < 0.01) |
| Diferenciación | (1−B)(1−B¹²) |
| Modelo escogido | **SARIMA(0,1,1)(0,1,1)₁₂**: θ = −0.35, Θ = −0.85, σ = 0.29 ppm |
| Residuales | Ruido blanco normal y homocedástico (todas las pruebas con p > 0.05) |

Candidatos comparados:

| Modelo | Origen | Observación |
|---|---|---|
| M1: SARIMA(0,1,1)(0,1,1)₁₂ | ACF/PACF y mínimo BIC | Escogido |
| M2: SARIMA(1,1,1)(0,1,1)₁₂ | Segundo por BIC | φ₁ no significativo (p = 0.09) |
| M3: SARIMA(1,1,1)(1,1,2)₁₂ | `auto.arima()` | Parámetros estacionales no significativos |
| M4: SARIMA(1,1,1)(2,1,2)₁₂ | Mínimo AIC | MA estacional al borde de no ser invertible |

Validación cruzada (10 pliegues, bloques de 23 meses), RMSE promedio en ppm:

| Modelo | RMSE |
|---|---|
| M3 (`auto.arima`) | 0.532 |
| **M1** | **0.539** |
| M2 | 0.543 |
| M4 | 0.564 |
| ETS | 0.612 |
| Holt-Winters | 0.626 |
| Ingenuo estacional | 2.415 |

Los SARIMA quedan prácticamente empatados entre sí y por delante de los métodos de
suavizamiento; M1 logra lo mismo con solo dos parámetros. A un paso, todos tienen un
error cercano a 0.30 ppm.

## Cómo reproducir el notebook

Requisitos: R 4.x, Jupyter con el kernel de R (`IRkernel`) y los paquetes

```r
install.packages(c("astsa", "forecast", "tseries", "nortest", "fBasics", "car",
                   "MLmetrics", "ggplot2", "gridExtra", "pastecs", "lmtest"))
```

El notebook instala solo los paquetes que falten. En VS Code basta con abrir
`notebooks/co2_sarima.ipynb` y escoger el kernel de R.

## Cómo ejecutar la app

```r
install.packages(c("shiny", "bslib", "forecast", "tseries", "nortest", "car",
                   "ggplot2", "plotly", "DT", "gridExtra"))
```

Desde PowerShell, parado en la raíz del repositorio:

```powershell
Rscript.exe -e "shiny::runApp('aplicacion-co2/app', launch.browser=TRUE)"
```

O abrir `aplicacion-co2/app/app.R` en RStudio y usar el botón **Run App**.

### Qué tiene la app

| Pestaña | Contenido |
|---|---|
| 1. Serie | Resumen, gráfico interactivo, boxplot por mes, perfil estacional y tabla Buys-Ballot |
| 2. Suavizamiento | Medias móviles centradas o unilaterales de cualquier orden, crecimiento de la tendencia por década, descomposición aditiva/multiplicativa y Holt-Winters con parámetros manuales u optimizados |
| 3. Estacionariedad | Levene por ventanas, Box-Cox (tres métodos y perfil de verosimilitud), diferenciación con d y D a escoger, ADF/KPSS, `ndiffs`/`nsdiffs` y ACF/PACF |
| 4. Identificación | Búsqueda de órdenes por AIC/BIC/AICc y `auto.arima()` |
| 5. Estimación | Criterios de información, coeficientes con significancia y raíces inversas en el círculo unitario |
| 6. Diagnóstico | Panel de seis gráficos y tabla de pruebas para el modelo escogido, más un resumen de p-valores de todos los modelos |
| 7. Pronóstico | Pronósticos de todos los modelos con bandas, tabla y descarga en CSV |
| 8. Capacidad predictiva | Métricas dentro de muestra, partición entrenamiento/prueba ajustable, validación cruzada con ventana expansiva y pronóstico a un paso |
| Metodología | Resumen de la ruta y de la conclusión |

En la barra lateral se escogen los modelos a comparar (los cuatro candidatos, un
SARIMA personalizado, `auto.arima`, Holt-Winters, ETS y el ingenuo estacional), la
transformación de Box-Cox que se aplica a todos y, si se quiere, otra serie desde un CSV.

La validación cruzada tarda alrededor de un minuto con los modelos por defecto; si se
incluye `auto.arima`, tarda bastante más porque vuelve a escoger el modelo en cada pliegue.
