# Más allá de ARMA y ARIMA: modelos ARFIMA para series con memoria larga

Seminario del curso **Series de Tiempo I** (Estadística, Universidad de Antioquia, 2026-2).

**Integrantes:** Juan Camilo Melo Ortíz · Santiago Duque Pérez
**Profesor:** Duván Cataño

## De qué se trata

Los modelos ARMA suponen memoria corta (la ACF cae exponencialmente) y los ARIMA
suponen una raíz unitaria (los choques son permanentes). Los modelos **ARFIMA**
permiten que el orden de diferenciación `d` tome valores fraccionarios y así
describen series con **memoria larga**: autocorrelaciones pequeñas que bajan muy
despacio y se mantienen durante muchísimos rezagos.

Como aplicación usamos la serie clásica de la memoria larga: los **niveles mínimos
anuales del río Nilo entre los años 622 y 1284** (663 observaciones), y comparamos
un ARFIMA contra los mejores ARMA y ARIMA ajustados a la misma serie.

## Estructura

```
seminario-arfima/
├── notebooks/
│   └── arfima_nilo.ipynb   análisis completo (kernel de R)
├── figuras/                figuras que genera el notebook, para el documento
├── documento/              fuente LaTeX y PDF del entregable escrito
└── presentacion/           diapositivas en Beamer
```

## Resultados principales

| | Valor |
|---|---|
| Parámetro de memoria `d` (máxima verosimilitud) | 0.393 (IC 95 %: 0.334 a 0.453) |
| Exponente de Hurst `H = d + 0.5` | 0.89 |
| Prueba ADF (H0: raíz unitaria) | se rechaza, p < 0.01 |
| Prueba KPSS (H0: estacionariedad) | se rechaza, p < 0.01 |

Que las dos pruebas rechacen al mismo tiempo indica que la serie no se comporta
ni como ARMA ni como ARIMA. El modelo escogido es un **ARFIMA(0,d,0)**, con un solo
parámetro, frente a un ARMA(2,1) y un ARIMA(1,1,3).

Pronóstico de las últimas 50 observaciones (RMSE):

| Modelo | A un paso | A 50 pasos |
|---|---|---|
| ARFIMA(0,d,0) | **50.1** | 62.6 |
| ARMA(2,1) | 50.8 | 68.9 |
| ARIMA(1,1,3) | 50.7 | 74.8 |
| Media (referencia) | 60.7 | **60.7** |

## Cómo reproducir el análisis

Requisitos: R 4.x y los paquetes

```r
install.packages(c("fracdiff", "forecast", "tseries", "longmemo"))
```

El notebook los instala solo si faltan.

1. Abrir `notebooks/arfima_nilo.ipynb` en VS Code (extensión Jupyter) o en Jupyter.
2. Seleccionar el kernel de R.
3. Ejecutar todas las celdas. Las figuras se regeneran en `figuras/`.

El notebook debe ejecutarse desde la carpeta `notebooks/`, porque guarda las
figuras en la ruta relativa `../figuras`.

## Referencias

- Geweke, J. y Porter-Hudak, S. (1983). The estimation and application of long memory time series models. *Journal of Time Series Analysis*, 4(4), 221–238.
- Granger, C. W. J. y Joyeux, R. (1980). An introduction to long-memory time series models and fractional differencing. *Journal of Time Series Analysis*, 1(1), 15–29.
- Haslett, J. y Raftery, A. E. (1989). Space-time modelling with long-memory dependence: Assessing Ireland's wind power resource. *Applied Statistics*, 38(1), 1–50.
- Hosking, J. R. M. (1981). Fractional differencing. *Biometrika*, 68(1), 165–176.
- Hurst, H. E. (1951). Long-term storage capacity of reservoirs. *Transactions of the American Society of Civil Engineers*, 116, 770–808.
- Maechler, M. et al. (2026). *fracdiff: Fractionally Differenced ARIMA aka ARFIMA(P,d,q) Models*. Paquete de R, versión 1.5-4.
- Shumway, R. H. y Stoffer, D. S. (2017). *Time Series Analysis and Its Applications: With R Examples* (4a ed.). Springer.
- Toussoun, O. (1925). Mémoire sur l'histoire du Nil. *Mémoires de l'Institut d'Égypte*, 18, 366–404.
