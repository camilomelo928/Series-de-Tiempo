# Series de Tiempo — UdeA

Trabajos del curso de **Series de Tiempo I** (Estadística, Universidad de Antioquia).

| Proyecto | Carpeta | Qué es |
|---|---|---|
| Seminario ARFIMA | [`seminario-arfima/`](seminario-arfima/) | Documento escrito, notebook reproducible en R y presentación Beamer sobre modelos de memoria larga, con aplicación al río Nilo |
| App de suavizamiento | [`src/suavizados.R`](src/suavizados.R) | App Shiny para comparar medias móviles, SES y Holt-Winters |
| Copias de clase | [`copias-clase/`](copias-clase/) | Notebooks del curso ejecutados como práctica |

## App de suavizamiento: métodos de suavizamiento para series de tiempo

Aplicación interactiva en R (Shiny) para ajustar y comparar tres métodos clásicos
de suavizamiento sobre una serie de tiempo cargada desde un archivo CSV.

Desarrollada como trabajo del curso de **Análisis de Series de Tiempo**,
Universidad de Antioquia.

### Qué hace

Cargas un CSV con dos columnas —tiempo y observaciones— y la app ajusta el método
que elijas, grafica el resultado con su intervalo del 95 %, reporta las métricas de
error dentro de muestra y te deja pronosticar a `h` períodos.

#### Métodos implementados

| Método | Componentes que modela | Parámetros |
|---|---|---|
| Media móvil unilateral | Nivel local | Orden `k` |
| Suavización exponencial simple (SES) | Nivel | `α` (manual u optimizado) |
| Holt-Winters | Nivel, tendencia y estación | `α`, `β`, `γ` (manuales u optimizados) |

#### Pestañas

- **Descomposición** — separa la serie en tendencia, estacionalidad y residuo, con
  selector aditivo/multiplicativo y la tabla de índices estacionales estimados.
- **Pronóstico** — serie observada, valores ajustados, pronóstico con banda del 95 %,
  tabla de RMSE / MAE / MAPE y el resumen completo del modelo.
- **Residuales** — residuales en el tiempo, su ACF y su histograma, para verificar
  que el método capturó la estructura de la serie.
- **Comparación** — los tres métodos ajustados a la misma serie y una tabla de errores
  para contrastarlos de un vistazo.

### Requisitos

R 4.x y tres paquetes:

```r
install.packages(c("shiny", "forecast", "ggplot2"))
```

Para la serie de ejemplo Johnson & Johnson —ganancias trimestrales por acción,
1960–1980— hace falta además:

```r
install.packages("astsa")
```

Sin ella la app funciona igual, con `AirPassengers` como serie de ejemplo.

### Cómo ejecutarla

Desde la terminal, parado en la raíz del proyecto:

```bash
Rscript -e "shiny::runApp('src/suavizados.R', launch.browser=TRUE)"
```

En Windows con PowerShell, usa `Rscript.exe` (a secas, `R` es un alias de
`Invoke-History` y no funciona).

También puedes abrir `src/suavizados.R` en RStudio y usar el botón **Run App**.

### Formato del CSV

Dos columnas, con encabezado:

1. **Tiempo.** Acepta valores numéricos (`1960`, `1960.25`) o fechas (`1960-03-31`).
   Si no logra interpretar ninguno de los dos, usa el índice de fila.
2. **Observaciones.** Numéricas.

La frecuencia se deduce del salto entre tiempos consecutivos, pero se puede fijar a
mano (1 anual, 4 trimestral, 12 mensual) desde la barra lateral. Holt-Winters
estacional necesita frecuencia mayor que 1 y al menos dos ciclos completos; si la
serie no cumple, la app cae automáticamente a Holt (nivel y tendencia) y lo avisa.

### Notas sobre la implementación

**Alineación de la media móvil.** El suavizado en `t` promedia `y_t … y_{t-k+1}` y se
usa como predicción de `t+1`, así que los valores ajustados van corridos un período.
Sin ese corrimiento el error incluiría la observación que se intenta predecir y el
RMSE saldría artificialmente bajo. El pronóstico a varios pasos es plano: una media
móvil no extrapola tendencia.

**Intervalos de la media móvil.** Son aproximados y constantes (±1.96·s). SES y
Holt-Winters sí tienen intervalos propios derivados de su formulación en espacio de
estados, que se ensanchan con el horizonte; la media móvil no es un modelo
estadístico completo y no tiene una fórmula estándar equivalente.
