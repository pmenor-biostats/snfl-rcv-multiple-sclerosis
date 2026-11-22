# sNfL como biomarcador de actividad inflamatoria en la EMRR

**pMenor · Estadística clínica en R**

Código del análisis estadístico del artículo:

> de las Heras Flórez S, García de la Rosa G, Martín García V, Hernández Pérez
> MÁ, Oramas Padrón P, Burruezo Rodríguez C, Nuevo García J, Carretero Pérez M.
> Serum neurofilament light chain as a biomarker of inflammatory activity in
> relapsing-remitting multiple sclerosis: a real-world longitudinal cohort
> study. *Clin Chem Lab Med.* 2026.
> [doi:10.1515/cclm-2026-0944](https://doi.org/10.1515/cclm-2026-0944)

## La pregunta

La cadena ligera de neurofilamentos en suero (sNfL) sube cuando hay daño
axonal, por ejemplo en un brote de esclerosis múltiple. Pero varía mucho de un
paciente a otro y depende de la edad y de la función renal, así que un mismo
valor puede ser normal en una persona y alto en otra. ¿Es más útil comparar a
cada paciente con la población (umbral poblacional, **PBT**) o consigo mismo
(valor de referencia del cambio, **RCV**)?

## El análisis

967 determinaciones de sNfL de 295 pacientes con esclerosis múltiple
remitente-recurrente, emparejadas en el tiempo con RM y brotes clínicos
(ventana de ±90 días).

| Objetivo | Método |
|---|---|
| 1. ¿Cuánto sube la sNfL con la actividad? | Modelo lineal mixto sobre log(sNfL), con intercepto aleatorio por paciente y ajuste por edad, tiempo de evolución, sexo, tabaquismo, eGFR y diabetes |
| 2. ¿RCV o PBT se asocian mejor con la actividad? | Regresión logística con errores estándar robustos por paciente; acuerdo RCV-PBT con kappa |
| 3. ¿Cuál detecta mejor la actividad? | Curvas ROC, test de DeLong, sensibilidad, especificidad, VPP y VPN con IC de Wilson; punto de corte de Youden con IC bootstrap |
| 4. ¿Refleja la sNfL el número de lesiones? | Spearman, regresión log-log y puntos de corte por número de lesiones |

**Análisis de sensibilidad:** efecto de la edad sin actividad, eGFR
independiente de la edad, interacciones (edad y tabaquismo con la actividad),
solo determinaciones con RM dentro de la ventana, ventanas de ±60 y ±120 días
y ROC con una determinación por paciente.

Algunos detalles que pueden servir en otros estudios:

- **Emparejamiento temporal** de cada determinación con la RM más cercana
  (unión *nearest* de `data.table`) y con los brotes de la ventana, en
  [`R/funciones.R`](R/funciones.R).
- **RCV log-normal** (ver el apartado siguiente).
- **Datos repetidos por paciente** en todos los análisis: modelo mixto,
  errores robustos por paciente y ROC por paciente como comprobación.

### El valor de referencia del cambio (RCV): fórmula simétrica y log-normal

El RCV indica cuánto tiene que cambiar un resultado respecto al anterior del
mismo paciente para que el cambio no se explique solo por la variación
analítica (CVA) y biológica intraindividual (CVI). Hay dos formas de
calcularlo:

| | Fórmula | Con Z = 1,645 (unilateral, 95 %), CVA = 6,0 % y CVI = 8,64 % |
|---|---|---|
| **Simétrica** (clásica) | RCV = Z × √2 × √(CVA² + CVI²) | ±24,5 % |
| **Log-normal** (asimétrica) | RCV↑ = exp(Z × √2 × σ) − 1, con σ = √[ln(CVA² + 1) + ln(CVI² + 1)] | +27,7 % (subida) |

La simétrica supone que el resultado se distribuye de forma normal, así que la
subida y la bajada significativas son iguales. La sNfL, como muchas magnitudes
de laboratorio, se distribuye de forma log-normal: hace falta subir algo más
que bajar para que el cambio sea significativo. Por eso la log-normal es la
adecuada aquí, y es la que da el 27,7 % usado en el estudio.

Este repositorio usa la log-normal, en la función `rcv_increase()` de
[`R/funciones.R`](R/funciones.R).

### Resultados publicados

| Medida | Resultado |
|---|---|
| CCI ajustado | 0,536 |
| Brote clínico (±90 días) | +72,1 % de sNfL (IC 95 %: +53,7 % a +92,7 %) |
| RM activa (±90 días) | +24,5 % (IC 95 %: +13,7 % a +36,3 %) |
| eGFR G3-G5 frente a G1 | +67,8 % |
| RCV elevado → RM activa / brote | OR 4,93 / 6,48 |
| PBT elevado → RM activa / brote | OR 1,52 / 2,38 (no significativos) |
| Acuerdo RCV-PBT | kappa 0,457 |
| AUC para RM activa: sNfL / RCV / PBT | 0,707 / 0,669 / 0,649 |
| sNfL frente a n.º de lesiones nuevas | ρ = 0,498 |

**Conclusión del artículo:** la sNfL es un biomarcador válido de actividad
inflamatoria en la práctica real. Por la gran variabilidad entre pacientes,
interpretarla con el RCV (cada paciente frente a sí mismo) es más informativo
que con umbrales poblacionales, y la función renal debe tenerse siempre en
cuenta.

## Datos

Este repositorio **no contiene datos de pacientes**. `R/simular_datos.R` genera
tres tablas ficticias (determinaciones de sNfL, RM y brotes) con la misma
estructura y efectos parecidos a los del estudio. Las cifras se parecen a las
publicadas, pero no son las mismas.

El umbral poblacional (PBT) usa los límites superiores de referencia
(percentil 95) por décadas de edad establecidos con Lumipulse G en 282 sujetos
sanos: Varelli M, Fasano S, Varriale F, Lavano A. Establishment of normal
reference values for serum neurofilament light based on apparently healthy
subjects aged 18–90 years. *Neurol Sci.* 2024;45:37.

| Edad (años) | 18-29 | 30-39 | 40-49 | 50-59 | 60-69 | 70-90 |
|---|---|---|---|---|---|---|
| Límite superior (pg/mL) | 9,4 | 9,6 | 13,3 | 21,2 | 28,1 | 46,8 |

## Cómo reproducirlo

1. Clona el repositorio y abre la carpeta como proyecto de RStudio.
2. Instala los paquetes con las mismas versiones que se usaron (quedan
   registradas en `renv.lock`):

   ```r
   install.packages("renv")
   renv::restore()
   ```

3. Ejecuta `rmarkdown::render("analisis.Rmd")`.

El resultado, con código, salidas y figuras, se puede leer directamente en
[`analisis.md`](analisis.md). Las figuras usan las fuentes IBM Plex Sans e IBM
Plex Mono; si no están instaladas, R usa otras y avisa.

## Estructura

```
├── analisis.Rmd          análisis completo
├── analisis.md           el mismo análisis ya ejecutado (se lee en GitHub)
├── R/
│   ├── funciones.R       emparejamiento temporal, RCV, IC de Wilson, Youden bootstrap
│   ├── simular_datos.R   genera los datos ficticios
│   └── theme_pmenor.R    estilo de los gráficos
├── data/                 datos simulados (nfl.csv, mri.csv, relapses.csv)
└── output/               figuras
```

## Licencia

El código se publica con licencia [MIT](LICENSE): puedes usarlo y adaptarlo citando la fuente. Si lo usas en un trabajo, cita también el artículo.

---

**pMenor** · Estadística clínica en R · [pmenor.es](https://pmenor.es)
