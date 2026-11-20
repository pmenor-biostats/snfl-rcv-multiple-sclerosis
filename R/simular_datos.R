# simular_datos.R ---------------------------------------------------------
# Genera datos FICTICIOS con la misma estructura que los del estudio.
# Los datos reales no se publican: son de pacientes.
#
# Tres tablas, como en la práctica clínica:
#   data/nfl.csv       una fila por determinación de sNfL, con las covariables
#   data/mri.csv       una fila por RM: activa o no y número de lesiones nuevas
#   data/relapses.csv  una fila por brote clínico confirmado
#
# Qué reproduce (de forma aproximada):
#   - ~295 pacientes con EMRR y ~970 determinaciones (varias por paciente).
#   - ~9 % de determinaciones con RM activa y ~6 % con brote a ±90 días.
#   - Un "punto de partida" de sNfL propio de cada paciente (CCI ~0,5).
#   - Los efectos del artículo: brote ~+72 %, RM activa ~+25 % (más cuantas
#     más lesiones), edad +1 %/año, eGFR baja, exfumadores algo más bajos.
#   - Valores ausentes en tabaquismo, diabetes y eGFR.
# Las cifras no serán idénticas a las publicadas.
#
# Uso: source("R/simular_datos.R") crea los tres archivos en data/
# -------------------------------------------------------------------------

source("R/funciones.R")
set.seed(2026)

n_patients <- 295
study_start <- as.Date("2023-03-01")

# Pacientes -----------------------------------------------------------------
patients <- tibble(
  patient_id = sprintf("P%03d", seq_len(n_patients)),
  sex        = sample(c("Female", "Male"), n_patients, TRUE, c(0.78, 0.22)),
  age0       = pmin(pmax(rnorm(n_patients, 44, 10), 18), 72),
  dx_years0  = pmin(rgamma(n_patients, shape = 1.6, rate = 0.14), age0 - 12),
  smoking    = sample(c("Never", "Former", "Current", NA), n_patients, TRUE,
                      c(0.58, 0.17, 0.20, 0.05)),
  diabetes   = sample(c("No", "Yes", NA), n_patients, TRUE, c(0.93, 0.03, 0.04)),
  egfr0      = pmin(128 - 0.75 * age0 + rnorm(n_patients, 0, 13), 120),
  u          = rnorm(n_patients, 0, 0.30),   # nivel basal propio de cada paciente
  slope      = rnorm(n_patients, 0, 0.22),   # deriva lenta propia de cada paciente (log/año)
  start      = study_start + sample(0:420, n_patients, TRUE)
)

# Brotes ---------------------------------------------------------------------
relapses <- patients |>
  rowwise() |>
  reframe(patient_id, relapse_date = start + sort(sample(-60:700, rpois(1, 0.16))))

# RM de rutina (cada 9-14 meses) y RM por sospecha de brote -------------------
routine_mri <- patients |>
  rowwise() |>
  reframe(patient_id, mri_date = start - 150 + cumsum(sample(270:420, 3, TRUE)))
relapse_mri <- relapses |>
  filter(runif(n()) < 0.7) |>
  transmute(patient_id, mri_date = relapse_date + sample(0:40, n(), TRUE), after_relapse = TRUE)

mri <- bind_rows(routine_mri, relapse_mri) |>
  mutate(after_relapse = coalesce(after_relapse, FALSE),
         mri_active = as.integer(runif(n()) < ifelse(after_relapse, 0.75, 0.085)),
         n_lesions  = ifelse(mri_active == 1, pmin(1L + rnbinom(n(), size = 1.6, mu = 2.6), 20L), 0L),
         # El recuento de lesiones no consta en algunas RM activas
         n_lesions  = ifelse(mri_active == 1 & runif(n()) < 0.12, NA, n_lesions)) |>
  select(-after_relapse) |>
  arrange(patient_id, mri_date)

# Determinaciones de sNfL: de rutina (cada 4-8 meses) y ante sospecha de actividad
routine_nfl <- patients |>
  rowwise() |>
  reframe(patient_id,
          nfl_date = start + cumsum(c(0, sample(120:240, min(rpois(1, 2.2), 7), TRUE))))
event_nfl <- bind_rows(
  relapses |> filter(runif(n()) < 0.6) |>
    transmute(patient_id, nfl_date = relapse_date + sample(5:45, n(), TRUE)),
  mri |> filter(mri_active == 1, runif(n()) < 0.35) |>
    transmute(patient_id, nfl_date = mri_date + sample(-30:30, n(), TRUE))
)
nfl <- bind_rows(routine_nfl, event_nfl) |>
  distinct(patient_id, nfl_date) |>
  arrange(patient_id, nfl_date)

# Emparejar con los eventos para generar la sNfL de cada determinación
matched <- match_events(nfl, mri, relapses, window_days = 90) |>
  left_join(patients, by = "patient_id") |>
  mutate(
    years_from_start = as.numeric(nfl_date - start) / 365.25,
    age            = round(age0 + years_from_start, 1),
    years_since_dx = round(dx_years0 + years_from_start, 1),
    egfr           = round(egfr0 - 0.8 * years_from_start + rnorm(n(), 0, 4)),
    egfr_g = case_when(egfr >= 90 ~ 0, egfr >= 60 ~ 0.15, TRUE ~ 0.52),
    mri_effect = ifelse(mri_active_w == 1,
                        0.66 * log(coalesce(n_lesions_w, 3L) + 1) - 0.70, 0),
    # Durante la actividad, los exfumadores suben menos
    active = mri_active_w == 1 | relapse_w == 1,
    smoke_effect = case_when(smoking %in% "Former" ~ -0.10 - 0.30 * active,
                             smoking %in% "Current" ~ -0.03, TRUE ~ 0),
    log_nfl = log(10.4) + u + 0.011 * (age - 45) + egfr_g + smoke_effect +
      slope * years_from_start + 0.50 * relapse_w + mri_effect + rnorm(n(), 0, 0.13),
    nfl = round(exp(log_nfl), 2),
    # eGFR "> 90" informado sin valor exacto; algún valor ausente
    egfr = ifelse(runif(n()) < 0.03, NA, egfr)
  )

nfl_out <- matched |>
  select(patient_id, nfl_date, nfl, age, years_since_dx, sex, smoking, egfr, diabetes) |>
  arrange(patient_id, nfl_date)

dir.create("data", showWarnings = FALSE)
write.csv(nfl_out, "data/nfl.csv", row.names = FALSE, na = "")
write.csv(mri, "data/mri.csv", row.names = FALSE, na = "")
write.csv(relapses, "data/relapses.csv", row.names = FALSE, na = "")
