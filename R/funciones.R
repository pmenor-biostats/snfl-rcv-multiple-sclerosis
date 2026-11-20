# funciones.R -------------------------------------------------------------
# Funciones auxiliares del análisis de sNfL. Las usan analisis.Rmd y
# simular_datos.R.
# -------------------------------------------------------------------------

library(dplyr)
library(data.table)
library(pROC)

# Valor de referencia del cambio (RCV) ------------------------------------
# Versión log-normal (asimétrica): el aumento y la disminución significativos
# no son simétricos porque la sNfL se distribuye de forma log-normal.
#   RCV_subida = exp(Z * sqrt(2) * sigma) - 1
#   sigma = sqrt(ln(CVA^2 + 1) + ln(CVI^2 + 1))
# Con Z = 1,645 (unilateral, 95 %), CVA = 6,0 % y CVI = 8,64 % sale un 27,7 %.
rcv_increase <- function(cva, cvi, z = qnorm(0.95)) {
  sigma <- sqrt(log(cva^2 + 1) + log(cvi^2 + 1))
  exp(z * sqrt(2) * sigma) - 1
}

# Emparejamiento temporal de cada determinación con RM y brotes -----------
# Para cada determinación de sNfL:
#   - RM: la más cercana del mismo paciente (antes o después). Si está a más
#     de `window_days`, la determinación se clasifica como "RM inactiva".
#   - Brote: positivo si hay un brote confirmado a ±window_days.
match_events <- function(nfl, mri, relapses, window_days = 90) {
  nfl_dt <- as.data.table(nfl)
  mri_dt <- as.data.table(mri)[, mri_date_match := mri_date]

  # Unión "nearest" de data.table: la RM más cercana en el tiempo
  matched <- mri_dt[nfl_dt, on = .(patient_id, mri_date_match = nfl_date), roll = "nearest"]
  setnames(matched, "mri_date_match", "nfl_date")

  matched[, days_nfl_mri := as.integer(nfl_date - mri_date)]
  matched[, mri_within_w := as.integer(!is.na(days_nfl_mri) & abs(days_nfl_mri) <= window_days)]
  matched[, mri_active_w := as.integer(mri_within_w == 1L & mri_active %in% 1)]
  matched[, n_lesions_w  := fifelse(mri_active_w == 1L, n_lesions, NA_integer_)]

  # Brote a ±window_days de la extracción
  rel <- as.data.table(relapses)
  rel_w <- rel[matched[, .(patient_id, nfl_date)], on = "patient_id", allow.cartesian = TRUE,
               nomatch = NA][
    , .(relapse_w = as.integer(any(!is.na(relapse_date) &
                                     abs(as.integer(relapse_date - nfl_date)) <= window_days))),
    by = .(patient_id, nfl_date)]

  as.data.frame(rel_w[matched, on = .(patient_id, nfl_date)]) |>
    select(-any_of("i.nfl_date"))
}

# IC de Wilson para una proporción ---------------------------------------
wilson_ci <- function(x, n, conf = 0.95) {
  if (is.na(x) || is.na(n) || n == 0) return(c(NA_real_, NA_real_))
  z      <- qnorm(1 - (1 - conf) / 2)
  p      <- x / n
  denom  <- 1 + z^2 / n
  center <- (p + z^2 / (2 * n)) / denom
  margin <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denom
  c(max(0, center - margin), min(1, center + margin))
}

# Sensibilidad, especificidad, VPP y VPN con IC de Wilson -----------------
operating_point <- function(outcome, test_pos, label = "") {
  y <- as.integer(outcome); x <- as.integer(test_pos)
  TP <- sum(y == 1 & x == 1, na.rm = TRUE); TN <- sum(y == 0 & x == 0, na.rm = TRUE)
  FP <- sum(y == 0 & x == 1, na.rm = TRUE); FN <- sum(y == 1 & x == 0, na.rm = TRUE)
  fmt <- function(k, n) {
    ci <- wilson_ci(k, n)
    sprintf("%.1f (%.1f-%.1f)", 100 * k / n, 100 * ci[1], 100 * ci[2])
  }
  data.frame(label = label, TP, FP, TN, FN,
             sensitivity = fmt(TP, TP + FN), specificity = fmt(TN, TN + FP),
             ppv = fmt(TP, TP + FP), npv = fmt(TN, TN + FN))
}

# Punto de corte de Youden con IC bootstrap -------------------------------
# Remuestrea determinaciones con reemplazo y recalcula el punto de corte.
youden_boot <- function(outcome, predictor, boot_n = 2000, seed = 42) {
  fit <- function(y, x) {
    r <- roc(y, x, quiet = TRUE, direction = "<", levels = c(0, 1))
    coords(r, "best", best.method = "youden", transpose = FALSE,
           ret = c("threshold", "sensitivity", "specificity", "ppv", "npv"))[1, ]
  }
  observed <- fit(outcome, predictor)
  set.seed(seed)
  boot <- replicate(boot_n, {
    idx <- sample(length(outcome), replace = TRUE)
    if (length(unique(outcome[idx])) < 2) return(rep(NA_real_, 5))
    unlist(tryCatch(fit(outcome[idx], predictor[idx]), error = \(e) rep(NA_real_, 5)))
  })
  ci <- apply(boot, 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
  data.frame(metric = names(observed), estimate = unlist(observed),
             ci_low = ci[1, ], ci_high = ci[2, ], row.names = NULL)
}
