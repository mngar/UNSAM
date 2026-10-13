# ============================================================
# Script de instalación y carga de paquetes
# ============================================================

# 1. Definir la lista de paquetes requeridos
paquetes <- c(
  "simulMGF",
  "tidyverse",
  "nortest",
  "rMVP",
  "data.table",
  "CJAMP"
)

# 2. Identificar cuáles NO están instalados todavía
nuevos <- paquetes[!(paquetes %in% installed.packages()[, "Package"])]

# 3. Instalar únicamente los que faltan (con sus dependencias)
if (length(nuevos) > 0) {
  install.packages(nuevos, dependencies = TRUE)
} else {
  message("Todos los paquetes ya están instalados.")
}

# 4. Cargar todos los paquetes en la sesión actual
sapply(paquetes, function(pkg) {
  library(pkg, character.only = TRUE)
})
