# Exemplo: Múltiplos Blocos do Mesmo Tipo

## Caso de Uso: 2 Blocos KDE + 3 Blocos Composicionais

Este exemplo mostra como usar `prepare_data_blocks()` com múltiplos blocos de dados do mesmo tipo.

### Cenário

Você tem:
- **2 blocos KDE**: Zircão e Apatita (idades de minerais detríticos)
- **3 blocos Composicionais**: Óxidos maiores, Elementos traço, Mineralogia modal

### Dados de Exemplo

```r
# ============================================================================
# 1. Dados KDE - Zircão (idades U-Pb + temperatura Ti)
# ============================================================================
zircon_ages <- data.frame(
  sample = rep(c("S1", "S2", "S3", "S4", "S5"), each = 100),
  age_concordia = c(
    rnorm(100, 500, 50),   # S1
    rnorm(100, 1000, 100), # S2
    rnorm(100, 1500, 150), # S3
    rnorm(100, 800, 80),   # S4
    rnorm(100, 1200, 120)  # S5
  ),
  ti_temp = rnorm(500, 750, 50)
)

# ============================================================================
# 2. Dados KDE - Apatita (idades U-Pb + fission track)
# ============================================================================
apatite_ages <- data.frame(
  sample = rep(c("S1", "S2", "S3", "S4", "S5"), each = 80),
  age_u_pb = c(
    rnorm(80, 400, 40),
    rnorm(80, 600, 60),
    rnorm(80, 800, 80),
    rnorm(80, 500, 50),
    rnorm(80, 700, 70)
  ),
  age_fission_track = c(
    rnorm(80, 100, 10),
    rnorm(80, 150, 15),
    rnorm(80, 200, 20),
    rnorm(80, 120, 12),
    rnorm(80, 180, 18)
  )
)

# ============================================================================
# 3. Dados Composicionais - Óxidos Maiores (já em proporções)
# ============================================================================
major_oxides <- data.frame(
  sample = c("S1", "S2", "S3", "S4", "S5"),
  SiO2 = c(0.65, 0.60, 0.70, 0.62, 0.68),
  Al2O3 = c(0.15, 0.18, 0.12, 0.17, 0.14),
  FeO = c(0.10, 0.12, 0.08, 0.11, 0.09),
  MgO = c(0.05, 0.06, 0.04, 0.06, 0.05),
  CaO = c(0.05, 0.04, 0.06, 0.04, 0.04)
)

# ============================================================================
# 4. Dados Composicionais - Elementos Traço (já em proporções)
# ============================================================================
trace_elements <- data.frame(
  sample = c("S1", "S2", "S3", "S4", "S5"),
  Zr = c(0.30, 0.25, 0.35, 0.28, 0.32),
  Y = c(0.20, 0.22, 0.18, 0.21, 0.19),
  Nb = c(0.15, 0.18, 0.12, 0.17, 0.14),
  La = c(0.20, 0.20, 0.20, 0.19, 0.20),
  Ce = c(0.15, 0.15, 0.15, 0.15, 0.15)
)

# ============================================================================
# 5. Dados Composicionais - Mineralogia Modal (já em proporções)
# ============================================================================
modal_mineralogy <- data.frame(
  sample = c("S1", "S2", "S3", "S4", "S5"),
  Qtz = c(0.40, 0.35, 0.45, 0.38, 0.42),
  Fsp = c(0.30, 0.35, 0.25, 0.32, 0.28),
  Mica = c(0.15, 0.15, 0.15, 0.15, 0.15),
  Amph = c(0.10, 0.10, 0.10, 0.10, 0.10),
  Others = c(0.05, 0.05, 0.05, 0.05, 0.05)
)
```

### Preparação dos Blocos

```r
library(geomix)

# ============================================================================
# Preparar todos os 5 blocos de uma vez
# ============================================================================
prepared <- prepare_data_blocks(
  # 2 blocos KDE como lista nomeada
  kde_raw = list(
    Zircon = zircon_ages,
    Apatite = apatite_ages
  ),

  # Variáveis para cada bloco KDE
  kde_vars = list(
    Zircon = c("age_concordia", "ti_temp"),
    Apatite = c("age_u_pb", "age_fission_track")
  ),

  # 3 blocos composicionais como lista nomeada
  compositional_raw = list(
    MajorOxides = major_oxides,
    TraceElements = trace_elements,
    ModalMin = modal_mineralogy
  ),

  # Variáveis para cada bloco compositional
  compositional_vars = list(
    MajorOxides = c("SiO2", "Al2O3", "FeO", "MgO", "CaO"),
    TraceElements = c("Zr", "Y", "Nb", "La", "Ce"),
    ModalMin = c("Qtz", "Fsp", "Mica", "Amph", "Others")
  ),

  # Aplicar CLR nos composicionais
  apply_clr_compositional = TRUE,

  # Número de pontos para discretização KDE
  n_points_kde = 129,

  verbose = TRUE
)

# ============================================================================
# Verificar o resultado
# ============================================================================
names(prepared$data_list)
# [1] "Zircon"         "Apatite"        "MajorOxides"
# [4] "TraceElements"  "ModalMin"

prepared$data_types
#        Zircon        Apatite   MajorOxides TraceElements      ModalMin
#  "continuous"  "continuous"         "clr"         "clr"         "clr"

# Dimensões de cada bloco
sapply(prepared$data_list, dim)
#      Zircon Apatite MajorOxides TraceElements ModalMin
# [1,]      5       5           5             5        5
# [2,]    258     258           5             5        5
#
# Zircon: 5 amostras x 258 features (2 variáveis x 129 pontos)
# Apatite: 5 amostras x 258 features (2 variáveis x 129 pontos)
# MajorOxides: 5 amostras x 5 componentes (CLR)
# TraceElements: 5 amostras x 5 componentes (CLR)
# ModalMin: 5 amostras x 5 componentes (CLR)
```

### Unmixing

```r
# ============================================================================
# Realizar unmixing com todos os 5 blocos
# ============================================================================
result <- provenance_unmix(
  data_list = prepared$data_list,
  data_types = prepared$data_types,
  K = 3,  # 3 fontes (end-members)
  max_iter = 2000,
  verbose = TRUE
)

# ============================================================================
# Examinar resultados
# ============================================================================
print(result)
# Provenance Unmixing Results
# ===========================
#
# Data blocks: 5
#   - Zircon: continuous (258 features)
#   - Apatite: continuous (258 features)
#   - MajorOxides: clr (5 features)
#   - TraceElements: clr (5 features)
#   - ModalMin: clr (5 features)
#
# Number of sources (K): 3
# Number of samples (N): 5
# Iterations: 487
# Converged: TRUE
# Final loss: 1234.56

# Matriz de mistura A (proporções de cada fonte em cada amostra)
result$A
#     Source1 Source2 Source3
# S1    0.50    0.30    0.20
# S2    0.20    0.60    0.20
# S3    0.10    0.20    0.70
# S4    0.40    0.40    0.20
# S5    0.30    0.30    0.40

# Assinaturas das fontes para cada bloco
names(result$B_list)
# [1] "Zircon"         "Apatite"        "MajorOxides"
# [4] "TraceElements"  "ModalMin"

# Assinatura da Fonte 1 para Zircão
result$B_list$Zircon[1, 1:10]  # primeiros 10 valores

# Assinatura da Fonte 1 para Óxidos Maiores
result$B_list$MajorOxides[1, ]
#      SiO2     Al2O3       FeO       MgO       CaO
# 0.1234567 0.2345678 0.3456789 0.4567890 0.5678901
```

### Bootstrap para Incerteza

```r
# ============================================================================
# Estimar incerteza via bootstrap
# ============================================================================
boot_result <- bootstrap_provunmix(
  data_list = prepared$data_list,
  data_types = prepared$data_types,
  K = 3,
  n_boot = 100,
  verbose = TRUE
)

# Intervalos de confiança para a matriz de mistura
boot_result$A_mean    # médias bootstrap
boot_result$A_lower   # limite inferior 95%
boot_result$A_upper   # limite superior 95%
boot_result$A_se      # erro padrão

# Visualizar incerteza
plot(boot_result, type = "uncertainty")
```

### Visualizações

```r
# ============================================================================
# Gráficos
# ============================================================================

# 1. Proporções de mistura
plot(result, type = "mixing")

# 2. Assinaturas das fontes (escolher um bloco)
plot(result, type = "sources", block = 1)  # Zircon
plot(result, type = "sources", block = 3)  # MajorOxides

# 3. Histórico de convergência
plot(result, type = "loss")

# 4. Incerteza (se fez bootstrap)
plot(boot_result, type = "uncertainty")
```

## Alternativa: Aplicar CLR Seletivamente

Se você quiser aplicar CLR apenas em alguns blocos composicionais:

```r
prepared <- prepare_data_blocks(
  kde_raw = list(
    Zircon = zircon_ages,
    Apatite = apatite_ages
  ),
  kde_vars = list(
    Zircon = c("age_concordia", "ti_temp"),
    Apatite = c("age_u_pb", "age_fission_track")
  ),
  compositional_raw = list(
    MajorOxides = major_oxides,
    TraceElements = trace_elements,
    ModalMin = modal_mineralogy
  ),
  compositional_vars = list(
    MajorOxides = c("SiO2", "Al2O3", "FeO", "MgO", "CaO"),
    TraceElements = c("Zr", "Y", "Nb", "La", "Ce"),
    ModalMin = c("Qtz", "Fsp", "Mica", "Amph", "Others")
  ),
  # Aplicar CLR seletivamente (lista nomeada)
  apply_clr_compositional = list(
    MajorOxides = TRUE,
    TraceElements = TRUE,
    ModalMin = FALSE  # deixar em simplex
  )
)
```

## Flexibilidade Total

Você pode misturar qualquer combinação:

```r
# Exemplo: 1 KDE + 2 Point Counting + 3 Compositional
prepared <- prepare_data_blocks(
  kde_raw = zircon_ages,  # data.frame único
  point_counting_raw = list(
    Petrography = petrography_counts,
    HeavyMinerals = heavy_mineral_counts
  ),
  compositional_raw = list(
    MajorOxides = major_oxides,
    TraceElements = trace_elements,
    ModalMin = modal_mineralogy
  ),
  kde_vars = c("age_concordia", "ti_temp"),
  point_counting_vars = list(
    Petrography = c("Qtz", "Fsp", "Lith", "Musc"),
    HeavyMinerals = c("Zrn", "Ap", "Ttn", "Grt")
  ),
  compositional_vars = list(
    MajorOxides = c("SiO2", "Al2O3", "FeO"),
    TraceElements = c("Zr", "Y", "Nb"),
    ModalMin = c("Qtz", "Fsp", "Mica")
  )
)
```

## Resumo

**Sintaxe para múltiplos blocos:**
- Passe uma **lista nomeada** de data.frames
- Passe **listas nomeadas** de parâmetros (com mesmos nomes)
- Se parâmetro for valor único, aplica a todos os blocos
- Se parâmetro for lista nomeada, aplica valor específico a cada bloco

**Benefícios:**
- Interface limpa e intuitiva
- Suporta qualquer número de blocos de cada tipo
- Flexibilidade para parâmetros globais ou específicos
- Compatível com código existente (data.frame único ainda funciona)
