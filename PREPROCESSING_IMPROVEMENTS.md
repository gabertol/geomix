# Melhorias no Sistema de Preprocessing

## Resumo das Mudanças

As funções de preprocessing foram refatoradas para separar claramente os **três tipos distintos de entrada de dados**:

### 1. KDE (Kernel Density Estimation)
- **Entrada**: Valores contínuos (ex: idades de zircão, um valor por grão)
- **Função**: `build_dz_kde_block()` (mantida sem mudanças)
- **Tipo de saída**: `"continuous"`

### 2. Point Counting (Contagens Discretas)
- **Entrada**: Contagens de grãos (ex: contei 100 grãos, 40 Qtz, 30 Fsp, 30 Lith)
- **Função**: `build_point_counting_block()` (NOVA)
- **Processamento**: converte contagens → proporções → opcionalmente CLR
- **Tipo de saída**: `"compositional"` ou `"clr"`

### 3. Compositional (Já Composições)
- **Entrada**: Dados já em forma de proporções (ex: 0.40 Qtz, 0.30 Fsp, 0.30 Lith, somam 1.0)
- **Função**: `build_compositional_block()` (NOVA)
- **Processamento**: valida closure → re-normaliza → opcionalmente CLR
- **Tipo de saída**: `"compositional"` ou `"clr"`

## Novas Funções

### `build_point_counting_block()`
```r
# Para contagens discretas (point counting)
pc_block <- build_point_counting_block(
  counts_df = data.frame(
    sample = c("S1", "S2", "S3"),
    Qtz = c(80, 70, 90),      # CONTAGENS (inteiros)
    Fsp = c(60, 55, 65),
    Lith = c(40, 35, 45)
  ),
  mineral_cols = c("Qtz", "Fsp", "Lith"),
  apply_clr = FALSE
)

# Retorna:
# - data_mat: matriz de proporções (ou CLR se apply_clr=TRUE)
# - samples: nomes das amostras
# - components: nomes dos componentes
# - total_counts: contagens totais por amostra
# - is_clr: lógico
# - data_type: "compositional" ou "clr"
```

### `build_compositional_block()`
```r
# Para dados já composicionais
comp_block <- build_compositional_block(
  comp_df = data.frame(
    sample = c("S1", "S2", "S3"),
    Qtz = c(0.40, 0.35, 0.45),    # PROPORÇÕES (decimais, somam 1)
    Fsp = c(0.30, 0.35, 0.30),
    Lith = c(0.30, 0.30, 0.25)
  ),
  comp_cols = c("Qtz", "Fsp", "Lith"),
  apply_clr = TRUE,
  check_closure = TRUE  # verifica se somam ~1.0
)

# Retorna:
# - data_mat: matriz de proporções (ou CLR se apply_clr=TRUE)
# - samples: nomes das amostras
# - components: nomes dos componentes
# - row_sums: somas das linhas (devem ser ~1.0)
# - is_clr: lógico
# - data_type: "compositional" ou "clr"
```

## Função Atualizada

### `prepare_data_blocks()`
Parâmetros renomeados para clareza:

**Novos parâmetros** (recomendados):
- `kde_raw` (antes: `dz_raw`)
- `point_counting_raw` (novo, para contagens)
- `compositional_raw` (novo, para composições)
- `kde_vars` (antes: `dz_vars`)
- `point_counting_vars` (antes: `bp_vars`)
- `compositional_vars` (novo)
- `n_points_kde` (antes: `n_points_dz`)
- `apply_clr_point_counting` (antes: `apply_clr_bp`)
- `apply_clr_compositional` (novo)

**Parâmetros deprecados** (ainda funcionam com warnings):
- `dz_raw`, `dz_vars`, `n_points_dz` → usar `kde_*`
- `bp_raw`, `bp_vars`, `apply_clr_bp` → usar `point_counting_*`
- `hm_raw`, `hm_vars`, `apply_clr_hm` → usar `point_counting_*` ou `compositional_*`

## Compatibilidade Retroativa

### `build_bp_block()` - DEPRECADA
Esta função agora está **deprecada** e emite um warning direcionando para:
- `build_point_counting_block()` se os dados são contagens
- `build_compositional_block()` se os dados já são composições

Ainda funciona internamente chamando `build_point_counting_block()`.

## Exemplos de Uso

### Exemplo 1: KDE + Point Counting
```r
# Dados de zircão (valores contínuos) + petrografia (contagens)
prepared <- prepare_data_blocks(
  kde_raw = my_zircon_ages,
  point_counting_raw = my_petrography_counts,
  kde_vars = c("age_concordia", "ti_temp"),
  point_counting_vars = c("Qtz", "Fsp", "Lith", "Musc"),
  n_points_kde = 129
)

result <- provenance_unmix(
  data_list = prepared$data_list,
  data_types = prepared$data_types,
  K = 3
)
```

### Exemplo 2: Compositional + Point Counting
```r
# Geoquímica (já proporções) + minerais pesados (contagens)
prepared <- prepare_data_blocks(
  compositional_raw = my_oxide_data,
  point_counting_raw = my_heavy_minerals,
  compositional_vars = c("SiO2", "Al2O3", "FeO", "MgO"),
  point_counting_vars = c("Zrn", "Ap", "Ttn", "Grt"),
  apply_clr_compositional = TRUE
)
```

### Exemplo 3: Uso direto das funções
```r
# Point counting
pc <- build_point_counting_block(
  counts_df = my_counts,
  mineral_cols = c("Qtz", "Fsp", "Lith")
)

# Compositional
comp <- build_compositional_block(
  comp_df = my_proportions,
  comp_cols = c("SiO2", "Al2O3"),
  check_closure = TRUE
)

# KDE
kde <- build_dz_kde_block(
  ages_df = my_ages,
  age_vars = "age_concordia"
)

# Unmix
result <- provenance_unmix(
  data_list = list(
    PC = pc$data_mat,
    Comp = comp$data_mat,
    KDE = kde$DZ_mat
  ),
  data_types = c(
    PC = pc$data_type,
    Comp = comp$data_type,
    KDE = "continuous"
  ),
  K = 3
)
```

## Diferenças Principais

### Point Counting vs Compositional

| Aspecto | Point Counting | Compositional |
|---------|---------------|---------------|
| **Entrada** | Contagens (inteiros: 80, 60, 40) | Proporções (decimais: 0.4, 0.3, 0.3) |
| **Exemplo** | "Contei 100 grãos" | "40% é Quartzo" |
| **Valores típicos** | 10-500 por categoria | 0.0-1.0, soma=1.0 |
| **Validação** | `min_count_per_sample` | `check_closure` |
| **Processamento** | divide por soma das linhas | verifica se já soma 1.0 |
| **Campo extra retornado** | `total_counts` | `row_sums` |

## Benefícios

1. **Clareza conceitual**: Funções com nomes que refletem o tipo de dado
2. **Validação apropriada**: Cada tipo tem suas próprias checagens
3. **Documentação explícita**: Exemplos claros de quando usar cada função
4. **Mensagens de erro melhores**: Validações específicas por tipo
5. **Compatibilidade retroativa**: Código antigo ainda funciona (com warnings)

## Migração

Para migrar código existente:

1. **Se usava `build_bp_block()` com CONTAGENS**:
   ```r
   # Antes
   bp <- build_bp_block(counts_df, mineral_cols, apply_clr)

   # Depois
   bp <- build_point_counting_block(counts_df, mineral_cols, apply_clr)
   ```

2. **Se usava `build_bp_block()` com PROPORÇÕES**:
   ```r
   # Antes
   bp <- build_bp_block(comp_df, mineral_cols, apply_clr)

   # Depois
   bp <- build_compositional_block(comp_df, mineral_cols, apply_clr)
   ```

3. **Em `prepare_data_blocks()`**:
   ```r
   # Antes
   prepare_data_blocks(dz_raw, bp_raw, dz_vars, bp_vars)

   # Depois (para contagens)
   prepare_data_blocks(kde_raw, point_counting_raw, kde_vars, point_counting_vars)

   # Depois (para composições)
   prepare_data_blocks(kde_raw, compositional_raw, kde_vars, compositional_vars)
   ```

## NAMESPACE

Funções exportadas adicionadas:
- `build_point_counting_block`
- `build_compositional_block`

Função mantida (deprecada):
- `build_bp_block`

## Múltiplos Blocos do Mesmo Tipo (NOVO!)

A função `prepare_data_blocks()` agora suporta **múltiplos blocos do mesmo tipo** através de listas nomeadas.

### Exemplo: 2 KDE + 3 Composicionais

```r
prepared <- prepare_data_blocks(
  # 2 blocos KDE
  kde_raw = list(
    Zircon = zircon_ages,
    Apatite = apatite_ages
  ),
  kde_vars = list(
    Zircon = c("age_concordia", "ti_temp"),
    Apatite = c("age_u_pb", "age_fission_track")
  ),

  # 3 blocos composicionais
  compositional_raw = list(
    MajorOxides = major_oxides,
    TraceElements = trace_elements,
    ModalMin = modal_mineralogy
  ),
  compositional_vars = list(
    MajorOxides = c("SiO2", "Al2O3", "FeO", "MgO"),
    TraceElements = c("Zr", "Y", "Nb", "La", "Ce"),
    ModalMin = c("Qtz", "Fsp", "Mica", "Amph")
  ),

  apply_clr_compositional = TRUE
)

# Resultado: 5 blocos preparados
names(prepared$data_list)
# [1] "Zircon"         "Apatite"        "MajorOxides"
# [4] "TraceElements"  "ModalMin"
```

### Como Funciona

**Sintaxe:**
- Para **um único bloco**: passe um `data.frame`
- Para **múltiplos blocos**: passe uma **lista nomeada** de `data.frames`

**Parâmetros:**
- Se parâmetro é **valor único**: aplica a todos os blocos
- Se parâmetro é **lista nomeada**: aplica valor específico a cada bloco

**Exemplos:**

```r
# Parâmetro global (aplica a todos)
apply_clr_compositional = TRUE

# Parâmetro específico por bloco
apply_clr_compositional = list(
  MajorOxides = TRUE,
  TraceElements = TRUE,
  ModalMin = FALSE
)

# Número de pontos KDE global
n_points_kde = 129

# Número de pontos KDE específico
n_points_kde = list(
  Zircon = 129,
  Apatite = 64
)
```

### Compatibilidade

- ✅ **Data.frame único** (comportamento original): funciona normalmente
- ✅ **Lista nomeada** (novo): processa múltiplos blocos
- ✅ **Misto**: pode combinar data.frame único para um tipo e lista para outro

```r
# Exemplo misto: 1 KDE + 3 Composicionais
prepared <- prepare_data_blocks(
  kde_raw = zircon_ages,  # data.frame único
  compositional_raw = list(  # lista nomeada
    MajorOxides = major_oxides,
    TraceElements = trace_elements,
    ModalMin = modal_mineralogy
  ),
  kde_vars = c("age_concordia", "ti_temp"),
  compositional_vars = list(
    MajorOxides = c("SiO2", "Al2O3"),
    TraceElements = c("Zr", "Y"),
    ModalMin = c("Qtz", "Fsp")
  )
)
```

Veja `MULTIPLE_BLOCKS_EXAMPLE.md` para exemplos completos e detalhados!
