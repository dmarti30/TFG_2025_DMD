# Pipeline automatizado para detección de HFOs, concordancia temporal con puntas interictales y análisis mediante source imaging en epilepsia focal farmacorresistente

Repositorio del Trabajo de Fin de Grado centrado en el desarrollo y aplicación de un pipeline automatizado para:

- detectar **oscilaciones de alta frecuencia (HFOs)** en EEG de superficie,
- estudiar su **concordancia temporal con puntas interictales**,
- generar grupos experimentales para **source imaging**,
- y comparar espacialmente las reconstrucciones obtenidas respecto a la lesión anatómica.

---

## Objetivo del proyecto

Las oscilaciones de alta frecuencia (HFOs) son marcadores prometedores de la zona epileptógena, pero su detección en EEG de superficie es compleja debido a la presencia de actividad rápida no epileptiforme y artefactos.

Este proyecto desarrolla un flujo de trabajo que permite:

1. detectar automáticamente candidatos a HFO,
2. excluir segmentos malos (`BAD`),
3. identificar HFOs que co-ocurren con puntas interictales,
4. generar grupos experimentales comparables,
5. aplicar **source imaging** a cada grupo,
6. extraer scouts y comparar espacialmente su relación con la lesión.

---

## Estructura general del workflow

El flujo de trabajo completo va desde el **RAW** hasta los **resultados espaciales finales**.

### Etapas principales

1. **Preparación manual del RAW**
2. **Detección automática de HFOs**
3. **Análisis de overlap entre spikes y HFOs**
4. **Importación de epochs y averaging**
5. **Source imaging**
6. **Extracción manual de scouts**
7. **Análisis espacial y generación de resultados**
8. **Integración de resultados a nivel de cohorte**

---

## Componentes del proyecto

### 1. Scripts de detección de HFOs
Incluyen distintas versiones del detector, desde pruebas por canal hasta la versión final integrada sobre el RAW original.

La versión final recomendada es:

- `process_evt_detect_hfos_original_raw_no_bads`

Esta versión:
- parte del **raw original**,
- aplica automáticamente:
  - **notch filter**
  - **band-pass filter**
- detecta HFOs,
- excluye candidatos en segmentos `BAD`,
- guarda los eventos en el RAW original,
- y puede fusionar candidatos temporalmente cercanos entre canales.

---

### 2. Scripts de overlap entre spikes y HFOs
Permiten detectar la coincidencia temporal entre eventos extendidos.

La versión final recomendada es:

- `process_evt_detect_overlap_keep_spike`

Esta versión:
- compara spikes y HFOs,
- evalúa el solapamiento temporal según un porcentaje definido,
- y genera un nuevo set de eventos manteniendo la **ventana temporal completa de la spike**.

Este set se utiliza después como grupo experimental `Spike/HFO`.

---

### 3. Pipelines automáticos
El repositorio incluye scripts de MATLAB que llaman procesos de Brainstorm para automatizar parte del flujo de trabajo.

Entre ellos:

- `epilepsy_pipeline`
- `source_imaging_script`

Estos scripts permiten reducir la intervención manual en:
- detección,
- importación de eventos,
- averaging,
- y generación de reconstrucciones de fuente.

---

### 4. Scripts de resultados por paciente
Se desarrollaron distintos scripts para cuantificar la relación espacial entre:

- scout de la lesión (`mask`)
- scout de `Spike`
- scout de `Spike+HFO`
- scout de `HFO`

Estos scripts generan:
- tablas numéricas,
- figuras comparativas,
- archivos `.mat`,
- y reportes `.html`.

Ejemplos:
- `results_tfg_7.mlx`
- otras versiones previas de análisis espacial

---

### 5. Scripts de resultados globales de cohorte
Una vez generados los resultados individuales por paciente, se integran mediante scripts globales que:

- leen los HTML de todos los pacientes,
- extraen sus tablas,
- construyen una tabla común,
- calculan medias por grupo experimental,
- y permiten comparaciones a nivel de cohorte.

Ejemplo:
- `results_tfg_global_1.mlx`

---

## Grupos experimentales

A lo largo del proyecto se comparan principalmente tres grupos:

- **Spike**
- **Spike+HFO**
- **HFO**

Estos grupos se utilizan tanto en source imaging como en el análisis espacial posterior.

---

## Entorno de trabajo

Este proyecto combina:

- **Brainstorm** para:
  - visualización de RAW,
  - marcado de eventos,
  - importación de epochs,
  - averaging,
  - source imaging,
  - creación y exportación de scouts

- **MATLAB** para:
  - ejecución de procesos automatizados,
  - implementación de scripts personalizados,
  - análisis espacial,
  - generación de figuras y reportes

---

## Workflow resumido

### A. Preparación inicial
- comprobar frecuencia de muestreo del RAW
- duplicar y renombrar el grupo de spikes como `spikes`

### B. Detección automática
- ejecutar el pipeline de detección
- generar:
  - `HFO_candidate`
  - `Spike/HFO`

### C. Importación y averaging
- importar eventos como epochs
- generar los averages de:
  - `HFO`
  - `Spike/HFO`

### D. Source imaging
- aplicar source imaging a los grupos experimentales

### E. Extracción de scouts
- generar scouts manualmente a partir del máximo de actividad
- ajustar su extensión según la actividad visible

### F. Resultados
- comparar cada scout experimental con la lesión
- generar métricas espaciales y figuras
- guardar reportes individuales y globales

---

## Métricas espaciales utilizadas

Según la versión del script de resultados, se pueden calcular métricas como:

- número de vértices del scout
- área del scout (`cm²`)
- distancia mínima a la lesión
- distancia media a la lesión
- distancia máxima a la lesión
- mediana de distancias mínimas
- dispersión espacial
- desviación estándar de distancias

Estas métricas se calculan comparando los vértices de cada scout experimental con los vértices del scout de la lesión.

---

## Resultados generados

Dependiendo del script ejecutado, el repositorio puede generar:

- archivos `.mat`
- reportes `.html`
- figuras `.png`
- tablas resumen por paciente
- tablas globales de cohorte

---

## Requisitos

- MATLAB
- Brainstorm
- archivos RAW importados en Brainstorm
- eventos de spikes previamente marcados
- scouts exportados en formato `.mat`
- cortex exportado en formato `.mat`

---

## Uso recomendado

1. preparar el RAW en Brainstorm
2. ejecutar el pipeline de detección
3. generar los grupos experimentales
4. aplicar source imaging
5. extraer scouts
6. ejecutar el script de resultados del paciente
7. integrar todos los pacientes con el script global

---

## Organización sugerida del repositorio

```text
├── detection/
│   ├── process_evt_detect_hfo_candidates_by_channel.m
│   ├── process_evt_detect_hfo_allch.m
│   ├── process_evt_detect_hfo_candidates.m
│   ├── process_evt_detect_hfos_without_bad.m
│   └── process_evt_detect_hfos_original_raw_no_bads.m
│
├── overlap/
│   ├── process_evt_detect_overlap_extended.m
│   └── process_evt_detect_overlap_keep_spike.m
│
├── pipelines/
│   ├── epilepsy_pipeline.m
│   └── source_imaging_script.m
│
├── results/
│   ├── results_tfg_7.mlx
│   ├── results_tfg_global_1.mlx
│   └── otros scripts de resultados
│
├── docs/
│   ├── README de scripts específicos
│   └── documentación metodológica
│
└── README.md
