# results_tfg_7.mlx

Script de MATLAB para el análisis espacial de los scouts experimentales respecto a la lesión en pacientes con epilepsia focal farmacorresistente.

Este script compara tres grupos experimentales reconstruidos mediante source imaging:

- **Spike**
- **Spike+HFO**
- **HFO**

tomando como referencia anatómica el scout de la lesión (**mask**).

---

## Objetivo

`results_tfg_7.mlx` calcula métricas espaciales entre cada scout experimental y la lesión, y genera automáticamente:

- una **tabla de resultados**
- varias **figuras comparativas**
- un **archivo HTML** con los resultados visualizables
- un **archivo `.mat`** con las variables numéricas exportadas

---

## Archivos necesarios

Antes de ejecutar el script, debes disponer de los siguientes archivos exportados desde Brainstorm:

- **cortex** del paciente
- scout de la **lesión**: `mask`
- scout del grupo **Spike**: `spike`
- scout del grupo **Spike+HFO**: `spikeHFO`
- scout del grupo **HFO**: `HFO`

Todos estos archivos deben estar en formato `.mat`.

---

## Qué calcula el script

Para cada grupo experimental, el script calcula:

- `NumVertices`: número de vértices del scout
- `ScoutArea_cm2`: área del scout en cm²
- `MinDistance_mm`: media de las distancias mínimas de cada vértice del scout a la lesión
- `MeanDistance_mm`: media real de todas las distancias entre vértices del scout y vértices de la lesión
- `MaxDistance_mm`: media de las distancias máximas de cada vértice del scout a la lesión
- `MedianMinDistance_mm`: mediana de las distancias mínimas
- `SpatialDispersionMin_RMS_mm`: dispersión espacial calculada sobre las distancias mínimas
- `StdMinDistance_mm`: desviación estándar de las distancias mínimas

---

## Figuras generadas

El script genera figuras para comparar los tres grupos experimentales respecto a la lesión:

### Para `MinDistance_mm`
- gráfico tipo **raincloud** por grupo
- gráfico comparativo tipo **violin plot + boxplot + scatter**

### Para `MeanDistance_mm`
- gráfico tipo **raincloud** por grupo
- gráfico comparativo tipo **violin plot + boxplot + scatter**

### Para `MaxDistance_mm`
- gráfico tipo **raincloud** por grupo
- gráfico comparativo tipo **violin plot + boxplot + scatter**

Todas las figuras se generan usando la misma escala de distancia para facilitar la comparación visual entre métricas.

---

## Cómo usar el script

### 1. Exportar los scouts necesarios
Desde Brainstorm, exporta:

- `mask`
- `spike`
- `spikeHFO`
- `HFO`

y también el archivo del **cortex** del paciente.

---

### 2. Indicar el identificador del paciente
En el script, edita la variable:

```matlab
patient_id = 'FHP';
