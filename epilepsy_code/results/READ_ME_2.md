# results_tfg_global_1.mlx

Script de MATLAB para generar un **resumen global de cohorte** a partir de los resultados individuales obtenidos previamente con los scripts de análisis por paciente.

Este script no recalcula source imaging ni vuelve a analizar scouts directamente.  
Su función es **leer los HTML de resultados de múltiples pacientes**, extraer sus tablas numéricas y construir una tabla conjunta con los valores medios por grupo experimental.

---

## Objetivo

`results_tfg_global_1.mlx` permite integrar los resultados individuales de varios pacientes y obtener una visión global de la cohorte.

A partir de los archivos HTML generados previamente para cada paciente, el script:

- extrae la tabla numérica de resultados
- une todos los pacientes en una tabla común
- calcula la **media aritmética por grupo experimental**
- guarda un resumen general en formato **HTML** y **MAT**

---

## Qué archivos utiliza

Este script usa como entrada los archivos HTML de resultados individuales generados anteriormente para cada paciente.

Cada HTML debe contener la tabla de resultados con los grupos experimentales:

- **Spike**
- **Spike+HFO**
- **HFO**

y las variables numéricas asociadas.

El script busca automáticamente estos archivos dentro de la carpeta indicada.

---

## Qué calcula

A partir de todos los HTML encontrados, el script construye:

### 1. Tabla común de todos los pacientes
Incluye una fila por paciente y por grupo experimental.

### 2. Tabla media por grupo experimental
Calcula la **media aritmética** de las variables numéricas disponibles en la tabla de resultados.

Las variables que intenta promediar son:

- `NumVertices`
- `ScoutArea_cm2`
- `DistanceMaxPoints_mm`
- `MeanDistance_mm`
- `MedianDistance_mm`
- `SpatialDispersion_RMS_mm`
- `StdDistance_mm`

> El script solo promedia las variables que realmente existan en los HTML leídos.

---

## Cómo funciona

### 1. Indicar la carpeta raíz
En la variable `root_folder` debes indicar la carpeta donde están almacenados los resultados de los pacientes:

```matlab
root_folder = 'C:\Users\prego\Documents\datos_brainstorm\pacientes_rafa';
