# Workflow completo: desde RAW hasta Source Imaging + Scouts

Este documento explica, paso a paso, cómo obtener los grupos experimentales y sus reconstrucciones de fuente a partir de un archivo RAW de EEG.

El flujo combina:

- pasos **manuales** en Brainstorm
- pasos **automatizados** mediante scripts de MATLAB/Brainstorm

El objetivo final es obtener:

- los grupos experimentales:
  - **HFO**
  - **Spike/HFO**
- sus **averages**
- sus **source imaging**
- y los **scouts finales** para el análisis espacial

---

# Resumen del flujo

El proceso completo sigue este orden:

1. **Preparación manual del RAW**
2. **Aplicación del pipeline automático de detección**
3. **Copia manual de archivos necesarios para source imaging**
4. **Aplicación del pipeline automático de source imaging**
5. **Extracción manual de scouts finales**

---

# 1. Preparación inicial del RAW (manual)

Antes de ejecutar cualquier script, hay que preparar el archivo RAW en Brainstorm.

## 1.1 Comprobar la frecuencia de muestreo
Visualiza el RAW y comprueba si está registrado a:

- **256 Hz**
- **512 Hz**

Esto es importante porque los parámetros del band-pass cambian según la frecuencia del RAW.

## 1.2 Duplicar los eventos de spikes
Si ya tienes marcadas las poblaciones de spikes en el RAW:

- duplica ese grupo de eventos
- renómbralo como:

```text
spikes
