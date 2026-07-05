# Scripts de detección de HFOs y análisis de solapamiento

Este documento resume los scripts desarrollados en Brainstorm/MATLAB para la detección automática de oscilaciones de alta frecuencia (HFOs) y el análisis de solapamiento con puntas interictales.  
Se presentan **ordenados de más simples a más complejos**, indicando su función principal, su utilidad práctica y cuál fue la **versión final** utilizada en el flujo de trabajo.

---

## 1. `process_evt_detect_hfo_candidates_by_channel.m`
### Descripción
Script básico de detección de HFOs aplicado a **un único canal seleccionado**.

### Qué hace
- Detecta candidatos a HFO en un solo canal.
- Trabaja sobre una señal previamente filtrada en **band-pass 60–250 Hz**.
- Tiene en cuenta los **segmentos malos**.
- Genera un **set de eventos** únicamente para ese canal.

### Utilidad
Fue útil como versión inicial para:
- comprobar que la detección funcionaba correctamente,
- inspeccionar manualmente la actividad en un canal concreto,
- validar parámetros de detección antes de ampliar el análisis al conjunto de canales.

### Limitaciones
- Solo analiza un canal cada vez.
- No permite estudiar la distribución espacial global de los HFOs.

---

## 2. `process_evt_detect_hfo_allch.m`
### Descripción
Versión extendida de la detección de HFOs aplicada a **todos los canales**, generando un set de eventos independiente para cada uno.

### Qué hace
- Detecta HFOs en los **72 canales** del EEG.
- Trabaja sobre señal filtrada en **band-pass 60–250 Hz**.
- Tiene en cuenta los **segmentos malos**.
- Genera un **set de eventos por canal**.

### Utilidad
Permite:
- estudiar la distribución canal a canal de los candidatos,
- comprobar en qué electrodos aparece la actividad de alta frecuencia,
- conservar la información espacial a nivel de canal.

### Limitaciones
- Genera muchos sets de eventos.
- Puede resultar menos práctico cuando se busca un análisis global de episodios HFO independientemente del canal.

---

## 3. `process_evt_detect_hfo_candidates.m`
### Descripción
Versión simplificada del análisis multicanal, en la que todos los HFOs detectados se guardan en **un único set de eventos**.

### Qué hace
- Detecta HFOs en todos los canales.
- Trabaja sobre señal filtrada en **band-pass 60–250 Hz**.
- Tiene en cuenta los **segmentos malos**.
- Agrupa todos los candidatos en **un solo set de eventos**.

### Utilidad
Facilita:
- disponer de un único grupo de HFOs candidatos,
- importar eventos más fácilmente,
- trabajar después con ellos en análisis de promediado o source imaging.

### Limitaciones
- Aunque agrupa todos los eventos en un solo set, sigue partiendo de una detección sobre señal ya filtrada.
- No resuelve automáticamente el flujo completo desde el **raw original**.

---

## 4. `process_evt_detect_hfos_without_bad.m`
### Descripción
Versión depurada de la detección global de HFOs, excluyendo explícitamente los candidatos detectados dentro de segmentos malos.

### Qué hace
- Detecta HFOs en todos los canales.
- Trabaja sobre señal filtrada en **band-pass 60–250 Hz**.
- **Excluye** los candidatos que coinciden con eventos `BAD`.
- Genera un único **set global de eventos**.

### Utilidad
Mejora la especificidad del pipeline al evitar:
- artefactos,
- ruido de movimiento,
- actividad espuria contenida en segmentos marcados como malos.

### Limitaciones
- Sigue requiriendo trabajar sobre una señal ya filtrada.
- El resultado todavía no se genera directamente sobre el raw original.

---

## 5. `process_evt_detect_hfos_original_raw_no_bads`
### Descripción
Versión más completa y robusta del detector de HFOs.  
Es la versión que integra el **preprocesado automático** y la **detección sobre el raw original**.

### Qué hace
- Parte del **raw original**.
- Aplica automáticamente:
  - **notch filter**
  - **band-pass 60–250 Hz**
- Detecta HFOs sobre la señal filtrada.
- Excluye candidatos en **segmentos malos**.
- Guarda los eventos resultantes **sobre el raw original**.
- Además, **fusiona candidatos cercanos** aunque provengan de canales distintos.

### Utilidad
Es la versión más práctica para el pipeline final, porque:
- automatiza el preprocesado,
- evita tener que importar manualmente los eventos desde un archivo filtrado al raw original,
- reduce redundancias cuando un mismo episodio HFO aparece simultáneamente en varios canales.

### Estado
✅ **Versión final del detector global de HFOs**

---

# Scripts de análisis de solapamiento (overlap)

---

## 6. `process_evt_detect_overlap_extended`
### Descripción
Primera versión del análisis de solapamiento entre dos grupos de eventos extendidos.

### Qué hace
- Compara dos sets de **eventos extendidos**.
- Evalúa si existe solapamiento temporal entre ellos según un **porcentaje mínimo definido por el usuario**.
- Genera un nuevo set que contiene **solo la parte combinada** entre ambos eventos.

### Utilidad
Fue útil para:
- detectar concordancia temporal entre dos tipos de eventos,
- estudiar coincidencias estrictas de actividad temporal.

### Limitaciones
- El evento resultante conserva solo la parte compartida.
- No siempre era la opción más adecuada cuando se quería conservar el evento spike completo.

---

## 7. `process_evt_detect_overlap_keep_spike`
### Descripción
Versión adaptada del análisis de overlap orientada al estudio de puntas con HFO.

### Qué hace
- Compara dos sets de **eventos extendidos**.
- Evalúa el solapamiento temporal según un **porcentaje mínimo**.
- Cuando existe overlap, genera un nuevo set cuyo **time-window coincide con el primer evento especificado** (por ejemplo, la spike completa).

### Utilidad
Fue diseñada específicamente para:
- identificar spikes que contienen HFOs,
- conservar la ventana temporal completa de la punta,
- generar el set `Spike/HFO` para su posterior importación y análisis.

### Estado
✅ **Versión final del análisis de overlap**

---

# Resumen de versiones finales

## Detector final de HFOs
**`process_evt_detect_hfos_original_raw_no_bads`**

Es la versión final porque:
- trabaja directamente desde el **raw original**,
- aplica el filtrado automáticamente,
- excluye segmentos malos,
- guarda los eventos sobre el raw original,
- y fusiona candidatos cercanos entre canales.

## Overlap final entre spikes y HFOs
**`process_evt_detect_overlap_keep_spike`**

Es la versión final porque:
- detecta el solapamiento temporal entre spikes y HFOs,
- y genera un nuevo set manteniendo la ventana temporal completa de la spike, lo cual es lo más útil para el análisis posterior.

---

# Orden de complejidad del desarrollo

1. `process_evt_detect_hfo_candidates_by_channel.m`
2. `process_evt_detect_hfo_allch.m`
3. `process_evt_detect_hfo_candidates.m`
4. `process_evt_detect_hfos_without_bad.m`
5. `process_evt_detect_hfos_original_raw_no_bads`
6. `process_evt_detect_overlap_extended`
7. `process_evt_detect_overlap_keep_spike`

---

# Uso recomendado en el pipeline final

Para el flujo de trabajo definitivo, se recomienda usar:

- **Detección de HFOs:**  
  `process_evt_detect_hfos_original_raw_no_bads`

- **Detección de spikes con HFOs:**  
  `process_evt_detect_overlap_keep_spike`

Estos dos scripts constituyen la base del pipeline final empleado para:
- detectar HFOs automáticamente,
- seleccionar spikes que contienen HFOs,
- generar eventos listos para importación,
- y continuar con el análisis mediante averaging y source imaging.
