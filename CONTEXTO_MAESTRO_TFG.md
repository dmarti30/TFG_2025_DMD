# CONTEXTO MAESTRO DEL PROYECTO TFG

## Guía de continuidad y replicación para ChatGPT, MATLAB y Brainstorm

> **Propósito**
>
> Este documento concentra el contexto científico, metodológico, técnico y operativo necesario para continuar este proyecto sin perder las decisiones acumuladas durante su desarrollo.
>
> Debe utilizarse como documento de referencia antes de modificar un proceso de Brainstorm, un pipeline de MATLAB, un script de resultados o la interpretación de las métricas.
>
> La prioridad es **reproducibilidad**. Un nuevo chat debe considerar que está continuando un proyecto existente y no rediseñándolo desde cero. Si se propone una mejora, debe distinguirse de la versión actual y no sustituir silenciosamente la lógica ya validada.

---

# 1. Resumen del proyecto

El proyecto es un Trabajo de Fin de Grado centrado en epilepsia focal farmacorresistente y en el análisis de EEG de superficie mediante MATLAB + Brainstorm.

Título de trabajo utilizado:

**Precise detection of the seizure onset zone in drug-resistant epilepsy using an EEG pipeline based on high-frequency oscillations correlated with interictal discharges.**

La hipótesis general es que las **oscilaciones de alta frecuencia (HFOs)** pueden aportar información útil para la localización de regiones epileptógenas y que los HFOs que presentan concordancia temporal con **descargas interictales / spikes / IEDs** pueden ser especialmente informativos.

El flujo general desarrollado es:

~~~text
EEG RAW
  ↓
preprocesado
  ↓
detección automática de HFOs
  ↓
concordancia temporal HFO ↔ Spike
  ↓
generación de grupos experimentales
  ↓
epoching
  ↓
average
  ↓
source imaging
  ↓
scouts corticales
  ↓
comparación espacial con lesión MRI
  ↓
métricas por paciente
  ↓
análisis de cohorte
~~~

La evolución futura principal consiste en automatizar los pasos todavía manuales y mejorar la especificidad del detector.

---

# 2. Contexto científico y terminología

## 2.1 Cohorte

En la fase de cohorte se ha trabajado con aproximadamente **40 pacientes** con epilepsia focal farmacorresistente.

Características generales del conjunto de datos:

- EEG de superficie de alta densidad.
- Aproximadamente 72 canales.
- Spikes / IEDs marcadas o revisadas por un epileptólogo.
- MRI estructural.
- Lesión visible utilizada como referencia anatómica cuando está disponible.
- Brainstorm como plataforma principal de EEG, epochs, averages, source imaging y scouts.
- MATLAB como entorno de automatización, procesos personalizados y análisis de resultados.

La lesión de MRI se utiliza como **proxy anatómico** de la región epileptógena. No debe describirse como una ground truth perfecta de la seizure onset zone.

## 2.2 Grupos experimentales

Los tres grupos conceptuales que deben mantenerse son:

1. **Spike**
2. **Spike+HFO** / **Spike_HFO** / **Spike/HFO** según el contexto de archivo o evento
3. **HFO**

Interpretación:

- **Spike**: puntas interictales de referencia.
- **Spike+HFO**: spikes que cumplen el criterio de concordancia temporal con un HFO.
- **HFO**: candidatos HFO detectados independientemente de la coincidencia con spike.

Los nombres concretos de evento y carpeta pueden variar por limitaciones de Brainstorm, pero no debe cambiarse su significado.

---

# 3. Repositorio y estructura actual

Repositorio principal:

~~~text
dmarti30/TFG_2025_DMD
~~~

Estructura relevante:

~~~text
TFG_2025_DMD/
├── README.md
├── process_evt_detect_hfos_original_raw_no_bads.m
├── process_evt_detect_overlap_keep_spike.m
└── epilepsy_code/
    ├── cristina_code/
    ├── pipelines/
    │   ├── READ_ME.md
    │   ├── READ_ME.txt
    │   ├── epilepsy_pipeline.m
    │   └── source_imaing_script.m
    ├── results/
    │   ├── READ_ME.md
    │   ├── READ_ME_2.md
    │   ├── result_tfg.mlx
    │   ├── result_tfg_2.mlx
    │   ├── result_tfg_3.mlx
    │   ├── result_tfg_4.mlx
    │   ├── result_tfg_5.mlx
    │   ├── results_tfg_6.mlx
    │   ├── results_tfg_7.mlx
    │   └── results_tfg_global_1.mlx
    └── scripts/
        ├── READ_ME.md
        ├── process_evt_detect_hfo_candidates_by_channel.m
        ├── process_evt_detect_hfo_allch.m
        ├── process_evt_detect_hfo_candidates.m
        ├── process_evt_detect_hfos_without_bad.m
        ├── process_evt_detect_hfos_original_raw_no_bads.m
        ├── process_evt_detect_overlap_extended.m
        └── process_evt_detect_overlap_keep_spike.m
~~~

**Importante:** el archivo del segundo pipeline está actualmente escrito como:

~~~text
source_imaing_script.m
~~~

con ese typo en el nombre. No asumir automáticamente que se llama `source_imaging_script.m` al buscarlo en el repositorio.

También existen copias de los dos procesos finales en la raíz y dentro de `epilepsy_code/scripts/`. Al modificar código, comprobar qué copia es la que realmente utiliza la instalación local de Brainstorm y mantener sincronización si se desea conservar ambas.

---

# 4. Diferencia fundamental: proceso Brainstorm vs script MATLAB

Esta distinción es crítica.

## 4.1 Procesos personalizados de Brainstorm

Son archivos MATLAB que siguen la arquitectura de Brainstorm y aparecen en su interfaz.

Nombres típicos:

~~~text
process_evt_detect_...
process_create_...
~~~

Arquitectura habitual:

~~~matlab
function varargout = process_xxx(varargin)
    eval(macro_method);
end

function sProcess = GetDescription()
    ...
end

function Comment = FormatComment(sProcess)
    ...
end

function OutputFiles = Run(sProcess, sInputs)
    ...
end
~~~

Se integran con:

- Process1 / Process2.
- `bst_process('CallProcess', ...)`.
- `bst_report`.
- las estructuras internas de la base de datos de Brainstorm.

La ubicación típica para procesos personalizados del usuario es:

~~~text
$HOME/.brainstorm/process/
~~~

Cuando el usuario pide una **función para Brainstorm**, **integrada en la interfaz**, o que aparezca como un proceso seleccionable, la solución correcta suele ser un proceso Brainstorm con esta arquitectura.

## 4.2 Scripts MATLAB de pipeline

Ejemplos:

~~~text
epilepsy_pipeline.m
source_imaing_script.m
~~~

Son scripts secuenciales generados o adaptados a partir del pipeline de Brainstorm. Llaman a procesos existentes mediante `bst_process`.

Un pipeline no debe reimplementar innecesariamente la lógica interna de un proceso custom.

### Regla de mantenimiento

Si se modifica un pipeline:

- conservar el orden lógico salvo petición expresa;
- no alterar parámetros no solicitados;
- no duplicar la lógica del detector dentro del pipeline;
- comprobar que el output de cada paso es compatible con el siguiente.

---

# 5. APIs y estructuras de Brainstorm relevantes

Al desarrollar código nuevo, preferir las APIs de Brainstorm frente a manipulación arbitraria de archivos MAT.

Funciones relevantes:

~~~matlab
in_bst(...)
in_bst_data(...)
in_bst_results(...)
in_bst_channel(...)
in_tess_bst(...)
in_fread(...)
bst_get(...)
bst_memory(...)
bst_save(...)
db_template(...)
db_add_data(...)
bst_process(...)
bst_report(...)
tess_vertconn(...)
panel_scout(...)
panel_record(...)
process_notch(...)
process_bandpass(...)
~~~

## 5.1 Estructura de un archivo de sources

Campos relevantes:

~~~matlab
ResultsMat.ImageGridAmp
ResultsMat.ImagingKernel
ResultsMat.Time
ResultsMat.SurfaceFile
ResultsMat.HeadModelType
ResultsMat.nComponents
ResultsMat.DataFile
ResultsMat.Atlas
ResultsMat.GridAtlas
ResultsMat.GoodChannel
ResultsMat.ZScore
~~~

Un resultado puede estar:

- almacenado completamente en `ImageGridAmp`, o
- almacenado como kernel `ImagingKernel` + datos asociados.

Por ello, para reconstruir el mapa completo cuando sea necesario, una opción robusta es:

~~~matlab
ResultsMat = in_bst(ResultFile, TimeWindow, 1);
~~~

No asumir que `ImageGridAmp` estará siempre materializado en disco.

## 5.2 Scout Brainstorm

Crear la estructura con:

~~~matlab
sScout = db_template('Scout');
~~~

Campos principales:

~~~matlab
sScout.Vertices
sScout.Seed
sScout.Color
sScout.Label
sScout.Function
sScout.Region
sScout.Handles
~~~

Significado:

- `Vertices`: índices de los vértices que forman el scout.
- `Seed`: vértice de referencia/origen.
- `Label`: nombre.
- `Function`: función de scout.
- `Region`: codificación de región de Brainstorm.

### Regla crítica sobre Seed

En un scout funcional creado desde el máximo, el `Seed` puede representar el vértice de máxima actividad si se construyó de esa forma.

En el scout anatómico de lesión `mask`, el `Seed` **NO es un máximo de actividad**. La lesión no contiene actividad eléctrica. Es únicamente un vértice de referencia del scout.

Por tanto, nunca describir una distancia:

~~~text
group_seed ↔ lesion_seed
~~~

como “distancia entre máximos”.

---

# 6. Preparación manual del RAW

Antes del pipeline:

1. visualizar el RAW;
2. comprobar la frecuencia de muestreo;
3. localizar la población de spikes;
4. duplicarla;
5. renombrar la copia exactamente como:

~~~text
spikes
~~~

## 6.1 Frecuencia de muestreo y banda HFO

Si el RAW está aproximadamente a:

~~~text
Fs = 256 Hz
~~~

Nyquist es 128 Hz, por lo que se ha planteado:

~~~text
high-pass = 60 Hz
low-pass  = 120 Hz
~~~

Si:

~~~text
Fs = 512 Hz
~~~

se ha utilizado:

~~~text
high-pass = 60 Hz
low-pass  = 250 Hz
~~~

**Nunca usar un low-pass de 250 Hz con Fs=256 Hz.**

---

# 7. Evolución de los detectores HFO

Orden de desarrollo de más simple a más completo:

~~~text
1. process_evt_detect_hfo_candidates_by_channel.m
2. process_evt_detect_hfo_allch.m
3. process_evt_detect_hfo_candidates.m
4. process_evt_detect_hfos_without_bad.m
5. process_evt_detect_hfos_original_raw_no_bads.m
~~~

## 7.1 `process_evt_detect_hfo_candidates_by_channel.m`

Primera aproximación.

- Un único canal.
- Detección sobre señal filtrada.
- Banda de trabajo histórica 60–250 Hz.
- Consideración de segmentos malos.
- Útil para validar visualmente el detector y ajustar parámetros.

No es la versión final.

## 7.2 `process_evt_detect_hfo_allch.m`

Segunda evolución.

- Detecta HFOs en todos los canales.
- Aproximadamente 72 canales.
- Genera un set de eventos independiente por canal.
- Útil para inspección espacial canal a canal.

Problema: genera demasiados sets para un workflow global.

## 7.3 `process_evt_detect_hfo_candidates.m`

Tercera evolución.

- Todos los canales.
- Un único set global de HFO.
- Más práctico para importación posterior.
- Sigue trabajando sobre un archivo ya filtrado.

## 7.4 `process_evt_detect_hfos_without_bad.m`

- Todos los canales.
- Un set global.
- Excluye candidatos detectados dentro de segmentos BAD.
- Sigue dependiendo de un band-pass previamente disponible.

## 7.5 `process_evt_detect_hfos_original_raw_no_bads.m`

**Versión final recomendada del detector.**

Parte del RAW original y realiza:

~~~text
RAW original
  ↓
notch en memoria
  ↓
band-pass en memoria
  ↓
detección HFO
  ↓
rechazo BAD
  ↓
merge intra-canal
  ↓
merge global entre canales
  ↓
eventos HFO guardados sobre RAW original
~~~

Ventajas:

- no obliga a trabajar manualmente con raws filtrados;
- mantiene los timestamps directamente sobre el raw original;
- evita duplicar un mismo episodio cuando aparece en varios canales;
- conserva notas con información de canal.

El propio código contiene una regla explícita: la lógica de detección se preservó desde la versión anterior y la principal evolución fue el preprocesado interno y el merge global multicanal.

---

# 8. Lógica del detector HFO final

La secuencia del detector actual es aproximadamente:

1. seleccionar canales EEG;
2. leer RAW continuo;
3. calcular Fs a partir del vector temporal;
4. aplicar ventana temporal opcional;
5. notch;
6. band-pass;
7. por cada canal:
   - Hilbert;
   - envolvente;
   - threshold robusto;
   - histéresis ON/OFF;
   - merge de segmentos cercanos;
   - duración mínima;
   - zero crossings;
   - ciclos mínimos;
   - máximo de zero crossings dinámico;
   - rechazo si solapa BAD;
8. recoger candidatos de todos los canales;
9. merge global multicanal;
10. guardar evento extendido en RAW original.

## 8.1 Threshold robusto

La lógica implementada es:

~~~matlab
env = abs(hilbert(x));

m = median(env);
madv = median(abs(env - m));
sigmaRob = 1.4826 * madv;

thr_on  = m + kMad * sigmaRob;
thr_off = m + (kOffRatio * kMad) * sigmaRob;
~~~

Si `sigmaRob` es prácticamente cero, el código puede recurrir a `std(env)`.

## 8.2 Duración y ciclos

La detección filtra candidatos por:

- duración mínima;
- número mínimo de ciclos;
- número máximo plausible de zero crossings.

Conceptualmente:

~~~matlab
zc = sum(diff(seg > 0) ~= 0);
cycles = zc / 2;
~~~

y el máximo dinámico:

~~~matlab
maxZC = floor(2 * fHighAss * dur * maxZCFact);
~~~

## 8.3 Segmentos BAD

Los BAD se leen del **RAW original**.

Un candidato que solape un intervalo BAD debe rechazarse.

Motivación:

- movimiento;
- EMG;
- saturación;
- ruido técnico;
- artefactos evidentes.

## 8.4 Merge multicanal

El detector puede generar varios candidatos para el mismo episodio fisiológico porque aparece en varios canales.

La versión final aplica un merge global de candidatos temporalmente solapados o suficientemente próximos y combina la información de canal.

La finalidad es:

> contar episodios temporales, no multiplicar el mismo HFO por el número de electrodos en que aparece.

---

# 9. Parámetros actuales del detector: distinguir defaults y pipeline

Este punto es importante porque el proceso tiene **defaults de interfaz**, mientras que `epilepsy_pipeline.m` puede sobrescribirlos.

## 9.1 Defaults definidos dentro del proceso

En la versión del repositorio:

~~~text
eventname       = HFO_candidate
kMad            = 5
kOffRatio       = 0.5
min duration    = 20 ms (valor mostrado por el proceso)
merge gap       = 10 ms (valor mostrado por el proceso)
min cycles      = 4
fHighAss        = 240 Hz
maxZCFactor     = 1.3
overwrite       = true
notchfreqs      = [50 100 150 200] Hz
highpass        = 60 Hz
lowpass         = 250 Hz
keepfiltered    = false
~~~

## 9.2 Valores actualmente escritos en `epilepsy_pipeline.m`

La llamada del pipeline del repositorio usa:

~~~matlab
'eventname',    'HFO_candidate'
'timewindow',   []
'kmad',         6
'koff',         0.5
'mindur',       0.040
'mergegap',     0.005
'mincycles',    4
'fhigh',        240
'maxzcfactor',  1
'overwrite',    1
'notchfreqs',   50
'highpass',     60
'lowpass',      250
'keepfiltered', 0
~~~

### ADVERTENCIA IMPORTANTE DE UNIDADES

El proceso custom y el script generado por Brainstorm deben revisarse conjuntamente antes de interpretar `0.040` y `0.005` como 40 ms y 5 ms.

En el código del proceso se leen las opciones y se convierten usando divisiones por 1000. Brainstorm también puede aplicar transformaciones de unidades en opciones `value`. Por tanto:

> **No modificar ni reinterpretar silenciosamente estos valores. Antes de una nueva ejecución de cohorte, verificar en Brainstorm qué valor efectivo llega a `Run` y documentarlo.**

Esta comprobación es especialmente importante para asegurar reproducibilidad.

## 9.3 Regla de conservación

Si el usuario pide, por ejemplo, “modifica solo el merge”, no cambiar:

- `kMad`;
- `kOffRatio`;
- duración;
- filtro;
- ciclos;
- zero crossings;
- BAD;
- Hilbert.

Cambiar exclusivamente la parte solicitada.

---

# 10. Procesos de overlap HFO–Spike

## 10.1 `process_evt_detect_overlap_extended.m`

Primera versión.

- Compara dos grupos de eventos extendidos.
- Usa porcentaje mínimo de overlap.
- Genera un nuevo evento basado en la región combinada/solapada según esa primera implementación.

No es la versión final recomendada.

## 10.2 `process_evt_detect_overlap_keep_spike.m`

**Versión final recomendada.**

Objetivo: determinar qué spikes contienen suficiente overlap con un HFO sin perder la ventana completa de la spike.

Criterio implementado:

~~~matlab
overlapDur = min(s2, h2) - max(s1, h1);

minRequired = minOverlapFrac * min(spikeDur, hfoDur);

if overlapDur >= minRequired
    keepSpike(iSpk) = true;
end
~~~

Es decir:

~~~text
overlap_duration >= porcentaje * duración_del_evento_más_corto
~~~

Si se cumple, la salida no es solo la intersección: se conserva el intervalo completo de la spike.

Ejemplo:

~~~text
Spike:    |-------------------------|
HFO:               |------|

Salida:   |-------------------------|
~~~

Esto es importante para que `Spike+HFO` mantenga una ventana comparable con `Spike` en epoching, average y source imaging.

Defaults del proceso:

~~~text
spikelabel   = IED
hfolabel     = HFO
newlabel     = Spike+HFO
minoverlap   = 33 %
overwrite    = true
~~~

Pero el pipeline actual lo llama con:

~~~text
spikelabel = spikes
hfolabel   = HFO_candidate
newlabel   = Spike/HFO
minoverlap = 10 %
overwrite  = false
~~~

La versión realmente usada por el pipeline debe describirse mediante estos parámetros de llamada, no únicamente mediante los defaults de la interfaz.

---

# 11. `epilepsy_pipeline.m`

Primer pipeline global.

## 11.1 Preparación manual

Antes de ejecutarlo:

1. editar la ruta al raw original;
2. indicar el nombre del paciente;
3. adaptar la banda a Fs;
4. asegurar que existe el evento `spikes`.

## 11.2 Secuencia exacta del pipeline actual

### A. Detectar HFO

Llama:

~~~text
process_evt_detect_hfos_original_raw_no_bads
~~~

y genera:

~~~text
HFO_candidate
~~~

### B. Convertir spikes a eventos extendidos

Llama a:

~~~text
process_evt_extended
~~~

con:

~~~text
eventname  = spikes
timewindow = [-0.20, 0]
~~~

Por tanto, cada spike se extiende hacia atrás respecto al marcador según esta configuración.

### C. Overlap

Llama:

~~~text
process_evt_detect_overlap_keep_spike
~~~

con:

~~~text
Spike label  = spikes
HFO label    = HFO_candidate
Output       = Spike/HFO
Min overlap  = 10 %
~~~

### D. Convertir a eventos simples

Llama:

~~~text
process_evt_simple
~~~

para:

~~~text
spikes
Spike/HFO
HFO_candidate
~~~

con método:

~~~text
end
~~~

Es decir, conserva el **final del evento** al convertir a evento simple.

### E. Importar epochs

Para `Spike/HFO`:

~~~text
epoch = [-0.20, +0.20] s
~~~

Para `HFO_candidate`:

~~~text
epoch = [-0.20, +0.20] s
~~~

Configuración actual relevante:

~~~text
split       = 0
createcond  = 1
ignoreshort = 1
usectfcomp  = 0
usessp      = 0
baseline    = []
~~~

### F. Average

Se utiliza:

~~~text
Average: By trial group (folder average)
Arithmetic mean
No weighting
~~~

para `Spike/HFO` y para `HFO_candidate`.

---

# 12. Paso manual intermedio antes de source imaging

Tras generar epochs/averages, históricamente se copian manualmente entre carpetas:

~~~text
Head model
Noise covariance
Data covariance
~~~

En Brainstorm:

~~~text
Click derecho → Copy to other folders
~~~

Este paso garantiza que las condiciones HFO y Spike/HFO tengan disponibles los objetos requeridos por el modelo inverso.

Es una de las fases futuras candidatas a automatización.

---

# 13. `source_imaing_script.m`

Segundo pipeline.

**Nombre actual del archivo en el repositorio:** `source_imaing_script.m`.

La versión actual hace dos cosas diferentes:

## 13.1 Input Spike

Carga un resultado de source imaging de Spike ya existente y aplica:

~~~text
Unconstrained to flat map
method = norm
~~~

Es decir, convierte orientaciones XYZ a:

~~~matlab
sqrt(x^2 + y^2 + z^2)
~~~

cuando corresponde.

## 13.2 Input HFO

Carga el average HFO y calcula sources mediante:

~~~text
process_inverse_2018
InverseMethod  = lcmv
InverseMeasure = nai
SourceOrient   = fixed
DataTypes      = EEG
output         = kernel only / shared
~~~

Parámetros presentes en el script:

~~~text
Loose          = 0.2
UseDepth       = 1
WeightExp      = 0.5
WeightLimit    = 10
NoiseMethod    = median
NoiseReg       = 0.1
SnrMethod      = rms
SnrRms         = 1e-06
SnrFixed       = 3
ComputeKernel  = 1
~~~

## 13.3 Input Spike/HFO

Aplica la misma configuración LCMV/NAI sobre el average de `Spike_HFO`.

### Regla

No asumir que este método es universalmente el mejor source imaging. Es la configuración utilizada por esta versión del proyecto.

Históricamente también se exploraron métodos como:

- sLORETA;
- cMEM;
- wMEM;
- beamformer;
- dipole modelling.

Si se cambia el método inverso, no se deben mezclar resultados nuevos con los antiguos sin recalcular de forma homogénea la cohorte.

---

# 14. Generación manual de scouts

Procedimiento original:

1. abrir el mapa de sources;
2. localizar el máximo de actividad;
3. crear un scout con ese máximo como origen;
4. crecer/reducir el scout hasta adaptarlo a la actividad visible;
5. repetir para:
   - Spike;
   - Spike+HFO;
   - HFO.

La lesión se conserva como scout anatómico `mask`.

Problemas del método manual:

- subjetividad;
- variabilidad entre pacientes;
- tiempo;
- difícil reproducibilidad;
- dependencia del operador.

Por ello, una evolución prioritaria es automatizar esta etapa.

---

# 15. Proceso automático de scouts

Se diseñó una primera versión:

~~~text
process_create_scout_from_sources_auto.m
~~~

Su objetivo es producir automáticamente un scout Brainstorm a partir del mapa cortical reconstruido.

## 15.1 Lógica propuesta

1. cargar el source map completo;
2. seleccionar una ventana temporal;
3. obtener una magnitud no negativa por vértice;
4. encontrar el máximo global;
5. usar ese vértice como `Seed`;
6. aplicar threshold relativo;
7. conservar solo el componente cortical conectado que contiene el seed;
8. imponer opcionalmente tamaño mínimo/máximo;
9. crear `db_template('Scout')`;
10. añadir/reemplazar el scout en un atlas;
11. guardar el atlas en la superficie cortical;
12. registrar QC.

## 15.2 Modo recomendado: mapa en el instante del máximo global

Buscar:

~~~matlab
[PeakValue, LinearIndex] = max(VertexActivity(:));
[SeedVertex, iPeakTime] = ind2sub(size(VertexActivity), LinearIndex);
SpatialScore = VertexActivity(:, iPeakTime);
~~~

Esto reproduce mejor el procedimiento manual de localizar primero el máximo y después delimitar la activación visible alrededor de ese instante.

## 15.3 Threshold relativo

Valor inicial de experimentación:

~~~text
70 %
~~~

Cálculo:

~~~matlab
AbsoluteThreshold = 0.70 * SpatialScore(SeedVertex);
CandidateMask = SpatialScore >= AbsoluteThreshold;
~~~

**70 % no es un umbral clínicamente validado.**

Debe hacerse sensibilidad, por ejemplo:

~~~text
50 %
60 %
70 %
80 %
90 %
~~~

comparando contra scouts manuales.

Después se fija un único criterio para toda la cohorte.

## 15.4 Componente conectado

No incluir todos los vértices suprathreshold del cerebro.

Conservar solo el componente conectado al máximo.

Conectividad cortical:

~~~matlab
VertConn = tess_vertconn(Vertices, Faces);
~~~

Esto evita que un segundo foco distante por encima del threshold entre artificialmente en el mismo scout.

## 15.5 Fuentes constrained vs unconstrained

Si:

~~~text
nComponents = 1
~~~

usar:

~~~matlab
VertexActivity = abs(ImageGridAmp);
~~~

Si:

~~~text
nComponents = 3
~~~

calcular norma XYZ:

~~~matlab
sqrt(X.^2 + Y.^2 + Z.^2)
~~~

No comparar componentes XYZ por separado con resultados constrained.

## 15.6 Inputs que deben rechazarse o manejarse explícitamente

- volume source models;
- mixed models;
- resultados ya reducidos a atlas;
- archivos sin `SurfaceFile`;
- dimensiones incompatibles;
- ventana temporal vacía;
- scout vacío.

---

# 16. Validación del scout automático

Antes de sustituir la delimitación manual:

1. escoger pacientes representativos;
2. generar scout manual;
3. generar scout automático;
4. repetir para distintos thresholds;
5. comparar:
   - Dice;
   - Jaccard;
   - área;
   - número de vértices;
   - seed;
   - distancia a lesión;
6. seleccionar criterio;
7. congelarlo para la cohorte.

Una validación más avanzada puede incluir bootstrap de eventos y mapas de estabilidad.

---

# 17. `results_tfg_7.mlx`: análisis por paciente

Script individual actual.

Inputs:

~~~text
cortex del paciente
mask
spike
spikeHFO
HFO
~~~

Los scouts se exportan desde Brainstorm como MAT.

El usuario configura:

~~~matlab
patient_id
cortex_file
mask_file
spike_file
spikeHFO_file
HFO_file
~~~

y ejecuta el Live Script.

---

# 18. Coordenadas y unidades

Brainstorm almacena normalmente los vértices corticales en metros.

Para trabajar en mm:

~~~matlab
unit_scale = 1000;
VerticesMM = VerticesRaw * unit_scale;
~~~

Nunca mezclar coordenadas en metros con distancias etiquetadas como milímetros.

---

# 19. Matriz de distancias

Para scout experimental:

~~~matlab
group_coords_mm
~~~

y lesión:

~~~matlab
lesion_coords_mm
~~~

calcular:

~~~matlab
D = pdist2(group_coords_mm, lesion_coords_mm);
~~~

Si el grupo tiene N vértices y la lesión M:

~~~text
D es N × M
~~~

Ejemplo 4 × 4:

~~~text
16 distancias
~~~

Toda la nomenclatura actual de `MinDistance`, `MeanDistance` y `MaxDistance` nace de esta matriz.

---

# 20. Definiciones EXACTAS de las métricas espaciales actuales

Este apartado es crítico para no volver a mezclar las versiones históricas.

## 20.1 `MinDistance_mm`

Anteriormente esta operación se llamaba erróneamente `MeanDistance_mm`.

Ahora:

~~~matlab
dmin = min(D, [], 2);
MinDistance_mm = mean(dmin);
~~~

Interpretación:

1. cada vértice experimental se compara con todos los vértices de la lesión;
2. se conserva la distancia más corta para ese vértice experimental;
3. se promedian esas distancias mínimas.

No es:

- la distancia mínima absoluta entre scouts;
- la media de todas las distancias.

## 20.2 `MeanDistance_mm`

Nueva definición:

~~~matlab
dmean = mean(D, 2);
MeanDistance_mm = mean(dmean);
~~~

Equivalente a:

~~~matlab
MeanDistance_mm = mean(D(:));
~~~

Interpretación:

> media real de **todas** las distancias vértice experimental ↔ vértice lesión.

En un ejemplo 4 × 4 se usan las 16 distancias.

## 20.3 `MaxDistance_mm`

Sustituye al antiguo `DistanceMaxPoints_mm`.

~~~matlab
dmax = max(D, [], 2);
MaxDistance_mm = mean(dmax);
~~~

Interpretación:

1. para cada vértice experimental se busca el vértice de lesión más lejano;
2. se obtiene un máximo por vértice experimental;
3. se hace la media de esos máximos.

No usar seed-to-seed para esta variable.

---

# 21. Métricas adicionales

El script individual calcula o puede calcular:

~~~text
NumVertices
ScoutArea_cm2
MinDistance_mm
MeanDistance_mm
MaxDistance_mm
MedianMinDistance_mm
SpatialDispersionMin_RMS_mm
StdMinDistance_mm
~~~

Definiciones:

~~~matlab
MedianMinDistance_mm = median(dmin);

SpatialDispersionMin_RMS_mm = sqrt(mean(dmin.^2));

StdMinDistance_mm = std(dmin);
~~~

Estas tres últimas se basan en las distancias mínimas por vértice.

---

# 22. Área del scout

Para cada triángulo cortical completamente incluido:

~~~matlab
tri_area = 0.5 * norm(cross(v2-v1, v3-v1));
~~~

Se suman todas las caras.

Si las coordenadas están en metros:

~~~text
m² → cm² = × 10,000
~~~

Una formulación usada:

~~~matlab
ScoutArea_cm2 = total_area_raw * (unit_scale^2) / 100;
~~~

con:

~~~matlab
unit_scale = 1000;
~~~

---

# 23. Figuras de resultados

Para cada familia:

~~~text
MinDistance
MeanDistance
MaxDistance
~~~

generar las mismas figuras:

1. raincloud por grupo;
2. comparación conjunta con:
   - violin;
   - boxplot;
   - scatter.

Grupos:

~~~text
Spike
Spike+HFO
HFO
~~~

## Regla de escala

Las tres métricas deben compartir una escala global comparable.

Construir el rango con todos los valores de:

~~~text
min_vertex_distances
mean_vertex_distances
max_vertex_distances
~~~

y aplicar límites coherentes entre figuras.

---

# 24. Output de `results_tfg_7.mlx`

Se genera una carpeta similar a:

~~~text
results_tfg_7_patient_<PATIENT_ID>/
│
├── results_tfg_7_patient_<PATIENT_ID>.html
├── results_tfg_7_patient_<PATIENT_ID>.mat
└── figures/
~~~

El HTML contiene tabla + figuras.

El MAT debe conservar las variables numéricas y distribuciones para permitir análisis global posterior sin necesidad de parsear HTML.

---

# 25. Error histórico: `DistanceMaxPoints_mm`

Una versión anterior calculaba:

~~~text
DistanceMaxPoints_mm
~~~

como distancia:

~~~text
group_seed ↔ lesion_seed
~~~

Esto es conceptualmente problemático porque el seed de lesión no es un máximo funcional.

Debe considerarse obsoleto para las nuevas métricas.

No volver a usarlo como “distancia entre máximos”.

---

# 26. Si se necesita “máximo funcional → lesión”

La métrica correcta es:

~~~matlab
source_peak_coord = VerticesMM(source_seed, :);

d = sqrt(sum((lesion_coords_mm - source_peak_coord).^2, 2));

PeakToLesionMinDistance_mm = min(d);
~~~

Interpretación:

> distancia desde el máximo funcional al punto más cercano de la máscara lesional.

No confundir con:

~~~matlab
min(dmin)
~~~

porque `min(dmin)` puede partir de cualquier vértice del scout, no necesariamente del máximo funcional.

---

# 27. Resultados históricos de cohorte y cambio de nomenclatura

Un resultado de cohorte discutido previamente fue aproximadamente:

~~~text
Spike       ≈ 30.9 mm
Spike+HFO   ≈ 28.5 mm
HFO         ≈ 66.9 mm
~~~

con:

~~~text
Friedman p < 0.001
Kendall W ≈ 0.4
~~~

Post hoc Wilcoxon + Holm:

~~~text
Spike vs Spike+HFO  ≈ p 0.4
Spike vs HFO        < 0.001
Spike+HFO vs HFO    < 0.001
~~~

y ANOVA de medidas repetidas aproximadamente:

~~~text
F ≈ 40
p < 0.001
~~~

### ADVERTENCIA

Esos valores provenían de la variable histórica llamada `MeanDistance_mm`, cuya operación era:

~~~matlab
mean(min(D, [], 2))
~~~

En la nomenclatura actual eso corresponde a:

~~~text
MinDistance_mm
~~~

No etiquetar esos resultados antiguos como el nuevo `MeanDistance_mm = mean(D(:))`.

---

# 28. `results_tfg_global_1.mlx`

Script de resumen global.

La versión actual/histórica:

1. busca HTML individuales;
2. extrae tablas;
3. añade paciente;
4. concatena;
5. calcula medias por grupo;
6. genera HTML y MAT globales.

Outputs típicos:

~~~text
general_results_summary/
├── results_general_mean.html
└── results_general_mean.mat
~~~

Variables:

~~~text
CommonTable
MeanTable
html_files
~~~

## Problema de compatibilidad a vigilar

El README histórico del global todavía enumera columnas antiguas:

~~~text
DistanceMaxPoints_mm
MeanDistance_mm
MedianDistance_mm
SpatialDispersion_RMS_mm
StdDistance_mm
~~~

mientras `results_tfg_7` actual usa:

~~~text
MinDistance_mm
MeanDistance_mm
MaxDistance_mm
MedianMinDistance_mm
SpatialDispersionMin_RMS_mm
StdMinDistance_mm
~~~

Por tanto, antes de usar el script global con los nuevos resultados hay que actualizarlo o verificar su compatibilidad.

### Mejora recomendada

Para futuras versiones, leer directamente los MAT individuales es más robusto que parsear HTML.

---

# 29. Estadística de cohorte

La unidad estadística debe ser **el paciente**, no cada vértice del scout.

Los vértices son útiles para visualizar la distribución intra-paciente, pero no deben tratarse como sujetos independientes.

Para tres condiciones emparejadas:

~~~text
Spike
Spike+HFO
HFO
~~~

se ha utilizado:

- Friedman;
- Kendall's W;
- Wilcoxon signed-rank post hoc;
- corrección de Holm.

ANOVA de medidas repetidas se ha utilizado como análisis adicional.

---

# 30. Resultados piloto históricos

En una fase inicial se trabajó con un paciente piloto.

Datos históricos:

~~~text
72 canales
referencia Cz
duración ≈ 20,785 s ≈ 5.77 h
70 spikes manuales
8,573 candidatos HFO globales
2,113 candidatos en canales próximos a lesión
8 HFO con overlap global
8/70 ≈ 11.4 %
3 overlaps en el conjunto restringido de canales
~~~

Estos datos son del piloto y no deben confundirse con la cohorte final.

Canales usados en aquella restricción exploratoria:

~~~text
FT9, T9, TP9, P9, T3, TP7, FT7, FC6, T5,
C5, CP5, P5, FC3, C3, CP3, P3, FC1, C1,
CP1, P1
~~~

No convertir esta lista en una regla general para todos los pacientes.

---

# 31. Workflow completo actual

## Fase A — manual

~~~text
1. Abrir RAW.
2. Confirmar Fs.
3. Duplicar spikes.
4. Renombrar como "spikes".
~~~

## Fase B — `epilepsy_pipeline.m`

~~~text
RAW
 ↓
notch
 ↓
band-pass
 ↓
HFO_candidate
 ↓
BAD rejection
 ↓
merge multicanal
 ↓
spikes extendidos
 ↓
overlap → Spike/HFO
 ↓
eventos simples
 ↓
epochs
 ↓
averages HFO y Spike/HFO
~~~

## Fase C — manual

Copiar:

~~~text
head model
noise covariance
data covariance
~~~

## Fase D — `source_imaing_script.m`

~~~text
Spike source preexistente
HFO average
Spike/HFO average
  ↓
LCMV / NAI
  ↓
source maps
~~~

## Fase E — scouts

Históricamente:

~~~text
máximo → scout manual → grow/shrink
~~~

Evolución deseada:

~~~text
source map → auto-scout
~~~

## Fase F — resultados

Exportar:

~~~text
cortex
mask
spike
spikeHFO
HFO
~~~

Ejecutar:

~~~text
results_tfg_7.mlx
~~~

## Fase G — cohorte

Integrar resultados individuales y ejecutar estadística emparejada.

---

# 32. Validación del detector HFO

Que un script ejecute sin error no implica que la detección sea válida.

## 32.1 Validación técnica

Comprobar:

- timestamps dentro del raw;
- inicio < fin;
- duración plausible;
- canales válidos;
- ausencia de candidatos en BAD;
- merge correcto;
- no duplicados no deseados;
- eventos guardados en raw original;
- ausencia de desplazamiento temporal.

## 32.2 Validación visual

Revisar una muestra de candidatos mostrando:

- raw;
- señal filtrada;
- morfología;
- número de ciclos;
- relación con spike;
- posible artefacto.

## 32.3 Validación clínica

Idealmente comparar con:

- revisión de epileptólogo;
- SEEG/SOZ cuando exista;
- resección;
- outcome postquirúrgico;
- lesión MRI;
- hipótesis clínica multidisciplinar.

La lesión MRI es una referencia anatómica, no una definición perfecta de epileptogenicidad.

---

# 33. Mejora futura de especificidad

Una evolución razonable es separar:

~~~text
detector sensible
   ↓
candidatos
   ↓
filtro de falsos positivos
   ↓
HFO de alta confianza
~~~

Variables potenciales:

- duración;
- potencia;
- amplitud relativa;
- line length;
- número de canales;
- distribución espacial;
- actividad excesivamente global;
- proximidad entre canales;
- espectro;
- mapa tiempo-frecuencia;
- sospecha EMG;
- broadband abrupto;
- número de ciclos.

No sustituir el detector actual sin validación comparativa.

---

# 34. Score de confianza futuro

Puede crearse un score que combine:

~~~text
duración
ciclos
potencia HFO
amplitud robusta
número de canales
consistencia espacial
overlap con spike
posición del HFO dentro de la spike
probabilidad de artefacto
~~~

Esto es una extensión futura, no parte de la definición actual del pipeline.

---

# 35. Bootstrap y estabilidad espacial

Para evaluar si una localización es robusta:

1. seleccionar subconjuntos de eventos;
2. repetir average;
3. repetir source imaging;
4. generar scout automático;
5. medir:
   - overlap entre scouts;
   - desplazamiento del máximo;
   - distancia a lesión;
6. repetir múltiples veces.

Outputs potenciales:

~~~text
media
SD
IC 95 %
mapa de probabilidad espacial
~~~

Esto ayuda a distinguir una reconstrucción estable de una altamente dependiente de los eventos incluidos.

---

# 36. Principios obligatorios al escribir o modificar código

## 36.1 No cambiar lógica no solicitada

Si el usuario pide:

> “cambia únicamente el merge”

no modificar:

- threshold;
- filtro;
- Hilbert;
- duración;
- ciclos;
- criterio BAD.

## 36.2 Entregar scripts completos cuando se pide una versión funcional

Un proceso “perfecto” debe incluir:

- archivo completo;
- `GetDescription`;
- `FormatComment`;
- `Run`;
- validación de inputs;
- mensajes de error;
- comentarios;
- outputs claros.

## 36.3 Mantener nombres importantes

Ejemplos:

~~~text
spikes
HFO_candidate
Spike/HFO
Spike
Spike+HFO
HFO
mask
~~~

No renombrarlos por estética si un pipeline depende de ellos.

## 36.4 Validar dimensiones de sources

Siempre comprobar:

~~~matlab
nVertices = size(SurfaceMat.Vertices, 1);
~~~

y la correspondencia con:

~~~text
nRowsSources = nVertices × nComponents
~~~

No asumir que todas las fuentes son constrained.

## 36.5 Usar `bst_report`

En procesos Brainstorm:

~~~matlab
bst_report('Info', ...)
bst_report('Warning', ...)
bst_report('Error', ...)
~~~

para que los problemas aparezcan en el report.

## 36.6 No inventar APIs Brainstorm

Si una llamada interna no se conoce con seguridad:

1. buscar un proceso oficial similar;
2. revisar la firma;
3. replicar su patrón.

---

# 37. Cómo crear correctamente un nuevo proceso Brainstorm

Cuando el usuario pida un proceso integrado:

1. buscar un proceso oficial comparable;
2. imitar su estructura;
3. definir:
   - `Comment`;
   - `Category`;
   - `SubGroup`;
   - `Index`;
   - `InputTypes`;
   - `OutputTypes`;
   - `nInputs`;
   - `nMinFiles`;
4. declarar opciones de interfaz;
5. cargar datos con loaders oficiales;
6. validar tipos;
7. ejecutar;
8. guardar mediante mecanismos Brainstorm;
9. registrar en report;
10. probar con un paciente antes de lanzar cohorte.

---

# 38. Matemáticamente correcto no significa Brainstorm-compatible

Un script puede fallar aunque su algoritmo sea correcto por:

- rutas relativas;
- results tipo link/kernel;
- cache de superficies;
- atlas activo;
- `SurfaceFile`;
- `nComponents`;
- raw vs data vs results;
- registro en database;
- eventos simples vs extendidos;
- unidades transformadas por opciones de GUI.

Validar siempre el nivel matemático y el nivel Brainstorm.

---

# 39. Distancia euclídea actual vs posibles métricas futuras

Las métricas actuales usan distancia euclídea 3D:

~~~matlab
pdist2(group_coords_mm, lesion_coords_mm)
~~~

No son distancias geodésicas sobre la superficie cortical.

Si en el futuro se añade distancia geodésica:

- crear una métrica nueva;
- no reemplazar silenciosamente la euclídea;
- recalcular toda la cohorte si se desea comparabilidad.

---

# 40. Razones metodológicas de decisiones actuales

## Detectar sobre señal filtrada y guardar sobre RAW original

Permite:

- sensibilidad a HFO;
- mantener timestamps;
- revisar sobre raw;
- evitar múltiples archivos derivados en el workflow clínico.

## Excluir BAD

Reduce falsos positivos obvios.

## Merge multicanal

Evita contar un episodio tantas veces como canales lo detectan.

## Conservar spike completa en Spike/HFO

Mantiene ventanas comparables para average y source imaging.

## Scout conectado al máximo

Evita incluir focos suprathreshold aislados.

## Umbral fijo de scout

Reduce subjetividad entre pacientes una vez validado.

## Lesión como máscara anatómica

Impide atribuir “actividad” a un objeto MRI que no la tiene.

---

# 41. Errores que otro chat NO debe repetir

No debe:

- llamar al seed de lesión “máximo de actividad”;
- usar seed-to-seed como `MaxDistance_mm`;
- llamar `MeanDistance` a `mean(min(D,[],2))` en la nomenclatura actual;
- cambiar lógica del detector cuando se solicita otra modificación;
- confundir proceso Brainstorm con script MATLAB;
- usar 250 Hz con Fs=256 Hz;
- incluir todos los clusters suprathreshold en un scout automático sin conectividad;
- tratar vértices como sujetos independientes en estadística de cohorte;
- cambiar nombres de eventos sin revisar dependencias;
- asumir que un script es clínicamente válido porque ejecuta;
- afirmar que lesión MRI = SOZ;
- mezclar métricas de versiones antiguas y nuevas bajo el mismo nombre;
- asumir que los defaults del proceso son los parámetros realmente usados por el pipeline.

---

# 42. Checklist para detector HFO

~~~text
[ ] RAW correcto
[ ] Fs verificada
[ ] filtro compatible con Nyquist
[ ] opciones efectivas verificadas
[ ] canales EEG correctos
[ ] BAD cargados
[ ] Hilbert/threshold sin cambios no solicitados
[ ] duración válida
[ ] ciclos mínimos válidos
[ ] max ZC válido
[ ] candidatos BAD rechazados
[ ] merge intra-canal correcto
[ ] merge multicanal correcto
[ ] timestamps preservados
[ ] eventos en RAW original
[ ] canales/notas conservados
~~~

---

# 43. Checklist para overlap

~~~text
[ ] spike label correcto
[ ] HFO label correcto
[ ] ambos eventos extendidos
[ ] porcentaje interpretado sobre evento más corto
[ ] salida conserva spike completa
[ ] duplicados eliminados
[ ] overwrite/combinación revisado
[ ] output label compatible con pipeline
~~~

---

# 44. Checklist para scout automático

~~~text
[ ] source map cortical
[ ] SurfaceFile válido
[ ] nComponents correcto
[ ] ventana temporal correcta
[ ] mapa completo reconstruido si era kernel
[ ] máximo global correcto
[ ] Seed = máximo funcional
[ ] threshold documentado
[ ] componente conectado al seed
[ ] scout no vacío
[ ] atlas correcto
[ ] superficie guardada
[ ] inspección visual de QC
~~~

---

# 45. Checklist para `results_tfg_7`

~~~text
[ ] cortex correcto
[ ] mask correcto
[ ] Spike correcto
[ ] SpikeHFO correcto
[ ] HFO correcto
[ ] coordenadas en mm
[ ] D = pdist2(group, lesion)
[ ] MinDistance = mean(min(D,[],2))
[ ] MeanDistance = mean(D(:))
[ ] MaxDistance = mean(max(D,[],2))
[ ] métricas adicionales coherentes
[ ] área en cm²
[ ] escalas gráficas comunes
[ ] HTML guardado
[ ] MAT guardado
~~~

---

# 46. Prioridades de desarrollo futuras

Orden recomendado:

1. **Validar y automatizar scouts**.
2. **Conectar source map → scout → métricas sin exportación manual**.
3. **Automatizar copia de head model/covariances**.
4. **Mejorar especificidad HFO con filtro de artefactos**.
5. **Añadir score de confianza**.
6. **Bootstrap / estabilidad**.
7. **Validación con SEEG, resección y outcome**.
8. **Informe clínico automático por paciente**.

Objetivo final:

~~~text
RAW
 ↓
detección
 ↓
Spike/HFO
 ↓
epochs + averages
 ↓
source imaging
 ↓
auto-scout
 ↓
lesion comparison
 ↓
metrics
 ↓
QC
 ↓
HTML clínico
~~~

---

# 47. Contrato de trabajo para un nuevo chat

Un nuevo chat que reciba este documento debe seguir estas reglas:

1. **Continuar el proyecto, no reinventarlo.**
2. Identificar primero si la tarea afecta:
   - proceso Brainstorm;
   - pipeline MATLAB;
   - source imaging;
   - scout;
   - resultados individuales;
   - cohorte.
3. Preservar lógica no solicitada.
4. Si se modifica un parámetro científico, explicarlo.
5. Mantener compatibilidad con Brainstorm.
6. Consultar código oficial de Brainstorm si la API no es segura.
7. Distinguir defaults de interfaz de parámetros realmente usados en el pipeline.
8. No usar el seed de la lesión como máximo.
9. Mantener las definiciones actuales de distancia.
10. Entregar código completo y ejecutable cuando se pida un proceso final.

Las definiciones canónicas actuales son:

~~~matlab
D = pdist2(group_coords_mm, lesion_coords_mm);

dmin  = min(D, [], 2);
dmean = mean(D, 2);
dmax  = max(D, [], 2);

MinDistance_mm  = mean(dmin);
MeanDistance_mm = mean(dmean);   % equivalente a mean(D(:))
MaxDistance_mm  = mean(dmax);
~~~

Versiones finales del pipeline de eventos:

~~~text
HFO:
process_evt_detect_hfos_original_raw_no_bads.m

Overlap:
process_evt_detect_overlap_keep_spike.m
~~~

Análisis individual actual:

~~~text
results_tfg_7.mlx
~~~

Objetivo de evolución:

~~~text
source imaging → scout automático → métricas → informe
~~~

---

# 48. Resumen ultracorto para recuperar contexto rápidamente

~~~text
PROYECTO
Epilepsia focal farmacorresistente.
EEG ~72 canales.
MATLAB + Brainstorm.
Grupos: Spike / Spike+HFO / HFO.
Lesión MRI = referencia anatómica, no actividad.

HFO FINAL
process_evt_detect_hfos_original_raw_no_bads
- raw original
- notch + band-pass en memoria
- Hilbert
- median + k*MAD
- histéresis
- duración/ciclos/ZC
- excluye BAD
- merge intra-canal
- merge multicanal
- guarda HFO_candidate sobre raw original

OVERLAP FINAL
process_evt_detect_overlap_keep_spike
- overlap relativo al evento más corto
- pipeline actual usa 10 %
- conserva ventana completa de la spike
- output: Spike/HFO

PIPELINE 1
epilepsy_pipeline.m
- HFO
- spikes extendidas [-0.20,0]
- overlap
- simple events usando end
- epochs ±0.20 s
- averages HFO y Spike/HFO

PIPELINE 2
source_imaing_script.m
- Spike source ya existente
- HFO average → LCMV/NAI
- Spike/HFO average → LCMV/NAI

SCOUTS
Antes manuales.
Objetivo actual: máximo global + threshold + componente conectado.
Seed de lesión NO es máximo.

RESULTADOS
results_tfg_7.mlx
D = pdist2(scout, lesion)

MinDistance  = mean(min(D,2))
MeanDistance = mean(D(:))
MaxDistance  = mean(max(D,2))

Figuras:
raincloud + violin/box/scatter,
misma escala entre métricas.

COHORTE
~40 pacientes.
Unidad estadística = paciente.
Friedman + Kendall W + Wilcoxon/Holm.

REGLA PRINCIPAL
No cambiar lógica que el usuario no haya pedido cambiar.
~~~

---

# 49. Nota final

Este documento debe actualizarse cada vez que una nueva versión se declare **final**.

Cuando exista discrepancia entre:

1. memoria/conversación,
2. README histórico,
3. código actual del repositorio,

la prioridad para reproducir ejecución debe ser:

> **código actual realmente ejecutado + parámetros de llamada del pipeline**, documentando cualquier diferencia con la intención metodológica.

No debe corregirse silenciosamente una discrepancia. Primero se identifica, se explica y luego se decide si se conserva por reproducibilidad o se corrige generando una nueva versión.
