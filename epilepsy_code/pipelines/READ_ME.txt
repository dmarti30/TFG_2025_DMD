Workflow completo (detallado paso a paso) para obtener desde raw-->source imaging+scouts

(MANUAL)-Pasos iniciales previos a cualquier script 
1)Visualizar si raw emplea 256hz o 512hz.
2)Duplicar los eventos poblaciones de spikes para renombrarlo "spikes"

(MANUAL)-Aplicar epilepsy_pipeline
Para poder aplicar el pipeline correctamente:
1)Importar en linea 5 el directorio del raw original y en linea 7 indicar nombre del paciente
2)Modificar lineas 25-27 en función de la frecuencia del raw(en caso de raw a 256, highpass 60 y lowpass 120) (en caso de raw a 512hz, highpass 60 y lowpass 250)

(AUTOMÁTICO)-Pipeline, lo que hace es aplicar un notch, un bandpass, calcular los hfos, reimportarlos al original y eliminar los archivos no deseados. Aplica el estudio de concordancia temporal, convierte los eventos en single para importar sus epochs y hace average de grupos experimentales HFO y Spikes_HFO.

(MANUAL)-Pasos intermedios copiar head model, noise covariance y data covariance (click derecho-->copy to toher folders)

(MANUAL)-Aplicar segundo pipeline; source_imaging_script
1) importar directorio de sourceimaging de spike en el primer input files
2) importar directorio de average hfos en el segundo input files
3) importar directorio de average spikes_hfo en el tercer input files

(AUTOMATICO)-Aplicar pipeline; genera las source imaging automaticamente.


(MANUAL)-Extraer scouts para cada grupo experimental; de manera manual, extraes valor máximo del scout (para cada grupo experimental) y recreas su tamaño y ajustas hasta la actividad visible (sección scout, sources)