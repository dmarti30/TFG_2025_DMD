process_evt_detect_hfo_candidates_by_channel.m
Genera los HFOs para un solo canal que decidas, coge los segmentos malos (1ªversion), se hace sobre bandpass (60-250)

process_evt_detect_hfo_allch.m
Genera los HFOs generando un set de eventos para cada uno de los 72 canales, coge los segmentos malos (2ªversion) se hace sobre bandpass (60-250)

process_evt_detect_hfo_candidates.m
Genera los HFOs de todos los canales en un solo set de eventos, coge los segmentos malos y se hace sobre band-pass (60-250) (3ªversion)

process_evt_detect_hfos_without_bad
Genera los HFOs de todos los canales en un solo set de eventos sin coger los candidatos en segmentos malos, se hace sobre bandpass (60-250)

process_evt_detect_hfos_original_raw_no_bads
Genera los HFOs de todos los canales en un set de eventos, no coge candidatos en bad segments, se hace sobre el raw original, (hace notch, badnpass, detecta y se importa sobre el original) además de unir candidatos cercanos aunque se produzcan en distintos canales (4ªversion)


process_evt_detect_overlap_extended
proceso para detectar el overlap entre dos eventos extendidos bajo un porcentaje y generando en un set solo los eventos combinados (1ªversion)

process_evt_detect_overlap_keep_spike
proceso para detectar el overlap entre dos eventos extendidos bajo un procentaje y genera un set con el time-window del primer evento especificado (spike)