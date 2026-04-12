%% ========================================================================
%% CALCULADORA DE DISTANCIAS ESI ? SCOUT
%% Para análisis multi-paciente y multi-método
%% ========================================================================
%
% INSTRUCCIONES DE USO:
% 1. Exportar desde Brainstorm:
%    - Superficie cortical ? "Export to MATLAB" (como Cortex_surface)
%    - Scout de lesión ? "Export to MATLAB" (anotar nombre)
%    - Usar "Find maximum" en resultado ESI ? anotar vertex number
%
% 2. Ejecutar este script en MATLAB
%
% 3. Los resultados se guardan automáticamente en:resultados = calcular_distancias_ESI_scout();
%    - Archivo .mat con todas las variables
%    - Archivo .txt con reporte legible
%    - Figura con histograma (opcional)
%
% Autor: [Tu nombre]
% Fecha: octubre 2025
%% ========================================================================

function resultados = calcular_distancias_ESI_scout(varargin)
    % Configuración por defecto
    cfg = struct();
    cfg.patient_id = 'PACIENTE_001';
    cfg.method = 'sLORETA';
    cfg.esi_vertex = [];
    cfg.scout_indices = [];
    cfg.cortex_var = 'Cortex_surface';
    cfg.output_dir = fullfile(getenv('HOME'), 'Desktop', 'ESI_GAMMA', 'resultados');
    cfg.save_results = true;
    cfg.plot_histogram = true;
    cfg.interactive = true;
    
    % Parsear argumentos opcionales
    if nargin > 0
        for i = 1:2:length(varargin)
            cfg.(varargin{i}) = varargin{i+1};
        end
    end
    
    %% ====================================================================
    %% MODO INTERACTIVO
    %% ====================================================================
    if cfg.interactive
        fprintf('\n');
        fprintf('========================================\n');
        fprintf('  CALCULADORA DISTANCIAS ESI ? SCOUT\n');
        fprintf('========================================\n\n');
        
        % Solicitar información del paciente
        if isempty(cfg.patient_id) || strcmp(cfg.patient_id, 'PACIENTE_001')
            cfg.patient_id = input('ID del paciente: ', 's');
            if isempty(cfg.patient_id)
                cfg.patient_id = 'PACIENTE_001';
            end
        end
        
        % Solicitar método ESI
        fprintf('\nMétodos ESI disponibles:\n');
        fprintf('  1. sLORETA\n');
        fprintf('  2. MNE\n');
        fprintf('  3. dSPM\n');
        fprintf('  4. LCMV\n');
        fprintf('  5. Otro\n');
        metodo_num = input('Seleccione método (1-5) [1]: ');
        if isempty(metodo_num)
            metodo_num = 1;
        end
        metodos = {'sLORETA', 'MNE', 'dSPM', 'LCMV', 'Otro'};
        cfg.method = metodos{metodo_num};
        if metodo_num == 5
            cfg.method = input('Especifique método: ', 's');
        end
        
        % Solicitar vertex del máximo
        fprintf('\n? Use "Find maximum" en Brainstorm para obtener el vertex number\n');
        if isempty(cfg.esi_vertex)
            cfg.esi_vertex = input('Vertex number del máximo ESI: ');
        end
        
        % Verificar variables en workspace
        fprintf('\n? Verificando variables en workspace...\n');
        vars = evalin('base', 'who');
        
        % Buscar superficie cortical
        cortex_found = false;
        for i = 1:length(vars)
            if contains(lower(vars{i}), 'cortex') || contains(lower(vars{i}), 'surface')
                fprintf('   ? Encontrada: %s\n', vars{i});
                use_this = input(sprintf('   ¿Usar %s como superficie cortical? (s/n) [s]: ', vars{i}), 's');
                if isempty(use_this) || lower(use_this) == 's'
                    cfg.cortex_var = vars{i};
                    cortex_found = true;
                    break;
                end
            end
        end
        
        if ~cortex_found
            error('? No se encontró superficie cortical. Exportar desde Brainstorm primero.');
        end
        
        % Buscar scout
        fprintf('\n? Buscando scout exportado...\n');
        scout_found = false;
        for i = 1:length(vars)
            if contains(lower(vars{i}), 'scout') || contains(lower(vars{i}), 'mask') || contains(lower(vars{i}), 'lesion')
                fprintf('   ? Encontrado: %s\n', vars{i});
                use_this = input(sprintf('   ¿Usar %s como scout? (s/n) [s]: ', vars{i}), 's');
                if isempty(use_this) || lower(use_this) == 's'
                    scout_var = evalin('base', vars{i});
                    if isfield(scout_var, 'Scouts') && isfield(scout_var.Scouts, 'Vertices')
                        cfg.scout_indices = scout_var.Scouts.Vertices;
                    elseif isfield(scout_var, 'Vertices')
                        cfg.scout_indices = scout_var.Vertices;
                    end
                    scout_found = true;
                    break;
                end
            end
        end
        
        % Si no se encontró scout, solicitar índices manualmente
        if ~scout_found || isempty(cfg.scout_indices)
            fprintf('\n??  No se encontró scout automáticamente.\n');
            fprintf('   Opción 1: Exportar scout desde Brainstorm\n');
            fprintf('   Opción 2: Pegar índices manualmente\n');
            opcion = input('Seleccione opción (1/2): ');
            if opcion == 2
                fprintf('\n? Pegue los índices del scout (formato: [1 2 3 ...]): ');
                cfg.scout_indices = input('');
            else
                error('? Scout no encontrado. Exportar desde Brainstorm.');
            end
        end
        
        fprintf('\n');
    end
    
    %% ====================================================================
    %% CARGAR DATOS
    %% ====================================================================
    fprintf('? Cargando datos...\n');
    
    % Cargar superficie cortical
    cortex_data = evalin('base', cfg.cortex_var);
    
    if ~isfield(cortex_data, 'Vertices')
        error('? La variable %s no contiene el campo Vertices', cfg.cortex_var);
    end
    
    fprintf('   ? Superficie cortical: %d vértices\n', size(cortex_data.Vertices, 1));
    fprintf('   ? Scout: %d vértices\n', length(cfg.scout_indices));
    fprintf('   ? Vértice máximo ESI: %d\n', cfg.esi_vertex);
    
    %% ====================================================================
    %% CALCULAR DISTANCIAS
    %% ====================================================================
    fprintf('\n? Calculando distancias...\n');
    
    % Obtener coordenadas
    scout_coords = cortex_data.Vertices(cfg.scout_indices, :);
    esi_coords = cortex_data.Vertices(cfg.esi_vertex, :);
    
    % Calcular distancias euclidianas (metros ? mm)
    distances = sqrt(sum((scout_coords - esi_coords).^2, 2)) * 1000;
    
    % Calcular métricas
    resultados = struct();
    resultados.patient_id = cfg.patient_id;
    resultados.method = cfg.method;
    resultados.date = datetime('now');
    resultados.esi_vertex = cfg.esi_vertex;
    resultados.esi_coords_m = esi_coords;
    resultados.scout_n_vertices = length(cfg.scout_indices);
    resultados.scout_indices = cfg.scout_indices;
    resultados.distances_mm = distances;
    
    % Métricas estadísticas
    resultados.min_distance_mm = min(distances);
    resultados.max_distance_mm = max(distances);
    resultados.mean_distance_mm = mean(distances);
    resultados.median_distance_mm = median(distances);
    resultados.std_distance_mm = std(distances);
    % Calcular cuartiles sin Statistics Toolbox
    sorted_dist = sort(distances);
    n = length(sorted_dist);
    q25_idx = round(n * 0.25);
    q75_idx = round(n * 0.75);
    resultados.q25_distance_mm = sorted_dist(max(1, q25_idx));
    resultados.q75_distance_mm = sorted_dist(min(n, q75_idx));
    
    % Verificar si el máximo está dentro del scout
    resultados.inside_scout = ismember(cfg.esi_vertex, cfg.scout_indices);
    
    %% ====================================================================
    %% MOSTRAR RESULTADOS
    %% ====================================================================
    fprintf('\n');
    fprintf('========================================\n');
    fprintf('         RESULTADOS DE ANÁLISIS\n');
    fprintf('========================================\n');
    fprintf('Paciente:       %s\n', resultados.patient_id);
    fprintf('Método ESI:     %s\n', resultados.method);
    fprintf('Fecha:          %s\n', datestr(resultados.date));
    fprintf('----------------------------------------\n');
    fprintf('Vértice máximo: %d\n', resultados.esi_vertex);
    fprintf('Coordenadas:    [%.4f, %.4f, %.4f] m\n', esi_coords);
    fprintf('Scout vértices: %d\n', resultados.scout_n_vertices);
    fprintf('========================================\n');
    fprintf('          DISTANCIAS (mm)\n');
    fprintf('========================================\n');
    fprintf('? MÍNIMA:      %.2f mm\n', resultados.min_distance_mm);
    fprintf('? MÁXIMA:      %.2f mm\n', resultados.max_distance_mm);
    fprintf('? MEDIA:       %.2f mm\n', resultados.mean_distance_mm);
    fprintf('? MEDIANA:     %.2f mm\n', resultados.median_distance_mm);
    fprintf('? Desv. Est.:  %.2f mm\n', resultados.std_distance_mm);
    fprintf('   Q25:         %.2f mm\n', resultados.q25_distance_mm);
    fprintf('   Q75:         %.2f mm\n', resultados.q75_distance_mm);
    fprintf('========================================\n\n');
    
    % Interpretación clínica
    if resultados.inside_scout
        fprintf('? RESULTADO: El máximo de %s está DENTRO del scout\n', cfg.method);
        fprintf('   ? Concordancia perfecta entre ESI y lesión estructural\n');
    elseif resultados.min_distance_mm < 5
        fprintf('? RESULTADO: El máximo está MUY CERCA del scout (%.2f mm)\n', resultados.min_distance_mm);
        fprintf('   ??  Posible concordancia con margen de error de localización\n');
    elseif resultados.min_distance_mm < 10
        fprintf('? RESULTADO: El máximo está CERCA del scout (%.2f mm)\n', resultados.min_distance_mm);
        fprintf('   ??  Revisar registración y métodos de inversión\n');
    else
        fprintf('? RESULTADO: El máximo está LEJOS del scout (%.2f mm)\n', resultados.min_distance_mm);
        fprintf('   ??  Posible discordancia - revisar análisis\n');
    end
    fprintf('\n');
    
    %% ====================================================================
    %% VISUALIZACIÓN
    %% ====================================================================
    if cfg.plot_histogram
        figure('Name', sprintf('%s - %s: Distribución de Distancias', ...
            resultados.patient_id, resultados.method), ...
            'Position', [100, 100, 800, 600]);
        
        histogram(distances, 20, 'FaceColor', [0.2 0.6 0.8], 'EdgeColor', 'k');
        hold on;
        
        % Líneas de referencia
        xline(resultados.min_distance_mm, 'r--', 'LineWidth', 2, ...
            'Label', sprintf('Mín: %.1f mm', resultados.min_distance_mm));
        xline(resultados.mean_distance_mm, 'g--', 'LineWidth', 2, ...
            'Label', sprintf('Media: %.1f mm', resultados.mean_distance_mm));
        xline(resultados.max_distance_mm, 'b--', 'LineWidth', 2, ...
            'Label', sprintf('Máx: %.1f mm', resultados.max_distance_mm));
        
        xlabel('Distancia (mm)', 'FontSize', 12);
        ylabel('Frecuencia', 'FontSize', 12);
        title(sprintf('Distribución de distancias: %s - %s\nVértice máximo %d ? Scout (%d vértices)', ...
            resultados.patient_id, resultados.method, resultados.esi_vertex, resultados.scout_n_vertices), ...
            'FontSize', 14, 'FontWeight', 'bold');
        grid on;
        hold off;
    end
    
    %% ====================================================================
    %% GUARDAR RESULTADOS
    %% ====================================================================
    if cfg.save_results
        % Crear directorio de salida si no existe
        if ~exist(cfg.output_dir, 'dir')
            mkdir(cfg.output_dir);
        end
        
        % Nombre base del archivo
        timestamp = datestr(now, 'yyyymmdd_HHMMSS');
        base_filename = sprintf('%s_%s_%s', resultados.patient_id, resultados.method, timestamp);
        
        % Guardar archivo .mat
        mat_filename = fullfile(cfg.output_dir, [base_filename '.mat']);
        save(mat_filename, 'resultados', 'cfg');
        fprintf('? Resultados guardados en: %s\n', mat_filename);
        
        % Guardar reporte de texto
        txt_filename = fullfile(cfg.output_dir, [base_filename '.txt']);
        fid = fopen(txt_filename, 'w');
        fprintf(fid, '========================================\n');
        fprintf(fid, 'ANÁLISIS DE DISTANCIAS ESI ? SCOUT\n');
        fprintf(fid, '========================================\n');
        fprintf(fid, 'Paciente:        %s\n', resultados.patient_id);
        fprintf(fid, 'Método ESI:      %s\n', resultados.method);
        fprintf(fid, 'Fecha:           %s\n', datestr(resultados.date));
        fprintf(fid, '----------------------------------------\n');
        fprintf(fid, 'Vértice máximo:  %d\n', resultados.esi_vertex);
        fprintf(fid, 'Coordenadas SCS: [%.4f, %.4f, %.4f] m\n', esi_coords);
        fprintf(fid, 'Scout vértices:  %d\n', resultados.scout_n_vertices);
        fprintf(fid, '========================================\n');
        fprintf(fid, 'DISTANCIAS (mm)\n');
        fprintf(fid, '========================================\n');
        fprintf(fid, 'Mínima:          %.2f mm\n', resultados.min_distance_mm);
        fprintf(fid, 'Máxima:          %.2f mm\n', resultados.max_distance_mm);
        fprintf(fid, 'Media:           %.2f mm\n', resultados.mean_distance_mm);
        fprintf(fid, 'Mediana:         %.2f mm\n', resultados.median_distance_mm);
        fprintf(fid, 'Desv. Est.:      %.2f mm\n', resultados.std_distance_mm);
        fprintf(fid, 'Q25:             %.2f mm\n', resultados.q25_distance_mm);
        fprintf(fid, 'Q75:             %.2f mm\n', resultados.q75_distance_mm);
        fprintf(fid, '========================================\n');
        fprintf(fid, 'Dentro del scout: %s\n', mat2str(resultados.inside_scout));
        fprintf(fid, '========================================\n');
        fclose(fid);
        fprintf('? Reporte guardado en: %s\n', txt_filename);
        
        % Guardar figura si se creó
        if cfg.plot_histogram
            fig_filename = fullfile(cfg.output_dir, [base_filename '.png']);
            saveas(gcf, fig_filename);
            fprintf('??  Figura guardada en: %s\n', fig_filename);
        end
        
        fprintf('\n? Todos los archivos guardados en: %s\n', cfg.output_dir);
    end
    
    fprintf('\n');
end

%% ========================================================================
%% EJEMPLOS DE USO
%% ========================================================================
%{

%% EJEMPLO 1: Modo interactivo (recomendado)
% 1. Exportar desde Brainstorm: superficie cortical y scout
% 2. Ejecutar:
resultados = calcular_distancias_ESI_scout();

%% EJEMPLO 2: Modo no interactivo con parámetros específicos
resultados = calcular_distancias_ESI_scout(...
    'patient_id', 'PAC_001', ...
    'method', 'sLORETA', ...
    'esi_vertex', 12080, ...
    'scout_indices', [11418 11507 11546 11595 ...], ...
    'cortex_var', 'Cortex_surface', ...
    'interactive', false);

%% EJEMPLO 3: Procesar múltiples pacientes
pacientes = {'PAC_001', 'PAC_002', 'PAC_003'};
vertices_max = [12080, 15234, 8956];

for i = 1:length(pacientes)
    fprintf('\n===== Procesando %s =====\n', pacientes{i});
    resultados{i} = calcular_distancias_ESI_scout(...
        'patient_id', pacientes{i}, ...
        'esi_vertex', vertices_max(i), ...
        'interactive', false);
end

%% EJEMPLO 4: Comparar múltiples métodos ESI para el mismo paciente
metodos = {'sLORETA', 'MNE', 'dSPM', 'LCMV'};
vertices = [12080, 12100, 12050, 12090];

for i = 1:length(metodos)
    resultados_metodos{i} = calcular_distancias_ESI_scout(...
        'patient_id', 'PAC_001', ...
        'method', metodos{i}, ...
        'esi_vertex', vertices(i), ...
        'interactive', false);
end

% Comparar resultados
comparacion = struct();
for i = 1:length(metodos)
    comparacion.(metodos{i}) = resultados_metodos{i}.min_distance_mm;
end
disp(comparacion);

%}