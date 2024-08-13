% Número de indivíduos na população
popSize = 100;

% Número de gerações
numGenerations = 100;

% Taxa de crossover
crossoverRate = 1;

% Taxa de mutação
mutationRate = 0.15;

% Número de Parâmetros
paramNum = 4;

% Limites dos parâmetros [Rs, Rr, Lm, Lls, Llr]
paramBounds = [0.1, 5; % Rs (0.6759)
    0.1, 5; % Rr (0.2615)
    0.001, 0.1; % Lm (0.0387)
    0.0001, 0.1]; % Lls (0.0028)

func = 'f_motor_high_data';
% func = 'f_motor_low_data';

% Tipos de seleção ->
% * Torneio (será feito o promp do número de integrantes)
% * Roleta
% * Rank
% * Estoc
% * Boltzmann (Será feito o prompt da temperatura inicial e final -> padrão 100 e 1

% Tipos de Crossover ->
% 1P
% 2P
% Arit
% Uniform
% Blend

% Tipos de Mutação ->
% Default
% Gauss
% Uniform
% Creep
% Inv


% Definindo os parâmetros
% selectionMethods = {'Torneio', 'Roleta', 'Rank', 'Estoc', 'Boltzmann'};
selectionMethods = {'Torneio', 'Roleta', 'Estoc', 'Boltzmann'};
% crossoverMethods = {'1P', '2P', 'Arit', 'Uniform', 'Blend'};
crossoverMethods = {'1P', '2P', 'Uniform', 'Blend'};
% mutationMethods = {'Default', 'Gauss', 'Uniform', 'Creep', 'Inv'};
mutationMethods = {'Gauss', 'Uniform', 'Creep'};

% Argumentos adicionais para os métodos de seleção
% selectionArgsArray = {3, [], [], 2, [100, 1]}; % Ajuste conforme necessário
selectionArgsArray = {3, [], 2, [100, 1]}; % Ajuste conforme necessário

numExecucoes = 5;

% Arrays para armazenar os dados acumulados
allBestFitness = zeros(numExecucoes, length(selectionMethods) * length(crossoverMethods) * length(mutationMethods));
allExecutionTimes = zeros(size(allBestFitness));
varNames = cell(1, size(allBestFitness, 2));

% Loop para número de execuções
for loopAtual = 1:numExecucoes
    % Resetar índice para nomes das variáveis
    varIndex = 1;

    % Loop para automatizar o processo
    for selIdx = 1:length(selectionMethods)
        for crossIdx = 1:length(crossoverMethods)
            for mutIdx = 1:length(mutationMethods)
                % Criar nome da variável
                varName = sprintf('exec%d_ga_%s_%s_%s', loopAtual, lower(selectionMethods{selIdx}), crossoverMethods{crossIdx}, lower(mutationMethods{mutIdx}));

                % Exibir informações
                disp(['**************** Run: ' num2str(loopAtual) ' -> '  selectionMethods{selIdx} ' + ' crossoverMethods{crossIdx} ' + ' mutationMethods{mutIdx} ' ****************']);

                % Criar e inicializar o objeto GA
                gaObj = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
                gaObj = gaObj.Initialize;

                % Preparar argumentos de seleção como uma célula
                if strcmp(selectionMethods{selIdx}, 'Estoc') || strcmp(selectionMethods{selIdx}, 'Boltzmann') || strcmp(selectionMethods{selIdx}, 'Torneio')
                    selectionArgs = {selectionMethods{selIdx}, selectionArgsArray{selIdx}};
                else
                    selectionArgs = selectionMethods(selIdx);
                end

                % Preparar argumentos de crossover como uma célula
                if strcmp(crossoverMethods{crossIdx}, 'Blend')
                    crossoverArgs = {crossoverMethods{crossIdx}, 0.5};
                else
                    crossoverArgs = crossoverMethods(crossIdx);
                end

                % Resolver com os métodos selecionados
                gaObj = gaObj.Solve(selectionArgs, crossoverArgs, mutationMethods{mutIdx});

                % Salvar o resultado na variável apropriada
                eval([varName ' = gaObj;']);

                % Salvar o nome da variável e coletar dados
                varNames{varIndex} = varName;
                allBestFitness(loopAtual, varIndex) = gaObj.bestFitness;
                allExecutionTimes(loopAtual, varIndex) = gaObj.executionTime;

                % Incrementar índice para nomes das variáveis
                varIndex = varIndex + 1;

                % Opcional: plotar resultados
                % gaObj.plotMotor;
                % gaObj.plotGA;
            end
        end
    end
end

% Gráfico de bestFitness
figure;
bar(allBestFitness);
set(gca, 'xticklabel', varNames);
xtickangle(45);
ylabel('Best Fitness');
title('Comparação de Best Fitness para cada Configuração do GA');
grid on;

% Gráfico de executionTime
figure;
bar(allExecutionTimes);
set(gca, 'xticklabel', varNames);
xtickangle(45);
ylabel('Execution Time (seconds)');
title('Comparação de Execution Time para cada Configuração do GA');
grid on;

% Calcular médias
meanBestFitness = median(allBestFitness, 1);
meanExecutionTimes = median(allExecutionTimes, 1);

% Calcular o valor mínimo e máximo de meanBestFitness
valorMinimo = min(meanBestFitness);
valorMaximo = max(meanBestFitness);
mediaGeral = median(meanBestFitness);

% % Gráfico de Média de Best Fitness
% figure;
% bar(meanBestFitness);
% set(gca, 'xticklabel', varNames);
% xtickangle(90);
% ylabel('Média do Best Fitness');
% title('Média do Best Fitness para cada Configuração do GA');
% grid on;
% 
% % Gráfico de Média de Execution Time
% figure;
% bar(meanExecutionTimes);
% set(gca, 'xticklabel', varNames);
% xtickangle(90);
% ylabel('Média do Execution Time (segundos)');
% title('Média do Execution Time para cada Configuração do GA');
% grid on;

%% Plotagem dos resultados por método de seleção
% Cores para cada tipo de seleção
colors = {'r', 'g', 'b', 'c', 'm'}; % Adicione ou altere as cores conforme necessário

% Número de combinações por tipo de seleção
numCombinacoesPorSelecao = length(crossoverMethods) * length(mutationMethods);

figure;
hold on; % Permite sobrepor várias séries de barras no mesmo gráfico

% Loop para cada método de seleção
for i = 1:length(selectionMethods)
    % Índices das barras correspondentes a este método de seleção
    idx = (i-1)*numCombinacoesPorSelecao + 1 : i*numCombinacoesPorSelecao;

    % Plotar barras para este método de seleção
    bar(idx, meanBestFitness(idx), colors{i});
end

hold off;
ylabel('Média do Best Fitness');
title(sprintf('SELEÇÃO -> gerações: %d população: %d ', numGenerations, popSize)); 
grid on;
hold on;
 plot(xlim, [valorMinimo valorMinimo], ':r', 'LineWidth', 2);  % Linha do mínimo
plot(xlim, [valorMaximo valorMaximo], ':g', 'LineWidth', 2);  % Linha do máximo
 plot(xlim, [mediaGeral mediaGeral], '--k', 'LineWidth', 2);  % Linha da média
legend([selectionMethods, 'Valor Mínimo', 'Valor Máximo', 'Média']);
hold off;

% Plotagem dos resultados por método de crossover
% Cores para cada tipo de crossover
colors = {'r', 'g', 'b', 'c', 'm', 'y'}; % Adicione ou altere as cores conforme necessário

% Número total de grupos (combinações)
numGrupos = length(meanBestFitness) / length(crossoverMethods);

% Largura das barras
larguraBarra = 1 / (length(crossoverMethods) + 1);

% Criar um gráfico de barras
figure;
hold on;

% Plotar as barras
for i = 1:length(crossoverMethods)
    % Posições das barras para este método de crossover
    posicoes = (1:numGrupos) - 0.5 + (i-1) * larguraBarra;
    
    % Dados das barras para este método de crossover
    dados = meanBestFitness((i-1)*numGrupos + 1:i*numGrupos);

    % Plotar as barras
    bar(posicoes, dados, larguraBarra, 'FaceColor', colors{i});
end

hold off;
ylabel('Média do Best Fitness');
title(sprintf('CROSSOVER -> gerações: %d população: %d ', numGenerations, popSize)); 
grid on;
hold on;
 plot(xlim, [valorMinimo valorMinimo], ':r', 'LineWidth', 2);  % Linha do mínimo
plot(xlim, [valorMaximo valorMaximo], ':g', 'LineWidth', 2);  % Linha do máximo
 plot(xlim, [mediaGeral mediaGeral], '--k', 'LineWidth', 2);  % Linha da média
legend([crossoverMethods, 'Valor Mínimo', 'Valor Máximo', 'Média']);
hold off;

% Plotagem dos métodos por mutação
% Cores para cada tipo de mutação
colors = {'r', 'g', 'b', 'c', 'm', 'y'}; % Adicione ou altere as cores conforme necessário

% Número de grupos (combinações de seleção e crossover)
numGrupos = length(selectionMethods) * length(crossoverMethods);

% Preparar dados para barras agrupadas
dadosBarras = zeros(length(mutationMethods), numGrupos);
for i = 1:length(mutationMethods)
    idx = i:length(mutationMethods):length(meanBestFitness);
    dadosBarras(i, :) = meanBestFitness(idx);
end

% Criar um gráfico de barras agrupadas
figure;
b = bar(dadosBarras', 'grouped');

% Ajustar a largura das barras
% Um valor maior que 1 faz as barras ficarem mais largas
for i = 1:length(b)
    b(i).BarWidth = 0.9;  % Ajustar a largura das barras
end

% Aplicar cores
for i = 1:length(b)
    b(i).FaceColor = colors{i};
end

% Configurar eixos e legendas
ylabel('Média do Best Fitness');
title(sprintf('MUTAÇÃO -> gerações: %d população: %d', numGenerations, popSize)); 
grid on;
hold on;
plot(xlim, [valorMinimo valorMinimo], ':r', 'LineWidth', 2);  % Linha do mínimo
 plot(xlim, [valorMaximo valorMaximo], ':g', 'LineWidth', 2);  % Linha do máximo
plot(xlim, [mediaGeral mediaGeral], '--k', 'LineWidth', 2);  % Linha da média
legend([mutationMethods, 'Valor Mínimo', 'Valor Máximo', 'Média']);
hold off;



%%
% %
% disp('**************** Torneio + 1P + Default ****************');
% ga_torneio_1p_default = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_torneio_1p_default = ga_torneio_1p_default.Initialize;
% ga_torneio_1p_default = ga_torneio_1p_default.Solve({'Torneio', 3}, '1P', 'Default');
% % ga_torneio_1p_default = ga_torneio_1p_default.plotMotor;
% % ga_torneio_1p_default.plotGA;

% %
% disp('**************** Roleta + 1P + Default ****************');
% ga_roleta_1p_default = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_roleta_1p_default = ga_roleta_1p_default.Initialize;
% ga_roleta_1p_default = ga_roleta_1p_default.Solve('Roleta', '1P', 'Default');
% % ga_roleta_1p_default = ga_roleta_1p_default.plotMotor;
% % ga_roleta_1p_default.plotGA;
%
% %
% disp('**************** Rank + 1P + Default ****************');
% ga_rank_1p_default = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_rank_1p_default = ga_rank_1p_default.Initialize;
% ga_rank_1p_default = ga_rank_1p_default.Solve('Rank', '1P', 'Default');
% % ga_rank_1p_default = ga_rank_1p_default.plotMotor;
% % ga_rank_1p_default.plotGA;
%
% %
% disp('**************** Estoc + 1P + Default ****************');
% ga_estoc_1p_default = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_estoc_1p_default = ga_estoc_1p_default.Initialize;
% ga_estoc_1p_default = ga_estoc_1p_default.Solve({'Estoc',1} , '1P', 'Default');
% % ga_estoc_1p_default = ga_estoc_1p_default.plotMotor;
% % ga_estoc_1p_default.plotGA;
%
% %
% disp('**************** Boltz + 1P + Default ****************');
% ga_boltz_1p_default = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_boltz_1p_default = ga_boltz_1p_default.Initialize;
% ga_boltz_1p_default = ga_boltz_1p_default.Solve({'Boltzmann',100 , 1}, '1P', 'Default');
% % ga_boltz_1p_default = ga_boltz_1p_default.plotMotor;
% % ga_boltz_1p_default.plotGA;
%
% %%
% %
% disp('**************** Torneio + 2P + Default ****************');
% ga_torneio_2p_default = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_torneio_2p_default = ga_torneio_2p_default.Initialize;
% ga_torneio_2p_default = ga_torneio_2p_default.Solve({'Torneio', 3}, '2P', 'Default');
% % ga_torneio_2p_default = ga_torneio_2p_default.plotMotor;
% % ga_torneio_2p_default.plotGA;
%
% %
% disp('**************** Roleta + 2P + Default ****************');
% ga_roleta_2p_default = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_roleta_2p_default = ga_roleta_2p_default.Initialize;
% ga_roleta_2p_default = ga_roleta_2p_default.Solve('Roleta', '2P', 'Default');
% % ga_roleta_2p_default = ga_roleta_2p_default.plotMotor;
% % ga_roleta_2p_default.plotGA;
%
% %
% disp('**************** Rank + 2P + Default ****************');
% ga_rank_2p_default = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_rank_2p_default = ga_rank_2p_default.Initialize;
% ga_rank_2p_default = ga_rank_2p_default.Solve('Rank', '2P', 'Default');
% % ga_rank_2p_default = ga_rank_2p_default.plotMotor;
% % ga_rank_2p_default.plotGA;
%
% %
% disp('**************** Estoc + 2P + Default ****************');
% ga_estoc_2p_default = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_estoc_2p_default = ga_estoc_2p_default.Initialize;
% ga_estoc_2p_default = ga_estoc_2p_default.Solve({'Estoc',2}, '2P', 'Default');
% % ga_estoc_2p_default = ga_estoc_2p_default.plotMotor;
% % ga_estoc_2p_default.plotGA;
%
% %
% disp('**************** Boltz + 2P + Default ****************');
% ga_boltz_2p_default = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_boltz_2p_default = ga_boltz_2p_default.Initialize;
% ga_boltz_2p_default = ga_boltz_2p_default.Solve({'Boltzmann',100 , 1}, '2P', 'Default');
% % ga_boltz_2p_default = ga_boltz_2p_default.plotMotor;
% % ga_boltz_2p_default.plotGA;
%
% %%
% %
% disp('**************** Torneio + Arit + Default ****************');
% ga_torneio_arit_default = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_torneio_arit_default = ga_torneio_arit_default.Initialize;
% ga_torneio_arit_default = ga_torneio_arit_default.Solve({'Torneio', 3}, 'Arit', 'Default');
% % ga_torneio_arit_default = ga_torneio_arit_default.plotMotor;
% % ga_torneio_arit_default.plotGA;
%
% %
% disp('**************** Roleta + Arit + Default ****************');
% ga_roleta_arit_default = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_roleta_arit_default = ga_roleta_arit_default.Initialize;
% ga_roleta_arit_default = ga_roleta_arit_default.Solve('Roleta', 'Arit', 'Default');
% % ga_roleta_arit_default = ga_roleta_arit_default.plotMotor;
% % ga_roleta_arit_default.plotGA;
%
% %
% disp('**************** Rank + Arit+ Default ****************');
% ga_rank_arit_default = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_rank_arit_default = ga_rank_arit_default.Initialize;
% ga_rank_arit_default = ga_rank_arit_default.Solve('Rank', 'Arit', 'Default');
% % ga_rank_arit_default = ga_rank_arit_default.plotMotor;
% % ga_rank_arit_default.plotGA;
%
% %
% disp('**************** Estoc + Arit + Default ****************');
% ga_estoc_arit_default = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_estoc_arit_default = ga_estoc_arit_default.Initialize;
% ga_estoc_arit_default = ga_estoc_arit_default.Solve({'Estoc',1}, 'Arit', 'Default');
% % ga_estoc_arit_default = ga_estoc_arit_default.plotMotor;
% % ga_estoc_arit_default.plotGA;
%
% %
% disp('**************** Boltz + Arit + Default ****************');
% ga_boltz_arit_default = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_boltz_arit_default = ga_boltz_arit_default.Initialize;
% ga_boltz_arit_default = ga_boltz_arit_default.Solve({'Boltzmann',100 , 1}, 'Arit', 'Default');
% % ga_boltz_arit_default = ga_boltz_arit_default.plotMotor;
% % ga_boltz_arit_default.plotGA;
%
% %%
% %
% disp('**************** Torneio + Unif + Default ****************');
% ga_torneio_unif_default = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_torneio_unif_default = ga_torneio_unif_default.Initialize;
% ga_torneio_unif_default = ga_torneio_unif_default.Solve({'Torneio', 3}, 'Uniform', 'Default');
% % ga_torneio_unif_default = ga_torneio_unif_default.plotMotor;
% % ga_torneio_unif_default.plotGA;
%
% %
% disp('**************** Roleta + Unif + Default ****************');
% ga_roleta_unif_default = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_roleta_unif_default = ga_roleta_unif_default.Initialize;
% ga_roleta_unif_default = ga_roleta_unif_default.Solve('Roleta', 'Uniform', 'Default');
% % ga_roleta_unif_default = ga_roleta_unif_default.plotMotor;
% % ga_roleta_unif_default.plotGA;
%
% %
% disp('**************** Rank + Unif+ Default ****************');
% ga_rank_unif_default = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_rank_unif_default = ga_rank_unif_default.Initialize;
% ga_rank_unif_default = ga_rank_unif_default.Solve('Rank', 'Uniform', 'Default');
% % ga_rank_unif_default = ga_rank_unif_default.plotMotor;
% % ga_rank_unif_default.plotGA;
%
% %
% disp('**************** Estoc + Unif + Default ****************');
% ga_estoc_unif_default = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_estoc_unif_default = ga_estoc_unif_default.Initialize;
% ga_estoc_unif_default = ga_estoc_unif_default.Solve({'Estoc',1}, 'Uniform', 'Default');
% % ga_estoc_unif_default = ga_estoc_unif_default.plotMotor;
% % ga_estoc_unif_default.plotGA;
%
% %
% disp('**************** Boltz + Unif + Default ****************');
% ga_boltz_unif_default = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_boltz_unif_default = ga_boltz_unif_default.Initialize;
% ga_boltz_unif_default = ga_boltz_unif_default.Solve({'Boltzmann',100 , 1}, 'Uniform', 'Default');
% % ga_boltz_unif_default = ga_boltz_unif_default.plotMotor;
% % ga_boltz_unif_default.plotGA;
%
% %%
% disp('**************** Torneio + Blend + Default ****************');
% ga_torneio_blend_default = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_torneio_blend_default = ga_torneio_blend_default.Initialize;
% ga_torneio_blend_default = ga_torneio_blend_default.Solve({'Torneio', 3}, {'Blend', 0.5}, 'Default');
% % ga_torneio_blend_default = ga_torneio_blend_default.plotMotor;
% % ga_torneio_blend_default.plotGA;
%
% disp('**************** Roleta + Blend + Default ****************');
% ga_roleta_blend_default = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_roleta_blend_default = ga_roleta_blend_default.Initialize;
% ga_roleta_blend_default = ga_roleta_blend_default.Solve('Roleta', {'Blend', 0.5}, 'Default');
% % ga_roleta_blend_default = ga_roleta_blend_default.plotMotor;
% % ga_roleta_blend_default.plotGA;
%
% %
% disp('**************** Rank + Blend + Default ****************');
% ga_rank_blend_default = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_rank_blend_default = ga_rank_blend_default.Initialize;
% ga_rank_blend_default = ga_rank_blend_default.Solve('Rank', {'Blend', 0.5}, 'Default');
% % ga_rank_blend_default = ga_rank_blend_default.plotMotor;
% % ga_rank_blend_default.plotGA;
%
% %
% disp('**************** Estoc + Blend + Default ****************');
% ga_estoc_blend_default = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_estoc_blend_default = ga_estoc_blend_default.Initialize;
% ga_estoc_blend_default = ga_estoc_blend_default.Solve({'Estoc',1}, {'Blend', 0.5}, 'Default');
% % ga_estoc_blend_default = ga_estoc_blend_default.plotMotor;
% % ga_estoc_blend_default.plotGA;
%
% %
% disp('**************** Boltz + Blend + Default ****************');
% ga_boltz_blend_default = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_boltz_blend_default = ga_boltz_blend_default.Initialize;
% ga_boltz_blend_default = ga_boltz_blend_default.Solve({'Boltzmann',100 , 1}, {'Blend', 0.5}, 'Default');
% % ga_boltz_blend_default = ga_boltz_blend_default.plotMotor;
% % ga_boltz_blend_default.plotGA;
%
%
%
%
%
%
% %%
% %
% disp('**************** Torneio + 1P + Gauss ****************');
% ga_torneio_1p_gauss = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_torneio_1p_gauss = ga_torneio_1p_gauss.Initialize;
% ga_torneio_1p_gauss = ga_torneio_1p_gauss.Solve({'Torneio', 3}, '1P', 'Gauss');
% % ga_torneio_1p_gauss = ga_torneio_1p_gauss.plotMotor;
% % ga_torneio_1p_gauss.plotGA;
%
% %
% disp('**************** Roleta + 1P + Gauss ****************');
% ga_roleta_1p_gauss = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_roleta_1p_gauss = ga_roleta_1p_gauss.Initialize;
% ga_roleta_1p_gauss = ga_roleta_1p_gauss.Solve('Roleta', '1P', 'Gauss');
% % ga_roleta_1p_gauss = ga_roleta_1p_gauss.plotMotor;
% % ga_roleta_1p_gauss.plotGA;
%
% %
% disp('**************** Rank + 1P + Gauss ****************');
% ga_rank_1p_gauss = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_rank_1p_gauss = ga_rank_1p_gauss.Initialize;
% ga_rank_1p_gauss = ga_rank_1p_gauss.Solve('Rank', '1P', 'Gauss');
% % ga_rank_1p_gauss = ga_rank_1p_gauss.plotMotor;
% % ga_rank_1p_gauss.plotGA;
%
% %
disp('**************** Estoc + 1P + Gauss ****************');
ga_estoc_1p_gauss = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
ga_estoc_1p_gauss = ga_estoc_1p_gauss.Initialize;
ga_estoc_1p_gauss = ga_estoc_1p_gauss.Solve({'Estoc',2} , '1P', 'Gauss');
% ga_estoc_1p_gauss = ga_estoc_1p_gauss.plotMotor;
% ga_estoc_1p_gauss.plotGA;
%
% %
% disp('**************** Boltz + 1P + Gauss ****************');
% ga_boltz_1p_gauss = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_boltz_1p_gauss = ga_boltz_1p_gauss.Initialize;
% ga_boltz_1p_gauss = ga_boltz_1p_gauss.Solve({'Boltzmann',100 , 1}, '1P', 'Gauss');
% % ga_boltz_1p_gauss = ga_boltz_1p_gauss.plotMotor;
% % ga_boltz_1p_gauss.plotGA;
%
% %%
% %
% disp('**************** Torneio + 2P + Gauss ****************');
% ga_torneio_2p_gauss = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_torneio_2p_gauss = ga_torneio_2p_gauss.Initialize;
% ga_torneio_2p_gauss = ga_torneio_2p_gauss.Solve({'Torneio', 3}, '2P', 'Gauss');
% % ga_torneio_2p_gauss = ga_torneio_2p_gauss.plotMotor;
% % ga_torneio_2p_gauss.plotGA;
%
% %
% disp('**************** Roleta + 2P + Gauss ****************');
% ga_roleta_2p_gauss = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_roleta_2p_gauss = ga_roleta_2p_gauss.Initialize;
% ga_roleta_2p_gauss = ga_roleta_2p_gauss.Solve('Roleta', '2P', 'Gauss');
% % ga_roleta_2p_gauss = ga_roleta_2p_gauss.plotMotor;
% % ga_roleta_2p_gauss.plotGA;
%
% %
% disp('**************** Rank + 2P + Gauss ****************');
% ga_rank_2p_gauss = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_rank_2p_gauss = ga_rank_2p_gauss.Initialize;
% ga_rank_2p_gauss = ga_rank_2p_gauss.Solve('Rank', '2P', 'Gauss');
% % ga_rank_2p_gauss = ga_rank_2p_gauss.plotMotor;
% % ga_rank_2p_gauss.plotGA;
%
% %
% disp('**************** Estoc + 2P + Gauss ****************');
% ga_estoc_2p_gauss = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_estoc_2p_gauss = ga_estoc_2p_gauss.Initialize;
% ga_estoc_2p_gauss = ga_estoc_2p_gauss.Solve({'Estoc',1}, '2P', 'Gauss');
% % ga_estoc_2p_gauss = ga_estoc_2p_gauss.plotMotor;
% % ga_estoc_2p_gauss.plotGA;
%
% %
% disp('**************** Boltz + 2P + Gauss ****************');
% ga_boltz_2p_gauss = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_boltz_2p_gauss = ga_boltz_2p_gauss.Initialize;
% ga_boltz_2p_gauss = ga_boltz_2p_gauss.Solve({'Boltzmann',100 , 1}, '2P', 'Gauss');
% % ga_boltz_2p_gauss = ga_boltz_2p_gauss.plotMotor;
% % ga_boltz_2p_gauss.plotGA;
%
% %%
% %
% disp('**************** Torneio + Arit + Gauss ****************');
% ga_torneio_arit_gauss = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_torneio_arit_gauss = ga_torneio_arit_gauss.Initialize;
% ga_torneio_arit_gauss = ga_torneio_arit_gauss.Solve({'Torneio', 3}, 'Arit', 'Gauss');
% % ga_torneio_arit_gauss = ga_torneio_arit_gauss.plotMotor;
% % ga_torneio_arit_gauss.plotGA;
%
% %
% disp('**************** Roleta + Arit + Gauss ****************');
% ga_roleta_arit_gauss = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_roleta_arit_gauss = ga_roleta_arit_gauss.Initialize;
% ga_roleta_arit_gauss = ga_roleta_arit_gauss.Solve('Roleta', 'Arit', 'Gauss');
% % ga_roleta_arit_gauss = ga_roleta_arit_gauss.plotMotor;
% % ga_roleta_arit_gauss.plotGA;
%
% %
% disp('**************** Rank + Arit+ Gauss ****************');
% ga_rank_arit_gauss = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_rank_arit_gauss = ga_rank_arit_gauss.Initialize;
% ga_rank_arit_gauss = ga_rank_arit_gauss.Solve('Rank', 'Arit', 'Gauss');
% % ga_rank_arit_gauss = ga_rank_arit_gauss.plotMotor;
% % ga_rank_arit_gauss.plotGA;
% 
% %
% disp('**************** Estoc + Arit + Gauss ****************');
% ga_estoc_arit_gauss = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_estoc_arit_gauss = ga_estoc_arit_gauss.Initialize;
% ga_estoc_arit_gauss = ga_estoc_arit_gauss.Solve({'Estoc',1}, 'Arit', 'Gauss');
% % ga_estoc_arit_gauss = ga_estoc_arit_gauss.plotMotor;
% % ga_estoc_arit_gauss.plotGA;
% 
% %
% disp('**************** Boltz + Arit + Gauss ****************');
% ga_boltz_arit_gauss = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_boltz_arit_gauss = ga_boltz_arit_gauss.Initialize;
% ga_boltz_arit_gauss = ga_boltz_arit_gauss.Solve({'Boltzmann',100 , 1}, 'Arit', 'Gauss');
% % ga_boltz_arit_gauss = ga_boltz_arit_gauss.plotMotor;
% % ga_boltz_arit_gauss.plotGA;
%
% %%
% %
% disp('**************** Torneio + Unif + Gauss ****************');
% ga_torneio_unif_gauss = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_torneio_unif_gauss = ga_torneio_unif_gauss.Initialize;
% ga_torneio_unif_gauss = ga_torneio_unif_gauss.Solve({'Torneio', 3}, 'Uniform', 'Gauss');
% % ga_torneio_unif_gauss = ga_torneio_unif_gauss.plotMotor;
% % ga_torneio_unif_gauss.plotGA;
%
% %
% disp('**************** Roleta + Unif + Gauss ****************');
% ga_roleta_unif_gauss = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_roleta_unif_gauss = ga_roleta_unif_gauss.Initialize;
% ga_roleta_unif_gauss = ga_roleta_unif_gauss.Solve('Roleta', 'Uniform', 'Gauss');
% % ga_roleta_unif_gauss = ga_roleta_unif_gauss.plotMotor;
% % ga_roleta_unif_gauss.plotGA;
%
% %
% disp('**************** Rank + Unif + Gauss ****************');
% ga_rank_unif_gauss = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_rank_unif_gauss = ga_rank_unif_gauss.Initialize;
% ga_rank_unif_gauss = ga_rank_unif_gauss.Solve('Rank', 'Uniform', 'Gauss');
% % ga_rank_unif_gauss = ga_rank_unif_gauss.plotMotor;
% % ga_rank_unif_gauss.plotGA;
%
% %
% disp('**************** Estoc + Unif + Gauss ****************');
% ga_estoc_unif_gauss = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_estoc_unif_gauss = ga_estoc_unif_gauss.Initialize;
% ga_estoc_unif_gauss = ga_estoc_unif_gauss.Solve({'Estoc',1}, 'Uniform', 'Gauss');
% % ga_estoc_unif_gauss = ga_estoc_unif_gauss.plotMotor;
% % ga_estoc_unif_gauss.plotGA;
%
% %
% disp('**************** Boltz + Unif + Gauss ****************');
% ga_boltz_unif_gauss = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_boltz_unif_gauss = ga_boltz_unif_gauss.Initialize;
% ga_boltz_unif_gauss = ga_boltz_unif_gauss.Solve({'Boltzmann',100 , 1}, 'Uniform', 'Gauss');
% % ga_boltz_unif_gauss = ga_boltz_unif_gauss.plotMotor;
% % ga_boltz_unif_gauss.plotGA;
%
% %%
% disp('**************** Torneio + Blend + Gauss ****************');
% ga_torneio_blend_gauss = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_torneio_blend_gauss = ga_torneio_blend_gauss.Initialize;
% ga_torneio_blend_gauss = ga_torneio_blend_gauss.Solve({'Torneio', 3}, {'Blend', 0.5}, 'Gauss');
% % ga_torneio_blend_gauss = ga_torneio_blend_gauss.plotMotor;
% % ga_torneio_blend_gauss.plotGA;
%
% disp('**************** Roleta + Blend + Gauss ****************');
% ga_roleta_blend_gauss = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_roleta_blend_gauss = ga_roleta_blend_gauss.Initialize;
% ga_roleta_blend_gauss = ga_roleta_blend_gauss.Solve('Roleta', {'Blend', 0.5}, 'Gauss');
% % ga_roleta_blend_gauss = ga_roleta_blend_gauss.plotMotor;
% % ga_roleta_blend_gauss.plotGA;
%
% %
% disp('**************** Rank + Blend + Gauss ****************');
% ga_rank_blend_gauss = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_rank_blend_gauss = ga_rank_blend_gauss.Initialize;
% ga_rank_blend_gauss = ga_rank_blend_gauss.Solve('Rank', {'Blend', 0.5}, 'Gauss');
% % ga_rank_blend_gauss = ga_rank_blend_gauss.plotMotor;
% % ga_rank_blend_gauss.plotGA;
%
% %
% disp('**************** Estoc + Blend + Gauss ****************');
% ga_estoc_blend_gauss = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_estoc_blend_gauss = ga_estoc_blend_gauss.Initialize;
% ga_estoc_blend_gauss = ga_estoc_blend_gauss.Solve({'Estoc',1}, {'Blend', 0.5}, 'Gauss');
% % ga_estoc_blend_gauss = ga_estoc_blend_gauss.plotMotor;
% % ga_estoc_blend_gauss.plotGA;
%
% %
% disp('**************** Boltz + Blend + Gauss ****************');
% ga_boltz_blend_gauss = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_boltz_blend_gauss = ga_boltz_blend_gauss.Initialize;
% ga_boltz_blend_gauss = ga_boltz_blend_gauss.Solve({'Boltzmann',100 , 1}, {'Blend', 0.5}, 'Gauss');
% % ga_boltz_blend_gauss = ga_boltz_blend_gauss.plotMotor;
% % ga_boltz_blend_gauss.plotGA;
%
%
%
%
% %% %%
% %
% disp('**************** Torneio + 1P + Uniform ****************');
% ga_torneio_1p_uniform = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_torneio_1p_uniform = ga_torneio_1p_uniform.Initialize;
% ga_torneio_1p_uniform = ga_torneio_1p_uniform.Solve({'Torneio', 3}, '1P', 'Uniform');
% % ga_torneio_1p_uniform = ga_torneio_1p_uniform.plotMotor;
% % ga_torneio_1p_uniform.plotGA;
%
% %
% disp('**************** Roleta + 1P + Uniform ****************');
% ga_roleta_1p_uniform = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_roleta_1p_uniform = ga_roleta_1p_uniform.Initialize;
% ga_roleta_1p_uniform = ga_roleta_1p_uniform.Solve('Roleta', '1P', 'Uniform');
% % ga_roleta_1p_uniform = ga_roleta_1p_uniform.plotMotor;
% % ga_roleta_1p_uniform.plotGA;
%
% %
% disp('**************** Rank + 1P + Uniform ****************');
% ga_rank_1p_uniform = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_rank_1p_uniform = ga_rank_1p_uniform.Initialize;
% ga_rank_1p_uniform = ga_rank_1p_uniform.Solve('Rank', '1P', 'Uniform');
% % ga_rank_1p_uniform = ga_rank_1p_uniform.plotMotor;
% % ga_rank_1p_uniform.plotGA;
%
% %
% disp('**************** Estoc + 1P + Uniform ****************');
% ga_estoc_1p_uniform = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_estoc_1p_uniform = ga_estoc_1p_uniform.Initialize;
% ga_estoc_1p_uniform = ga_estoc_1p_uniform.Solve({'Estoc',1} , '1P', 'Uniform');
% % ga_estoc_1p_uniform = ga_estoc_1p_uniform.plotMotor;
% % ga_estoc_1p_uniform.plotGA;
%
% %
% disp('**************** Boltz + 1P + Uniform ****************');
% ga_boltz_1p_uniform = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_boltz_1p_uniform = ga_boltz_1p_uniform.Initialize;
% ga_boltz_1p_uniform = ga_boltz_1p_uniform.Solve({'Boltzmann',100 , 1}, '1P', 'Uniform');
% % ga_boltz_1p_uniform = ga_boltz_1p_uniform.plotMotor;
% % ga_boltz_1p_uniform.plotGA;
%
% %%
% %
% disp('**************** Torneio + 2P + Uniform ****************');
% ga_torneio_2p_uniform = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_torneio_2p_uniform = ga_torneio_2p_uniform.Initialize;
% ga_torneio_2p_uniform = ga_torneio_2p_uniform.Solve({'Torneio', 3}, '2P', 'Uniform');
% % ga_torneio_2p_uniform = ga_torneio_2p_uniform.plotMotor;
% % ga_torneio_2p_uniform.plotGA;
%
% %
% disp('**************** Roleta + 2P + Uniform ****************');
% ga_roleta_2p_uniform = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_roleta_2p_uniform = ga_roleta_2p_uniform.Initialize;
% ga_roleta_2p_uniform = ga_roleta_2p_uniform.Solve('Roleta', '2P', 'Uniform');
% % ga_roleta_2p_uniform = ga_roleta_2p_uniform.plotMotor;
% % ga_roleta_2p_uniform.plotGA;
%
% %
% disp('**************** Rank + 2P + Uniform ****************');
% ga_rank_2p_uniform = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_rank_2p_uniform = ga_rank_2p_uniform.Initialize;
% ga_rank_2p_uniform = ga_rank_2p_uniform.Solve('Rank', '2P', 'Uniform');
% % ga_rank_2p_uniform = ga_rank_2p_uniform.plotMotor;
% % ga_rank_2p_uniform.plotGA;
%
% %
% disp('**************** Estoc + 2P + Uniform ****************');
% ga_estoc_2p_uniform = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_estoc_2p_uniform = ga_estoc_2p_uniform.Initialize;
% ga_estoc_2p_uniform = ga_estoc_2p_uniform.Solve({'Estoc',1}, '2P', 'Uniform');
% % ga_estoc_2p_uniform = ga_estoc_2p_uniform.plotMotor;
% % ga_estoc_2p_uniform.plotGA;
%
% %
% disp('**************** Boltz + 2P + Uniform ****************');
% ga_boltz_2p_uniform = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_boltz_2p_uniform = ga_boltz_2p_uniform.Initialize;
% ga_boltz_2p_uniform = ga_boltz_2p_uniform.Solve({'Boltzmann',100 , 1}, '2P', 'Uniform');
% % ga_boltz_2p_uniform = ga_boltz_2p_uniform.plotMotor;
% % ga_boltz_2p_uniform.plotGA;
%
% %%
% %
% disp('**************** Torneio + Arit + Uniform ****************');
% ga_torneio_arit_uniform = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_torneio_arit_uniform = ga_torneio_arit_uniform.Initialize;
% ga_torneio_arit_uniform = ga_torneio_arit_uniform.Solve({'Torneio', 3}, 'Arit', 'Uniform');
% % ga_torneio_arit_uniform = ga_torneio_arit_uniform.plotMotor;
% % ga_torneio_arit_uniform.plotGA;
%
% % %
% disp('**************** Roleta + Arit + Uniform ****************');
% ga_roleta_arit_uniform = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_roleta_arit_uniform = ga_roleta_arit_uniform.Initialize;
% ga_roleta_arit_uniform = ga_roleta_arit_uniform.Solve('Roleta', 'Arit', 'Uniform');
% % ga_roleta_arit_uniform = ga_roleta_arit_uniform.plotMotor;
% % ga_roleta_arit_uniform.plotGA;
%
% %
% disp('**************** Rank + Arit+ Uniform ****************');
% ga_rank_arit_uniform = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_rank_arit_uniform = ga_rank_arit_uniform.Initialize;
% ga_rank_arit_uniform = ga_rank_arit_uniform.Solve('Rank', 'Arit', 'Uniform');
% % ga_rank_arit_uniform = ga_rank_arit_uniform.plotMotor;
% % ga_rank_arit_uniform.plotGA;
%
% %
% disp('**************** Estoc + Arit + Uniform ****************');
% ga_estoc_arit_uniform = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_estoc_arit_uniform = ga_estoc_arit_uniform.Initialize;
% ga_estoc_arit_uniform = ga_estoc_arit_uniform.Solve({'Estoc',1}, 'Arit', 'Uniform');
% % ga_estoc_arit_uniform = ga_estoc_arit_uniform.plotMotor;
% % ga_estoc_arit_uniform.plotGA;
%
% %
% % disp('**************** Boltz + Arit + Uniform ****************');
% % ga_boltz_arit_uniform = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% % ga_boltz_arit_uniform = ga_boltz_arit_uniform.Initialize;
% % ga_boltz_arit_uniform = ga_boltz_arit_uniform.Solve({'Boltzmann',[100 , 1]}, 'Arit', 'Uniform');
% % ga_boltz_arit_uniform = ga_boltz_arit_uniform.plotMotor;
% % ga_boltz_arit_uniform.plotGA;
%
% %%
% %
% disp('**************** Torneio + Unif + Uniform ****************');
% ga_torneio_unif_uniform = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_torneio_unif_uniform = ga_torneio_unif_uniform.Initialize;
% ga_torneio_unif_uniform = ga_torneio_unif_uniform.Solve({'Torneio', 3}, 'Uniform', 'Uniform');
% % ga_torneio_unif_uniform = ga_torneio_unif_uniform.plotMotor;
% % ga_torneio_unif_uniform.plotGA;
%
% %
% disp('**************** Roleta + Unif + Uniform ****************');
% ga_roleta_unif_uniform = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_roleta_unif_uniform = ga_roleta_unif_uniform.Initialize;
% ga_roleta_unif_uniform = ga_roleta_unif_uniform.Solve('Roleta', 'Uniform', 'Uniform');
% % ga_roleta_unif_uniform = ga_roleta_unif_uniform.plotMotor;
% % ga_roleta_unif_uniform.plotGA;
%
% %
% disp('**************** Rank + Unif + Uniform ****************');
% ga_rank_unif_uniform = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_rank_unif_uniform = ga_rank_unif_uniform.Initialize;
% ga_rank_unif_uniform = ga_rank_unif_uniform.Solve('Rank', 'Uniform', 'Uniform');
% % ga_rank_unif_uniform = ga_rank_unif_uniform.plotMotor;
% % ga_rank_unif_uniform.plotGA;
%
% %
% disp('**************** Estoc + Unif + Uniform ****************');
% ga_estoc_unif_uniform = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_estoc_unif_uniform = ga_estoc_unif_uniform.Initialize;
% ga_estoc_unif_uniform = ga_estoc_unif_uniform.Solve({'Estoc',1}, 'Uniform', 'Uniform');
% % ga_estoc_unif_uniform = ga_estoc_unif_uniform.plotMotor;
% % ga_estoc_unif_uniform.plotGA;
%
% %
% disp('**************** Boltz + Unif + Uniform ****************');
% ga_boltz_unif_uniform = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_boltz_unif_uniform = ga_boltz_unif_uniform.Initialize;
% ga_boltz_unif_uniform = ga_boltz_unif_uniform.Solve({'Boltzmann',100 , 1}, 'Uniform', 'Uniform');
% % ga_boltz_unif_uniform = ga_boltz_unif_uniform.plotMotor;
% % ga_boltz_unif_uniform.plotGA;
%
% %%
% disp('**************** Torneio + Blend + Uniform ****************');
% ga_torneio_blend_uniform = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_torneio_blend_uniform = ga_torneio_blend_uniform.Initialize;
% ga_torneio_blend_uniform = ga_torneio_blend_uniform.Solve({'Torneio', 3}, {'Blend', 0.5}, 'Uniform');
% % ga_torneio_blend_uniform = ga_torneio_blend_uniform.plotMotor;
% % ga_torneio_blend_uniform.plotGA;
%
% disp('**************** Roleta + Blend + Uniform ****************');
% ga_roleta_blend_uniform = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_roleta_blend_uniform = ga_roleta_blend_uniform.Initialize;
% ga_roleta_blend_uniform = ga_roleta_blend_uniform.Solve('Roleta', {'Blend', 0.5}, 'Uniform');
% % ga_roleta_blend_uniform = ga_roleta_blend_uniform.plotMotor;
% % ga_roleta_blend_uniform.plotGA;
%
% %
% disp('**************** Rank + Blend + Uniform ****************');
% ga_rank_blend_uniform = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_rank_blend_uniform = ga_rank_blend_uniform.Initialize;
% ga_rank_blend_uniform = ga_rank_blend_uniform.Solve('Rank', {'Blend', 0.5}, 'Uniform');
% % ga_rank_blend_uniform = ga_rank_blend_uniform.plotMotor;
% % ga_rank_blend_uniform.plotGA;
%
% %
disp('**************** Estoc + Blend + Uniform ****************');
ga_estoc_blend_uniform = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
ga_estoc_blend_uniform = ga_estoc_blend_uniform.Initialize;
ga_estoc_blend_uniform = ga_estoc_blend_uniform.Solve({'Estoc',2}, {'Blend', 0.5}, 'Uniform');
% ga_estoc_blend_uniform = ga_estoc_blend_uniform.plotMotor;
% ga_estoc_blend_uniform.plotGA;
%
% %
% disp('**************** Boltz + Blend + Uniform ****************');
% ga_boltz_blend_uniform = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_boltz_blend_uniform = ga_boltz_blend_uniform.Initialize;
% ga_boltz_blend_uniform = ga_boltz_blend_uniform.Solve({'Boltzmann',100 , 1}, {'Blend', 0.5}, 'Uniform');
% % ga_boltz_blend_uniform = ga_boltz_blend_uniform.plotMotor;
% % ga_boltz_blend_uniform.plotGA;
%
%
%
%
% %% %%
% %
% disp('**************** Torneio + 1P + Creep ****************');
% ga_torneio_1p_creep = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_torneio_1p_creep = ga_torneio_1p_creep.Initialize;
% ga_torneio_1p_creep = ga_torneio_1p_creep.Solve({'Torneio', 3}, '1P', 'Creep');
% % ga_torneio_1p_creep = ga_torneio_1p_creep.plotMotor;
% % ga_torneio_1p_creep.plotGA;
%
% %
% disp('**************** Roleta + 1P + Creep ****************');
% ga_roleta_1p_creep = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_roleta_1p_creep = ga_roleta_1p_creep.Initialize;
% ga_roleta_1p_creep = ga_roleta_1p_creep.Solve('Roleta', '1P', 'Creep');
% % ga_roleta_1p_creep = ga_roleta_1p_creep.plotMotor;
% % ga_roleta_1p_creep.plotGA;
%
% %
% disp('**************** Rank + 1P + Creep ****************');
% ga_rank_1p_creep = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_rank_1p_creep = ga_rank_1p_creep.Initialize;
% ga_rank_1p_creep = ga_rank_1p_creep.Solve('Rank', '1P', 'Creep');
% % ga_rank_1p_creep = ga_rank_1p_creep.plotMotor;
% % ga_rank_1p_creep.plotGA;
%
% %
% disp('**************** Estoc + 1P + Creep ****************');
% ga_estoc_1p_creep = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_estoc_1p_creep = ga_estoc_1p_creep.Initialize;
% ga_estoc_1p_creep = ga_estoc_1p_creep.Solve({'Estoc',1} , '1P', 'Creep');
% % ga_estoc_1p_creep = ga_estoc_1p_creep.plotMotor;
% % ga_estoc_1p_creep.plotGA;
%
% %
% disp('**************** Boltz + 1P + Creep ****************');
% ga_boltz_1p_creep = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_boltz_1p_creep = ga_boltz_1p_creep.Initialize;
% ga_boltz_1p_creep = ga_boltz_1p_creep.Solve({'Boltzmann',100 , 1}, '1P', 'Creep');
% % ga_boltz_1p_creep = ga_boltz_1p_creep.plotMotor;
% % ga_boltz_1p_creep.plotGA;
%
% %%
% %
% disp('**************** Torneio + 2P + Creep ****************');
% ga_torneio_2p_creep = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_torneio_2p_creep = ga_torneio_2p_creep.Initialize;
% ga_torneio_2p_creep = ga_torneio_2p_creep.Solve({'Torneio', 3}, '2P', 'Creep');
% % ga_torneio_2p_creep = ga_torneio_2p_creep.plotMotor;
% % ga_torneio_2p_creep.plotGA;
%
% % %
% disp('**************** Roleta + 2P + Creep ****************');
% ga_roleta_2p_creep = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_roleta_2p_creep = ga_roleta_2p_creep.Initialize;
% ga_roleta_2p_creep = ga_roleta_2p_creep.Solve('Roleta', '2P', 'Creep');
% ga_roleta_2p_creep = ga_roleta_2p_creep.plotMotor;
% ga_roleta_2p_creep.plotGA;

 
% disp('**************** Rank + 2P + Creep ****************');
% ga_rank_2p_creep = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_rank_2p_creep = ga_rank_2p_creep.Initialize;
% ga_rank_2p_creep = ga_rank_2p_creep.Solve('Rank', '2P', 'Creep');
% % ga_rank_2p_creep = ga_rank_2p_creep.plotMotor;
% % ga_rank_2p_creep.plotGA;
%
% %
% disp('**************** Estoc + 2P + Creep ****************');
% ga_estoc_2p_creep = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_estoc_2p_creep = ga_estoc_2p_creep.Initialize;
% ga_estoc_2p_creep = ga_estoc_2p_creep.Solve({'Estoc',1}, '2P', 'Creep');
% % ga_estoc_2p_creep = ga_estoc_2p_creep.plotMotor;
% % ga_estoc_2p_creep.plotGA;
%
% %
% disp('**************** Boltz + 2P + Creep ****************');
% ga_boltz_2p_creep = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_boltz_2p_creep = ga_boltz_2p_creep.Initialize;
% ga_boltz_2p_creep = ga_boltz_2p_creep.Solve({'Boltzmann',100 , 1}, '2P', 'Creep');
% % ga_boltz_2p_creep = ga_boltz_2p_creep.plotMotor;
% % ga_boltz_2p_creep.plotGA;
%
% %%
% %
% disp('**************** Torneio + Arit + Creep ****************');
% ga_torneio_arit_creep = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_torneio_arit_creep = ga_torneio_arit_creep.Initialize;
% ga_torneio_arit_creep = ga_torneio_arit_creep.Solve({'Torneio', 3}, 'Arit', 'Creep');
% % ga_torneio_arit_creep = ga_torneio_arit_creep.plotMotor;
% % ga_torneio_arit_creep.plotGA;
%
% %
% disp('**************** Roleta + Arit + Creep ****************');
% ga_roleta_arit_creep = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_roleta_arit_creep = ga_roleta_arit_creep.Initialize;
% ga_roleta_arit_creep = ga_roleta_arit_creep.Solve('Roleta', 'Arit', 'Creep');
% % ga_roleta_arit_creep = ga_roleta_arit_creep.plotMotor;
% ga_roleta_arit_creep.plotGA;
%
% %
% disp('**************** Rank + Arit+ Creep ****************');
% ga_rank_arit_creep = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_rank_arit_creep = ga_rank_arit_creep.Initialize;
% ga_rank_arit_creep = ga_rank_arit_creep.Solve('Rank', 'Arit', 'Creep');
% % ga_rank_arit_creep = ga_rank_arit_creep.plotMotor;
% % ga_rank_arit_creep.plotGA;
%
% %
% disp('**************** Estoc + Arit + Creep ****************');
% ga_estoc_arit_creep = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_estoc_arit_creep = ga_estoc_arit_creep.Initialize;
% ga_estoc_arit_creep = ga_estoc_arit_creep.Solve({'Estoc',1}, 'Arit', 'Creep');
% % ga_estoc_arit_creep = ga_estoc_arit_creep.plotMotor;
% % ga_estoc_arit_creep.plotGA;
%
% %
% disp('**************** Boltz + Arit + Creep ****************');
% ga_boltz_arit_creep = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_boltz_arit_creep = ga_boltz_arit_creep.Initialize;
% ga_boltz_arit_creep = ga_boltz_arit_creep.Solve({'Boltzmann',100 , 1}, 'Arit', 'Creep');
% % ga_boltz_arit_creep = ga_boltz_arit_creep.plotMotor;
% % ga_boltz_arit_creep.plotGA;
%
% %%
% %
% disp('**************** Torneio + Unif + Creep ****************');
% ga_torneio_unif_creep = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_torneio_unif_creep = ga_torneio_unif_creep.Initialize;
% ga_torneio_unif_creep = ga_torneio_unif_creep.Solve({'Torneio', 3}, 'Uniform', 'Creep');
% % ga_torneio_unif_creep = ga_torneio_unif_creep.plotMotor;
% % ga_torneio_unif_creep.plotGA;
%
% %
% disp('**************** Roleta + Unif + Creep ****************');
% ga_roleta_unif_creep = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_roleta_unif_creep = ga_roleta_unif_creep.Initialize;
% ga_roleta_unif_creep = ga_roleta_unif_creep.Solve('Roleta', 'Uniform', 'Creep');
% % ga_roleta_unif_creep = ga_roleta_unif_creep.plotMotor;
% % ga_roleta_unif_creep.plotGA;
%
% %
% disp('**************** Rank + Unif + Creep ****************');
% ga_rank_unif_creep = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_rank_unif_creep = ga_rank_unif_creep.Initialize;
% ga_rank_unif_creep = ga_rank_unif_creep.Solve('Rank', 'Uniform', 'Creep');
% % ga_rank_unif_creep = ga_rank_unif_creep.plotMotor;
% % ga_rank_unif_creep.plotGA;
%
% %
% disp('**************** Estoc + Unif + Creep ****************');
% ga_estoc_unif_creep = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_estoc_unif_creep = ga_estoc_unif_creep.Initialize;
% ga_estoc_unif_creep = ga_estoc_unif_creep.Solve({'Estoc',1}, 'Uniform', 'Creep');
% % ga_estoc_unif_creep = ga_estoc_unif_creep.plotMotor;
% % ga_estoc_unif_creep.plotGA;
%
% %
% disp('**************** Boltz + Unif + Creep ****************');
% ga_boltz_unif_creep = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_boltz_unif_creep = ga_boltz_unif_creep.Initialize;
% ga_boltz_unif_creep = ga_boltz_unif_creep.Solve({'Boltzmann',100 , 1}, 'Uniform', 'Creep');
% % ga_boltz_unif_creep = ga_boltz_unif_creep.plotMotor;
% % ga_boltz_unif_creep.plotGA;
%
% %%
% disp('**************** Torneio + Blend + Creep ****************');
% ga_torneio_blend_creep = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_torneio_blend_creep = ga_torneio_blend_creep.Initialize;
% ga_torneio_blend_creep = ga_torneio_blend_creep.Solve({'Torneio', 3}, {'Blend', 0.5}, 'Creep');
% % ga_torneio_blend_creep = ga_torneio_blend_creep.plotMotor;
% % ga_torneio_blend_creep.plotGA;
%
% disp('**************** Roleta + Blend + Creep ****************');
% ga_roleta_blend_creep = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_roleta_blend_creep = ga_roleta_blend_creep.Initialize;
% ga_roleta_blend_creep = ga_roleta_blend_creep.Solve('Roleta', {'Blend', 0.5}, 'Creep');
% % ga_roleta_blend_creep = ga_roleta_blend_creep.plotMotor;
% % ga_roleta_blend_creep.plotGA;
%
% %
% disp('**************** Rank + Blend + Creep ****************');
% ga_rank_blend_creep = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_rank_blend_creep = ga_rank_blend_creep.Initialize;
% ga_rank_blend_creep = ga_rank_blend_creep.Solve('Rank', {'Blend', 0.5}, 'Creep');
% % ga_rank_blend_creep = ga_rank_blend_creep.plotMotor;
% % ga_rank_blend_creep.plotGA;
%
% %
% disp('**************** Estoc + Blend + Creep ****************');
% ga_estoc_blend_creep = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_estoc_blend_creep = ga_estoc_blend_creep.Initialize;
% ga_estoc_blend_creep = ga_estoc_blend_creep.Solve({'Estoc',1}, {'Blend', 0.5}, 'Creep');
% % ga_estoc_blend_creep = ga_estoc_blend_creep.plotMotor;
% % ga_estoc_blend_creep.plotGA;
%
% %
% disp('**************** Boltz + Blend + Creep ****************');
% ga_boltz_blend_creep = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_boltz_blend_creep = ga_boltz_blend_creep.Initialize;
% ga_boltz_blend_creep = ga_boltz_blend_creep.Solve({'Boltzmann',100 , 1}, {'Blend', 0.5}, 'Creep');
% % ga_boltz_blend_creep = ga_boltz_blend_creep.plotMotor;
% % ga_boltz_blend_creep.plotGA;
%
%
%
%
% % %%
% %
% disp('**************** Torneio + 1P + Inv ****************');
% ga_torneio_1p_inv = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_torneio_1p_inv = ga_torneio_1p_inv.Initialize;
% ga_torneio_1p_inv = ga_torneio_1p_inv.Solve({'Torneio', 3}, '1P', 'Inv');
% % ga_torneio_1p_inv = ga_torneio_1p_inv.plotMotor;
% % ga_torneio_1p_inv.plotGA;
%
% %
% disp('**************** Roleta + 1P + Inv ****************');
% ga_roleta_1p_inv = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_roleta_1p_inv = ga_roleta_1p_inv.Initialize;
% ga_roleta_1p_inv = ga_roleta_1p_inv.Solve('Roleta', '1P', 'Inv');
% % ga_roleta_1p_inv = ga_roleta_1p_inv.plotMotor;
% % ga_roleta_1p_inv.plotGA;
%
% %
% disp('**************** Rank + 1P + Inv ****************');
% ga_rank_1p_inv = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_rank_1p_inv = ga_rank_1p_inv.Initialize;
% ga_rank_1p_inv = ga_rank_1p_inv.Solve('Rank', '1P', 'Inv');
% % ga_rank_1p_inv = ga_rank_1p_inv.plotMotor;
% % ga_rank_1p_inv.plotGA;
%
% %
% disp('**************** Estoc + 1P + Inv ****************');
% ga_estoc_1p_inv = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_estoc_1p_inv = ga_estoc_1p_inv.Initialize;
% ga_estoc_1p_inv = ga_estoc_1p_inv.Solve({'Estoc',1} , '1P', 'Inv');
% % ga_estoc_1p_inv = ga_estoc_1p_inv.plotMotor;
% % ga_estoc_1p_inv.plotGA;
%
% %
% disp('**************** Boltz + 1P + Inv ****************');
% ga_boltz_1p_inv = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_boltz_1p_inv = ga_boltz_1p_inv.Initialize;
% ga_boltz_1p_inv = ga_boltz_1p_inv.Solve({'Boltzmann',100 , 1}, '1P', 'Inv');
% % ga_boltz_1p_inv = ga_boltz_1p_inv.plotMotor;
% % ga_boltz_1p_inv.plotGA;
%
% %%
% %
% disp('**************** Torneio + 2P + Inv ****************');
% ga_torneio_2p_inv = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_torneio_2p_inv = ga_torneio_2p_inv.Initialize;
% ga_torneio_2p_inv = ga_torneio_2p_inv.Solve({'Torneio', 3}, '2P', 'Inv');
% % ga_torneio_2p_inv = ga_torneio_2p_inv.plotMotor;
% % ga_torneio_2p_inv.plotGA;
%
% %
% disp('**************** Roleta + 2P + Inv ****************');
% ga_roleta_2p_inv = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_roleta_2p_inv = ga_roleta_2p_inv.Initialize;
% ga_roleta_2p_inv = ga_roleta_2p_inv.Solve('Roleta', '2P', 'Inv');
% % ga_roleta_2p_inv = ga_roleta_2p_inv.plotMotor;
% % ga_roleta_2p_inv.plotGA;
%
% %
% disp('**************** Rank + 2P + Inv ****************');
% ga_rank_2p_inv = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_rank_2p_inv = ga_rank_2p_inv.Initialize;
% ga_rank_2p_inv = ga_rank_2p_inv.Solve('Rank', '2P', 'Inv');
% % ga_rank_2p_inv = ga_rank_2p_inv.plotMotor;
% % ga_rank_2p_inv.plotGA;
%
% %
% disp('**************** Estoc + 2P + Inv ****************');
% ga_estoc_2p_inv = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_estoc_2p_inv = ga_estoc_2p_inv.Initialize;
% ga_estoc_2p_inv = ga_estoc_2p_inv.Solve({'Estoc',1}, '2P', 'Inv');
% % ga_estoc_2p_inv = ga_estoc_2p_inv.plotMotor;
% % ga_estoc_2p_inv.plotGA;
%
% %
% disp('**************** Boltz + 2P + Inv ****************');
% ga_boltz_2p_inv = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_boltz_2p_inv = ga_boltz_2p_inv.Initialize;
% ga_boltz_2p_inv = ga_boltz_2p_inv.Solve({'Boltzmann',100 , 1}, '2P', 'Inv');
% % ga_boltz_2p_inv = ga_boltz_2p_inv.plotMotor;
% % ga_boltz_2p_inv.plotGA;
%
% %%
% %
% disp('**************** Torneio + Arit + Inv ****************');
% ga_torneio_arit_inv = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_torneio_arit_inv = ga_torneio_arit_inv.Initialize;
% ga_torneio_arit_inv = ga_torneio_arit_inv.Solve({'Torneio', 3}, 'Arit', 'Inv');
% % ga_torneio_arit_inv = ga_torneio_arit_inv.plotMotor;
% % ga_torneio_arit_inv.plotGA;
%
% % %
% disp('**************** Roleta + Arit + Inv ****************');
% ga_roleta_arit_inv = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_roleta_arit_inv = ga_roleta_arit_inv.Initialize;
% ga_roleta_arit_inv = ga_roleta_arit_inv.Solve('Roleta', 'Arit', 'Inv');
% % ga_roleta_arit_inv = ga_roleta_arit_inv.plotMotor;
% % ga_roleta_arit_inv.plotGA;

% %
% disp('**************** Rank + Arit + Inv ****************');
% ga_rank_arit_inv = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_rank_arit_inv = ga_rank_arit_inv.Initialize;
% ga_rank_arit_inv = ga_rank_arit_inv.Solve('Rank', 'Arit', 'Inv');
% % ga_rank_arit_inv = ga_rank_arit_inv.plotMotor;
% % ga_rank_arit_inv.plotGA;
%
% %
% disp('**************** Estoc + Arit + Inv ****************');
% ga_estoc_arit_inv = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_estoc_arit_inv = ga_estoc_arit_inv.Initialize;
% ga_estoc_arit_inv = ga_estoc_arit_inv.Solve({'Estoc',1}, 'Arit', 'Inv');
% % ga_estoc_arit_inv = ga_estoc_arit_inv.plotMotor;
% % ga_estoc_arit_inv.plotGA;
%
% %
% disp('**************** Boltz + Arit + Inv ****************');
% ga_boltz_arit_inv = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_boltz_arit_inv = ga_boltz_arit_inv.Initialize;
% ga_boltz_arit_inv = ga_boltz_arit_inv.Solve({'Boltzmann',100 , 1}, 'Arit', 'Inv');
% % ga_boltz_arit_inv = ga_boltz_arit_inv.plotMotor;
% % ga_boltz_arit_inv.plotGA;
%
% %%
% %
% disp('**************** Torneio + Unif + Inv ****************');
% ga_torneio_unif_inv = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_torneio_unif_inv = ga_torneio_unif_inv.Initialize;
% ga_torneio_unif_inv = ga_torneio_unif_inv.Solve({'Torneio', 3}, 'Uniform', 'Inv');
% % ga_torneio_unif_inv = ga_torneio_unif_inv.plotMotor;
% % ga_torneio_unif_inv.plotGA;
%
% %
% disp('**************** Roleta + Unif + Inv ****************');
% ga_roleta_unif_inv = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_roleta_unif_inv = ga_roleta_unif_inv.Initialize;
% ga_roleta_unif_inv = ga_roleta_unif_inv.Solve('Roleta', 'Uniform', 'Inv');
% % ga_roleta_unif_inv = ga_roleta_unif_inv.plotMotor;
% % ga_roleta_unif_inv.plotGA;
%
% %
% disp('**************** Rank + Unif + Inv ****************');
% ga_rank_unif_inv = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_rank_unif_inv = ga_rank_unif_inv.Initialize;
% ga_rank_unif_inv = ga_rank_unif_inv.Solve('Rank', 'Uniform', 'Inv');
% % ga_rank_unif_inv = ga_rank_unif_inv.plotMotor;
% % ga_rank_unif_inv.plotGA;
%
% %
% disp('**************** Estoc + Unif + Inv ****************');
% ga_estoc_unif_inv = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_estoc_unif_inv = ga_estoc_unif_inv.Initialize;
% ga_estoc_unif_inv = ga_estoc_unif_inv.Solve({'Estoc',1}, 'Uniform', 'Inv');
% % ga_estoc_unif_inv = ga_estoc_unif_inv.plotMotor;
% % ga_estoc_unif_inv.plotGA;
%
% %
% disp('**************** Boltz + Unif + Inv ****************');
% ga_boltz_unif_inv = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_boltz_unif_inv = ga_boltz_unif_inv.Initialize;
% ga_boltz_unif_inv = ga_boltz_unif_inv.Solve({'Boltzmann',100 , 1}, 'Uniform', 'Inv');
% % ga_boltz_unif_inv = ga_boltz_unif_inv.plotMotor;
% % ga_boltz_unif_inv.plotGA;
%
% %%
% disp('**************** Torneio + Blend + Inv ****************');
% ga_torneio_blend_inv = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_torneio_blend_inv = ga_torneio_blend_inv.Initialize;
% ga_torneio_blend_inv = ga_torneio_blend_inv.Solve({'Torneio', 3}, {'Blend', 0.5}, 'Inv');
% % ga_torneio_blend_inv = ga_torneio_blend_inv.plotMotor;
% % ga_torneio_blend_inv.plotGA;
%
% disp('**************** Roleta + Blend + Inv ****************');
% ga_roleta_blend_inv = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_roleta_blend_inv = ga_roleta_blend_inv.Initialize;
% ga_roleta_blend_inv = ga_roleta_blend_inv.Solve('Roleta', {'Blend', 0.5}, 'Inv');
% % ga_roleta_blend_inv = ga_roleta_blend_inv.plotMotor;
% % ga_roleta_blend_inv.plotGA;
%
% %
% disp('**************** Rank + Blend + Inv ****************');
% ga_rank_blend_inv = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_rank_blend_inv = ga_rank_blend_inv.Initialize;
% ga_rank_blend_inv = ga_rank_blend_inv.Solve('Rank', {'Blend', 0.5}, 'Inv');
% % ga_rank_blend_inv = ga_rank_blend_inv.plotMotor;
% % ga_rank_blend_inv.plotGA;
%
% %
% disp('**************** Estoc + Blend + Inv ****************');
% ga_estoc_blend_inv = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_estoc_blend_inv = ga_estoc_blend_inv.Initialize;
% ga_estoc_blend_inv = ga_estoc_blend_inv.Solve({'Estoc',1}, {'Blend', 0.5}, 'Inv');
% % ga_estoc_blend_inv = ga_estoc_blend_inv.plotMotor;
% % ga_estoc_blend_inv.plotGA;
%
% %
% disp('**************** Boltz + Blend + Inv ****************');
% ga_boltz_blend_inv = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
% ga_boltz_blend_inv = ga_boltz_blend_inv.Initialize;
% ga_boltz_blend_inv = ga_boltz_blend_inv.Solve({'Boltzmann',100 , 1}, {'Blend', 0.5}, 'Inv');
% % ga_boltz_blend_inv = ga_boltz_blend_inv.plotMotor;
% % ga_boltz_blend_inv.plotGA;

% Número de indivíduos na população
popSize = 100;

% Número de gerações
numGenerations = 100;

% Taxa de crossover
crossoverRate = 1;

% Taxa de mutação
mutationRate = 0.15;

% Número de Parâmetros
paramNum = 4;

% Limites dos parâmetros [Rs, Rr, Lm, Lls, Llr]
paramBounds = [0.1, 5; % Rs (0.6759)
    0.1, 5; % Rr (0.2615)
    0.001, 0.1; % Lm (0.0387)
    0.0001, 0.1]; % Lls (0.0028)

% func = 'f_motor_high_data';
func = 'f_motor_low_data';

disp('**************** Estoc + Blend + Uniform ****************');
ga_estoc_blend_uniform = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func);
ga_estoc_blend_uniform = ga_estoc_blend_uniform.Initialize;
ga_estoc_blend_uniform = ga_estoc_blend_uniform.Solve({'Estoc',2}, {'Blend', 0.5}, 'Uniform');
% ga_estoc_blend_uniform = ga_estoc_blend_uniform.plotMotor;
% ga_estoc_blend_uniform.plotGA;

func = 'f_motor_high_data_complement';
% Número de Parâmetros
paramNum = 1;

% Número de indivíduos na população
popSize = 100;

% Número de gerações
numGenerations = 10;

% Limites dos parâmetros [Rs, Rr, Lm, Lls, Llr]
paramBounds = [0.0001, 0.1]; % Lls (0.0028)

Rs = ga_estoc_blend_uniform.bestSolution(1);
Rr = ga_estoc_blend_uniform.bestSolution(2);
Lm = ga_estoc_blend_uniform.bestSolution(3);

disp('**************** Estoc + Blend + Creep ****************');
ga_estoc_blend_creep = GA(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, func, Rs, Rr, Lm);
ga_estoc_blend_creep = ga_estoc_blend_creep.Initialize;
ga_estoc_blend_creep = ga_estoc_blend_creep.Solve({'Estoc',2}, {'Blend', 0.5}, 'Creep');
ga_estoc_blend_creep = ga_estoc_blend_creep.plotMotorComplement;
% ga_estoc_blend_creep.plotGA;
