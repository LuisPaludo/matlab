clear all;
close all;
clear;

% Número de indivíduos na população
popSize = 100;

% Número de gerações
numGenerations = 100;

% Taxa de crossover
crossoverRate = 1;

% Taxa de mutação
mutationRate = 0.15;

% ParamNum
paramNum = 6;

% Limites dos parâmetros
paramBounds = [0, 10000; % KP WR
    0, 10000; % KI WR
    0, 10000; % KP IQ/ID
    0, 10000;
    0, 10000; % KP IQ/ID
    0, 10000];

% Definindo os parâmetros
selectionMethods = {'Torneio', 'Roleta', 'Estoc'};
crossoverMethods = {'1P', '2P', 'Uniform', 'Blend'};
mutationMethods = {'Gauss', 'Uniform', 'Creep'};

% Argumentos adicionais para os métodos de seleção
selectionArgsArray = {3, [], 2}; % Ajuste conforme necessário

numExecucoes = 5;

% Arrays para armazenar os dados acumulados
allBestFitness = zeros(numExecucoes, length(selectionMethods) * length(crossoverMethods) * length(mutationMethods));
allExecutionTimes = zeros(size(allBestFitness));
varNames = cell(1, size(allBestFitness, 2));

Rs = 0.6759;
Rr = 0.2615;
Lm = 0.0387;
Lls = 0.00280;
buffer = 10;

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
                gaObj = GAPI_2(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, Rs, Rr, Lm, Lls, buffer);
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
                gaObj.plotMotor;
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
