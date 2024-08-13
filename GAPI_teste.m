% Número de indivíduos na população
popSize = 100;

% Número de gerações
numGenerations = 100;

% Taxa de crossover
crossoverRate = 1;

% Taxa de mutação
mutationRate = 0.15;

% Limites dos parâmetros
paramBounds = [0, 100000; % KP WR
    0, 100000; % KI WR
    0, 100000; % KP IQ/ID
    0, 100000;
    0, 100000; % KP IQ/ID
    0, 100000];

func = 'f_motor_low_data';

% Definindo os parâmetros
selectionMethods = {'Torneio', 'Roleta', 'Estoc', 'Boltzmann'};
crossoverMethods = {'1P', '2P', 'Uniform', 'Blend'};
mutationMethods = {'Gauss', 'Uniform', 'Creep'};

% Argumentos adicionais para os métodos de seleção
selectionArgsArray = {3, [], 2, [100, 1]}; % Ajuste conforme necessário

numExecucoes = 5;

% Arrays para armazenar os dados acumulados
allBestFitness = zeros(numExecucoes, length(selectionMethods) * length(crossoverMethods) * length(mutationMethods));
allExecutionTimes = zeros(size(allBestFitness));
varNames = cell(1, size(allBestFitness, 2));

Rs = 0.6759;
Rr = 0.2615;
Lm = 0.0387;
Lls = 0.00280;


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
                gaObj = GAPI(popSize, numGenerations, crossoverRate, mutationRate, paramBounds, Rs, Rr, Lm, Lls, 1);
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
