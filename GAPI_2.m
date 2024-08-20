classdef GAPI_2
    properties
        popSize {mustBePositive, mustBeInteger}
        generation {mustBeInteger}
        numGenerations {mustBePositive, mustBeInteger}
        crossoverRate {mustBePositive, mustBeFloat}
        mutationRate {mustBePositive, mustBeFloat}
        paramNum {mustBePositive, mustBeInteger}
        pop
        paramBounds
        lowerBounds
        upperBounds
        CostFunction
        fitnessScores
        bestFitnessHistory
        averageFitnessHistory
        bestIndividuals
        currentBestFitness
        bestFitness
        parent1
        parent2
        bestSolution
        wr_vetor
        iqs_vetor
        ids_vetor
        te_vetor
        t_vetor
        executionTime
        Rs
        Rr
        Lm
        Lls
        buffer;
    end

    methods
        % Construtor da classe
        function obj = GAPI_2(popSize, numGenerations, crossoverRate, mutationRate, paramNum, paramBounds, Rs, Rr, Lm, Lls, buffer)
            obj.popSize = popSize;
            obj.numGenerations = numGenerations;
            obj.crossoverRate = crossoverRate;
            obj.mutationRate = mutationRate;
            obj.paramNum = paramNum;
            obj.pop = zeros(popSize, paramNum);
            obj.paramBounds = paramBounds;
            obj.lowerBounds = repmat(paramBounds(:,1)', popSize, 1);
            obj.upperBounds = repmat(paramBounds(:,2)', popSize, 1);
            obj.fitnessScores = zeros(popSize, 1);
            obj.bestFitnessHistory = zeros(numGenerations, 1);
            obj.averageFitnessHistory = zeros(numGenerations, 1);
            obj.bestIndividuals = zeros(numGenerations, length(paramBounds));
            obj.generation = 0;

            obj.buffer = buffer;

            obj.Lm = Lm;
            obj.Rr = Rr;
            obj.Rs = Rs;
            obj.Lls = Lls;

            obj.CostFunction = @obj.f_motor;

        end

        % Inicializador
        function obj = Initialize(obj)

            % Gerar população inicial
            obj.pop = rand(obj.popSize, obj.paramNum) .* (obj.upperBounds - obj.lowerBounds) + obj.lowerBounds;
            % Inicializar vetor de aptidão
            for i = 1:obj.popSize
                obj.fitnessScores(i) = obj.CostFunction(obj.pop(i,:));
            end
        end

        %% Solvers
        % Solver
        function obj = Solve(obj, selectionArgs, crossoverArgs, mutation)
            tic;
            % Lógica de crossover
            if iscell(crossoverArgs)
                crossover = crossoverArgs{1};
                % additionalParams = crossoverArgs{2:end};
            else
                crossover = crossoverArgs;
                % additionalParams = [];
            end

            % Lógica de crossover
            if iscell(selectionArgs)
                selection = selectionArgs{1};
            else
                selection = selectionArgs;
            end
            % Loop Principal
            if strcmp(selection, 'Torneio')
                tournamentSize = selectionArgs{2};
            end
            if strcmp(selection, 'Estoc')
                estocSize  = selectionArgs{2};
            end
            if strcmp(selection, '')
                disp('Digite o tipo de seleção')
                return
            end

            if strcmp(crossover, 'Blend')
                % disp('Digite o valor de alpha')
                alpha = crossoverArgs{2};
            end
            if strcmp(crossover, '')
                disp('Selecione um método de crossover')
                return
            end

            if strcmp(mutation, '')
                disp('Selecione um método de mutação')
                return
            end

            for gen = 1:obj.numGenerations
                obj.generation = gen;
                % Aqui, uma nova matriz newPop é criada para armazenar a nova geração de indivíduos.
                % Ela tem o mesmo tamanho que a população atual (pop).
                newPop = zeros(size(obj.pop));  % Inicializar a nova população

                % Este loop itera sobre o tamanho da população (popSize). Em cada iteração, ele
                % seleciona dois pais, aplica o crossover e a mutação para criar um novo indivíduo, e então adiciona este indivíduo à newPop.
                for i = 1:obj.popSize
                    % Seleção
                    % A função de seleção é chamada duas vezes para selecionar dois pais da população atual.
                    if strcmp(selection, 'Torneio')
                        obj.parent1 = obj.tournamentSelection(tournamentSize);
                        obj.parent2 = obj.tournamentSelection(tournamentSize);
                    end

                    if strcmp(selection, 'Roleta')
                        obj.parent1 = obj.rouletteWheelSelection();
                        obj.parent2 = obj.rouletteWheelSelection();
                    end

                    if strcmp(selection, 'Estoc')
                        parents = obj.stochasticUniversalSelection(estocSize);
                        obj.parent1 = parents(1,:);
                        obj.parent2 = parents(2,:);
                    end

                    % Crossover
                    % Aqui, um teste aleatório é feito contra a crossoverRate para decidir se o crossover ocorrerá. Se o crossover ocorrer,
                    % a função onePointCrossover combina os genes dos pais para criar um novo indivíduo (offspring).
                    % Se não, um dos pais é copiado diretamente para a nova população.
                    if rand <= obj.crossoverRate
                        if strcmp(crossover, '1P')
                            offspring = obj.onePointCrossover();
                        end
                        if strcmp(crossover, '2P')
                            offspring = obj.twoPointCrossover();
                        end
                        if strcmp(crossover, 'Uniform')
                            offspring = obj.uniformCrossover();
                        end
                        if strcmp(crossover, 'Blend')
                            offspring = obj.blendCrossover(alpha);
                        end

                    else
                        offspring = obj.parent1; % Sem crossover, copiar o pai diretamente
                    end

                    % Mutação
                    % A função mutate aplica mutações aleatórias ao novo indivíduo, com base na mutationRate. Cada gene tem
                    % uma chance de ser mutado dentro dos limites definidos por paramBounds.
                    if strcmp(mutation, 'Default')
                        offspring = obj.mutate(offspring);
                    end
                    if strcmp(mutation, 'Gauss')
                        offspring = obj.gaussianMutate(offspring);
                    end
                    if strcmp(mutation, 'Uniform')
                        offspring = obj.uniformMutate(offspring);
                    end
                    if strcmp(mutation, 'Creep')
                        offspring = obj.creepMutate(offspring);
                    end

                    % Adicionar à nova população
                    % O novo indivíduo é adicionado à nova população.
                    newPop(i, :) = offspring;
                end

                % Avaliar a nova população
                % Após a criação da nova população, cada indivíduo é avaliado pela função de custo f_motor,
                % e os escores de aptidão são atualizados.
                for i = 1:obj.popSize
                    obj.fitnessScores(i) = obj.CostFunction(newPop(i,:));
                end

                % Atualizar a população
                % A população atual é substituída pela nova população.
                obj.pop = newPop;

                % Atualizar métricas
                obj.bestFitnessHistory(gen) = min(obj.fitnessScores);
                obj.averageFitnessHistory(gen) = mean(obj.fitnessScores);

                % Encontrando o melhor indivíduo da geração atual
                [obj.currentBestFitness, bestIdx] = min(obj.fitnessScores);
                obj.bestIndividuals(gen, :) = newPop(bestIdx, :);

                disp(['Geração: ' num2str(gen) ' Melhor Valor: ' num2str(min(obj.fitnessScores)) ' Melhor Individuo: ' num2str(newPop(bestIdx, :))])
            end
            % Encontre o melhor indivíduo após a última geração
            [obj.bestFitness, bestIdx] = min(obj.fitnessScores);
            obj.bestSolution = obj.pop(bestIdx, :);

            % ao final da função, você insere o toc para obter a duração
            elapsedTime = toc;

            % Armazenando o tempo de execução no objeto
            obj.executionTime = elapsedTime;

            disp(['Tempo total de execução: ' num2str(floor(elapsedTime/60)) ' minutos e ' num2str(elapsedTime - floor(elapsedTime/60)*60) ' segundos.']);

        end

        %% Métodos de Seleção
        % Seleção por torneio
        function selected = tournamentSelection(obj, tournamentSize)
            % Objetivo: Selecionar um indivíduo da população (pop) para reprodução.
            % pop: População atual
            % fitnessScores: Vetor de aptidões dos indivíduos
            % tournamentSize: Número de indivíduos em cada torneio

            % Explicação do Método de Seleção por Torneio
            % Em um torneio, um pequeno grupo de indivíduos é escolhido aleatoriamente da população.
            % Dentro deste grupo, o indivíduo com a melhor aptidão é selecionado.
            % Este método permite que indivíduos mais aptos tenham uma maior chance de serem selecionados,
            % mas também oferece a indivíduos menos aptos uma chance de serem escolhidos, mantendo a diversidade genética na população.

            % Escolher aleatoriamente indivíduos para o torneio
            indices = randi(size(obj.pop, 1), tournamentSize, 1);
            tournament = obj.pop(indices, :);
            tournamentFitness = obj.fitnessScores(indices);

            % Encontrar o indivíduo com a melhor aptidão no torneio
            [~, bestIdx] = min(tournamentFitness);
            selected = tournament(bestIdx, :);

            % Vantagens da Seleção por Torneio
            % Simples e Eficiente: É fácil de implementar e computacionalmente eficiente.
            % Pressão Seletiva Ajustável: Alterando o tamanho do torneio (tournamentSize), você pode ajustar a pressão seletiva.
            % Um torneio maior aumenta a pressão seletiva (favorecendo indivíduos mais aptos), enquanto um torneio menor reduz essa pressão.
            % Diversidade Genética: Permite manter uma diversidade genética saudável na população, o que é crucial
            % para a eficácia do algoritmo genético ao longo de múltiplas gerações.
        end

        % Seleção por Roleta
        function selected = rouletteWheelSelection(obj)
            % Objetivo: selecionar indivíduos da população atual de uma maneira que a probabilidade de ser
            % escolhido seja proporcional à aptidão do indivíduo.
            % pop: População total
            % fitnessScores: Aptidão de toda a população

            % Invertendo as aptidões
            % A inversão é feita para que menores custos resultem em maiores valores de aptidão.
            invertedFitnessScores = 1 ./ obj.fitnessScores;

            % Tratando valores NaN
            % Substituindo NaNs por um valor pequeno para evitar problemas na seleção
            nanIndices = isnan(invertedFitnessScores);
            invertedFitnessScores(nanIndices) = 1e-10;

            % Calculando a aptidão total
            % A soma das aptidões invertidas de todos os indivíduos na população.
            totalFitness = sum(invertedFitnessScores);

            % Probabilidades de seleção
            % A probabilidade de cada indivíduo ser escolhido é agora proporcional à sua aptidão invertida.
            selectionProbabilities = invertedFitnessScores / totalFitness;

            % Probabilidades cumulativas
            % As probabilidades de seleção de cada indivíduo, acumuladas em ordem.
            cumulativeProbabilities = cumsum(selectionProbabilities);

            % Selecionando um indivíduo
            % Um número aleatório entre 0 e 1 é gerado. O primeiro indivíduo cuja probabilidade
            % cumulativa é maior ou igual a esse número aleatório é escolhido.
            r = rand; % Número aleatório entre 0 e 1
            selectedIdx = find(cumulativeProbabilities >= r, 1, 'first');
            selected = obj.pop(selectedIdx, :);
        end

        % Seleção Estocástica Universal
        function selected = stochasticUniversalSelection(obj, numSelected)
            % pop: População atual
            % fitnessScores: Vetor de aptidões dos indivíduos
            % numSelected: Número de indivíduos a serem selecionados

            % A seleção estocástica universal mantém um bom equilíbrio entre a pressão
            % seletiva e a preservação da diversidade genética.

            % Invertendo as aptidões para problemas de minimização
            % Se menor aptidão é melhor, inverter os valores
            invertedFitnessScores = 1 ./ obj.fitnessScores;
            totalFitness = sum(invertedFitnessScores);

            % Probabilidades de seleção
            selectionProbabilities = invertedFitnessScores / totalFitness;

            % Probabilidades cumulativas
            cumulativeProbabilities = cumsum(selectionProbabilities);

            % Inicializando o vetor de indivíduos selecionados
            selected = zeros(numSelected, size(obj.pop, 2));

            % Passo da roleta e posição inicial
            step = 1 / numSelected;
            start = rand() * step;

            % Selecionando os indivíduos
            for i = 1:numSelected
                pointer = start + (i - 1) * step;
                idx = find(cumulativeProbabilities >= pointer, 1, 'first');
                selected(i, :) = obj.pop(idx, :);
            end


        end

        %% Métodos de Crossover
        % Crossover de um ponto
        function offspring = onePointCrossover(obj)
            % Objetivo: Produzir um novo indivíduo (descendente) a partir de dois indivíduos existentes (pais).
            % parent1, parent2: Dois indivíduos pais

            % Número de genes (parâmetros)
            numGenes = length(obj.parent1);

            % Escolher um ponto de crossover aleatório
            % Um ponto ao longo do vetor de genes é escolhido aleatoriamente. Este ponto determina onde o crossover acontecerá.
            try
                crossoverPoint = randi(numGenes - 1);
            catch
                disp('erro');
            end

            % Criar descendentes
            % Dois potenciais descendentes são criados. Para cada um deles, uma parte dos genes é herdada de parent1 e a
            % outra parte de parent2.
            % offspring1 recebe os genes de parent1 até o crossoverPoint e então os genes de parent2 do crossoverPoint até o final
            offspring1 = [obj.parent1(1:crossoverPoint), obj.parent2(crossoverPoint + 1:end)];
            % offspring2 recebe os genes de parent2 até o crossoverPoint e então os genes de parent1 do crossoverPoint até o final.
            offspring2 = [obj.parent2(1:crossoverPoint), obj.parent1(crossoverPoint + 1:end)];

            % Escolher um dos descendentes aleatoriamente para retorno
            % Um dos dois descendentes criados é escolhido aleatoriamente para ser retornado pela função.
            % A chance de escolher qualquer um dos dois é a mesma (50%).
            if rand < 0.5
                offspring = offspring1;
            else
                offspring = offspring2;
            end

        end

        % Crossover de dois pontos
        function offspring = twoPointCrossover(obj)
            % Objetivo: Produzir um novo indivíduo (descendente) a partir de dois indivíduos existentes (pais)
            % usando crossover de dois pontos.
            % parent1, parent2: Dois indivíduos pais

            % Número de genes (parâmetros)
            numGenes = length(obj.parent1);

            % Escolher dois pontos de crossover aleatórios
            crossoverPoints = sort(randperm(numGenes - 1, 2));

            % Criar descendentes
            % Os segmentos entre os dois pontos são trocados entre os pais.
            offspring1 = [obj.parent1(1:crossoverPoints(1)), ...
                obj.parent2(crossoverPoints(1) + 1:crossoverPoints(2)), ...
                obj.parent1(crossoverPoints(2) + 1:end)];
            offspring2 = [obj.parent2(1:crossoverPoints(1)), ...
                obj.parent1(crossoverPoints(1) + 1:crossoverPoints(2)), ...
                obj.parent2(crossoverPoints(2) + 1:end)];

            % Escolher um dos descendentes aleatoriamente para retorno
            if rand < 0.5
                offspring = offspring1;
            else
                offspring = offspring2;
            end
            % O crossover de dois pontos tende a introduzir mais diversidade genética nos
            % descendentes do que o crossover de um ponto, pois permite a troca de segmentos maiores e mais variados dos cromossomos dos
            % pais.
            % Esta técnica pode ser particularmente útil em espaços de busca onde a correlação entre genes adjacentes é importante.
        end

        % Crossover Uniforme
        function offspring = uniformCrossover(obj)
            % Objetivo: Produzir um novo indivíduo (descendente) a partir de dois indivíduos existentes (pais) usando crossover uniforme.
            % parent1, parent2: Dois indivíduos pais

            % Número de genes (parâmetros)
            numGenes = length(obj.parent1);

            % Criar máscara de crossover
            % Uma máscara binária (vetor de valores lógicos) é gerada aleatoriamente,
            % onde true indica que o gene deve ser herdado do parent1 e false do parent2.
            crossoverMask = rand(numGenes, 1) < 0.5;

            % Criar descendente
            % O descendente herda genes de parent1 ou parent2 com base na máscara de crossover.
            offspring = zeros(1, numGenes);
            offspring(crossoverMask) = obj.parent1(crossoverMask);
            offspring(~crossoverMask) = obj.parent2(~crossoverMask);

            % O crossover uniforme é eficaz em manter a diversidade genética na população,
            % pois cada gene tem a mesma chance de ser herdado de qualquer um dos pais.
            % Este método é vantajoso quando a posição dos genes no cromossomo não é crucial para a aptidão do indivíduo.
        end

        % Crossover Blend (BLX - α)
        function offspring = blendCrossover(obj, alpha)
            % Objetivo: Produzir um novo indivíduo (descendente) a partir de dois indivíduos existentes (pais) usando o método BLX-.
            % parent1, parent2: Dois indivíduos pais
            % alpha: Parâmetro que controla o grau de exploração fora dos valores dos pais

            % Número de genes (parâmetros)
            numGenes = length(obj.parent1);

            % Criar descendente
            offspring = zeros(1, numGenes);
            for i = 1:numGenes
                % Calculando os limites para o crossover BLX-α
                minGene = min(obj.parent1(i), obj.parent2(i));
                maxGene = max(obj.parent1(i), obj.parent2(i));
                range = maxGene - minGene;

                % Gerando o gene do descendente
                % Este parâmetro determina o quão longe dos valores dos pais os
                % genes dos descendentes podem ser gerados.
                % Um valor comum é 0.5
                % Para cada gene, o descendente pode ter um valor dentro de um intervalo estendido
                % definido pelos genes dos pais mais uma margem dada por alpha.
                % offspring(i) = (minGene - alpha * range) + rand * (range + 2 * alpha * range);

                if(offspring(i) < 0)
                    disp('Deu ruim')
                end

                geneValue = (minGene - alpha * range) + rand * (range + 2 * alpha * range);
                geneValue = max(geneValue, obj.lowerBounds(1, i));  % Garante que não seja menor que o limite inferior
                geneValue = min(geneValue, obj.upperBounds(1, i));  % Garante que não seja maior que o limite superior

                offspring(i) = geneValue;
            end

            % O BLX-α é bom para explorar o espaço ao redor dos valores dos pais, o que pode ajudar a evitar mínimos locais.
        end

        %% Métodos de Mutação
        % Mutação Individual
        function mutatedIndividual = mutate(obj, individual)
            % Objetivo: Introduzir mutações aleatórias em um indivíduo para criar variação genética.
            % individual: Um indivíduo (solução) que será potencialmente mutado.
            % mutationRate: Probabilidade de mutação de cada gene
            % paramBounds: Limites dos parâmetros do motor

            % Número de genes (parâmetros)
            numGenes = length(individual);

            % Copiar o indivíduo para a mutação
            mutatedIndividual = individual;

            % Aplicar mutação a cada gene
            for i = 1:numGenes
                % Para cada gene no indivíduo:
                % Um número aleatório entre 0 e 1 é gerado. Se esse número é menor ou igual à mutationRate,
                % o gene é selecionado para mutação.
                if rand <= obj.mutationRate
                    % Garantir que estamos acessando a dimensão correta de paramBounds
                    lowerBound = obj.paramBounds(i,1);
                    upperBound = obj.paramBounds(i,2);

                    % Gerar um novo valor para o gene dentro dos limites
                    % O valor do gene mutado é calculado gerando um número aleatório dentro dos limites
                    % especificados (lowerBound e upperBound). Isso é feito para garantir que o valor mutado seja válido.
                    mutatedIndividual(i) = rand * (upperBound - lowerBound) + lowerBound;
                end
            end

            % Explicação da Mutação
            % Variação Genética: A mutação introduz novas variações genéticas na população, o que é vital para a
            % exploração eficaz do espaço de busca.
            % Taxa de Mutação: A mutationRate controla o equilíbrio entre manter a estabilidade genética da população
            % e introduzir novas variações. Uma taxa muito alta pode levar a uma exploração excessiva (potencialmente destrutiva),
            % enquanto uma taxa muito baixa pode resultar em convergência prematura.
            % Respeito aos Limites: A mutação é realizada dentro dos limites especificados para cada gene,
            % garantindo que os valores mutados sejam sempre viáveis no contexto do problema.
        end

        % Mutação Gaussiana
        function mutatedIndividual = gaussianMutate(obj, individual)
            % Objetivo: Introduzir mutações aleatórias em um indivíduo usando uma perturbação gaussiana.
            % individual: Um indivíduo (solução) que será potencialmente mutado.
            % mutationRate: Probabilidade de mutação de cada gene
            % paramBounds: Limites dos parâmetros do motor

            % Número de genes (parâmetros)
            numGenes = length(individual);

            % Copiar o indivíduo para a mutação
            mutatedIndividual = individual;

            % Aplicar mutação gaussiana a cada gene
            for i = 1:numGenes
                if rand <= obj.mutationRate
                    % Escolhendo um desvio padrão para a mutação
                    % O desvio padrão pode ser ajustado ou ser uma fração do intervalo do gene
                    % Este é um parâmetro crucial na mutação gaussiana. Ele determina o alcance da perturbação que
                    % será adicionada ao gene. O desvio padrão é definido como uma porcentagem do intervalo
                    % total do gene, mas você pode ajustá-lo conforme necessário.
                    stdDev = (obj.paramBounds(i,2) - obj.paramBounds(i,1)) * 0.1;  % Ajuste conforme necessário

                    % Aplicando a mutação gaussiana
                    mutatedIndividual(i) = mutatedIndividual(i) + randn * stdDev;

                    % Garantir que o valor mutado respeite os limites
                    mutatedIndividual(i) = max(min(mutatedIndividual(i), obj.paramBounds(i,2)), obj.paramBounds(i,1));
                end
            end
        end

        % Mutação Uniforme
        function mutatedIndividual = uniformMutate(obj, individual)
            % Objetivo: Introduzir mutações aleatórias em um indivíduo usando uma distribuição uniforme.
            % individual: Um indivíduo (solução) que será potencialmente mutado.
            % mutationRate: Probabilidade de mutação de cada gene
            % paramBounds: Limites dos parâmetros do motor

            % Número de genes (parâmetros)
            numGenes = length(individual);

            % Copiar o indivíduo para a mutação
            mutatedIndividual = individual;

            % Aplicar mutação uniforme a cada gene
            % Para cada gene selecionado para mutação, um novo valor é gerado aleatoriamente dentro do intervalo especificado para esse gene
            for i = 1:numGenes
                if rand <= obj.mutationRate
                    % Gerar um novo valor para o gene dentro dos limites
                    mutatedIndividual(i) = obj.paramBounds(i,1) + rand * (obj.paramBounds(i,2) - obj.paramBounds(i,1));
                end
            end
            % A mutação uniforme é eficaz para problemas onde é benéfico explorar o espaço
            % de busca de maneira ampla e introduzir variações significativas.
        end

        % Mutação Creep
        function mutatedIndividual = creepMutate(obj, individual)
            % Objetivo: Introduzir pequenas mutações em um indivíduo (solução).
            % individual: Um indivíduo que será potencialmente mutado.
            % mutationRate: Probabilidade de mutação de cada gene
            % paramBounds: Limites dos parâmetros do motor

            % Número de genes (parâmetros)
            numGenes = length(individual);

            % Copiar o indivíduo para a mutação
            mutatedIndividual = individual;

            % Aplicar mutação creep a cada gene
            for i = 1:numGenes
                if rand <= obj.mutationRate
                    % Definindo o tamanho do passo da mutação creep
                    % Esta variável define o tamanho máximo da alteração que será aplicada ao valor do gene.
                    % Normalmente, é uma pequena fração do intervalo total do gene.
                    creepStepSize = (obj.paramBounds(i,2) - obj.paramBounds(i,1)) * 0.02;  % Ajuste conforme necessário

                    % Aplicando a mutação creep
                    % Para cada gene selecionado para mutação, adiciona-se ou subtrai-se uma pequena quantidade do seu valor.
                    % O rand - 0.5 gera um número aleatório entre -0.5 e 0.5, e multiplicá-lo por 2 * creepStepSize garante
                    % que a mudança esteja dentro do intervalo desejado.
                    mutatedIndividual(i) = mutatedIndividual(i) + (rand - 0.5) * 2 * creepStepSize;

                    % Garantir que o valor mutado respeite os limites
                    mutatedIndividual(i) = max(min(mutatedIndividual(i), obj.paramBounds(i,2)), obj.paramBounds(i,1));
                end
            end
            % O tamanho do passo da mutação creep pode ser ajustado para controlar a magnitude da mutação. Um tamanho de passo
            % maior pode resultar em mudanças mais significativas nos valores dos genes.
        end

        %% Funções de Custo
        function [cost] = f_motor(obj, x)
            %% Parâmetros da Simulação

            f = 10000;                                     % Frequencia de amostragem do sinal
            Tsc = 1/f;                                     % Periodo de amostragem do sinal
            p = 10;                                        % Numero de partes que o intervalo discreto e dividido
            h = Tsc/p;                                     % Passo de amostragem continuo
            Tsimu = 0.7;                                    % Tempo de Simulação
            Np = Tsimu/Tsc;                             % Número de Pontos (vetores)

            %% Parâmetros do Motor

            Polos = 8;                                   % Numero de polos
            frequencia = 60;                             % Frequência elétrica nominal

            Lr = obj.Lm + obj.Lls;                               % Indutância dos enrolamentos do rotor

            J =  0.1633;                                 % Inércia do motor
            K = 0.12;                                    % Coeficiente de atrito do eixo do motor
            weles = 2*pi*frequencia;                     % Velocidade sincrona em radianos elétricos por segundo
            wr = 0;                                      % Velocidade inicial
            P = 4*736;                                   % Potência do motor

            %% Reatâncias para espaço de estados

            Xm = weles*obj.Lm;                               % Reatância de magnetização
            Xls = weles*obj.Lls;                             % Reatância de dispersão do estator
            Xlr = weles*obj.Lls;                             % Reatância de dispersão do rotor
            Xml = 1/(1/Xm + 1/Xls + 1/Xlr);

            %% Constantes para solução mecânica do motor

            A =  - K/J;
            B1 =   Polos/(2*J);
            B2 = - Polos/(2*J);

            %% Ganhos Controladores

            KP_w = x(1);
            KI_w = x(2);

            KP_id = x(3);
            KI_id = x(4);

            KP_iq = x(5);
            KI_iq = x(6);

            %% Inicialização das variaveis

            Fqs = 0;
            Fds = 0;
            Fqr = 0;
            Fdr = 0;
            Ids = 0;
            Ids_ant = 0;
            Iqs = 0;
            lambda_dr_est = 0;
            theta = 0;
            UI_w = 0;
            UI_id = 0;
            UI_iq = 0;
            t = 0:Tsc:Np*Tsc-Tsc;
            custo_erros = 0;
            Iqs_atrasado = 0;
            Ids_atrasado = 0;
            

            %% buffer de atraso
            % Inicializando buffers para armazenar os valores anteriores
            buffer_size = obj.buffer;  % Tamanho do atraso (número de amostras)
            buffer_index = 1;
            Iqs_buffer = zeros(1, buffer_size);
            Ids_buffer = zeros(1, buffer_size);

            %% Torque de Carga

            Tn = P*Polos/(2*weles);

            Tl =  Tn*0.75*100* ((t-0.6).*(t >= 0.6) - (t-0.61).*(t >= 0.61));

            %% Corrente Id de referência

            lambda_nonminal = 127/(2*pi*frequencia)/(obj.Lm);

            Ids_ref = 10*lambda_nonminal*((t-0).*(t >= 0) - (t-0.1).*(t >= 0.1));
            ids_vetor_auxiliar = zeros(1, length(Ids_ref));

            %% Velocidade de Referência

            w_nom_ref = 2*pi*60;

            w_ref = 2*w_nom_ref * ((t-0.1).*(t >= 0.1) - (t-0.6).*(t >= 0.6));

            %% Constantes
            cTl = h*B2*Tl;
            cA = h*A;
            cB1 = h*B1;
            cXls = 1/Xls;
            cTe = 3/2*Polos/2*1/weles;
            c1 = Xml/Xls;
            c2 = Xml/Xlr;
            c3 = (2*Lr-Tsc*obj.Rr)/(2*Lr+Tsc*obj.Rr);
            c4 = (obj.Lm*obj.Rr*Tsc)/(2*Lr+obj.Rr*Tsc);
            c5 = (obj.Lm*obj.Rr*Tsc)/(2*Lr+obj.Rr*Tsc);
            c6 = h*weles;
            c7 = c6*obj.Rs/Xls;
            c8 = c6*obj.Rr/Xlr;
            c9 = obj.Lm*obj.Rr/Lr;

            limit = 127 * sqrt(2);

            %% Loop Simulação Motor
            for k = 1:Np
                % Estimador de fluxo rotórico e orientação do sist. de referência

                lambda_dr_est = lambda_dr_est*c3 + Ids_atrasado*c4 + Ids_ant;
                Ids_ant = Ids_atrasado*c5;

                % Define um valor mínimo para lambda_dr_est para evitar divisão por um valor pequeno
                lambda_min = max(lambda_dr_est, 0.1);
                
                % Calcula wsl de forma eficiente
                wsl = c9 * Iqs_atrasado / lambda_min * (lambda_dr_est > 0.1);

                w = wr + wsl;
                theta = theta + Tsc*w;

                %Velocidade
                e_w = w_ref(k) - wr;
                UI_w = UI_w + e_w*Tsc;
                U_w = KP_w*e_w + KI_w*UI_w;
                iqs_ref = U_w;

                %Servos de corrente
                e_id = Ids_ref(k) - Ids_atrasado;
                UI_id = UI_id + e_id*Tsc;
                U_Id = KP_id*e_id + KI_id*UI_id;

                e_iq = iqs_ref - Iqs_atrasado;
                UI_iq = UI_iq*e_iq*Tsc;
                U_Iq = KP_iq*e_iq + KI_iq*UI_iq;

                %% Solucionando a EDO eletrica (euler)                
                U_Iq = min(max(U_Iq, -limit), limit);
                U_Id = min(max(U_Id, -limit), limit);
                Vq = U_Iq;
                Vd = U_Id;

                aux = h*w;

                for ksuper=1:p
                    Fqm = c1*Fqs + c2*Fqr;
                    Fdm = c1*Fds + c2*Fdr;

                    Fqs = Fqs + c6*Vq - aux*Fds - c7*(Fqs-Fqm);
                    Fds = Fds + c6*Vd + aux*Fqs - c7*(Fds-Fdm);

                    Fqr = Fqr - (aux*Fdr-h*wr*Fdr) + c8*(Fqm-Fqr);
                    Fdr = Fdr - (h*wr*Fqr-aux*Fqr) + c8*(Fdm-Fdr);

                    Iqs = (Fqs-Fqm)*cXls;
                    Ids = (Fds-Fdm)*cXls;

                    Te = cTe*(Fds*Iqs - Fqs*Ids);

                    % Solução mecânica

                    wr = wr + cA*wr + cB1*Te + cTl(k);

                end

                % Extração do valor atrasado do buffer usando o índice circular
                Iqs_atrasado = Iqs_buffer(buffer_index);
                Ids_atrasado = Ids_buffer(buffer_index);
            
                % Atualiza os buffers nos índices circulares com os novos valores
                Iqs_buffer(buffer_index) = Iqs;
                Ids_buffer(buffer_index) = Ids;
            
                % Atualiza o índice circular
                buffer_index = mod(buffer_index, buffer_size) + 1;

                % Aplicando os pesos no cálculo do erro
                custo_erros = custo_erros + Tsc*sqrt(e_w^2 + 2*e_id^2);

                % ids_vetor_auxiliar(k) = Ids_atrasado;
            end
            % ids_vetor_auxiliar(k+1) = Ids_atrasado;
            % % Cálculo do erro entre Ids e a corrente de referência
            % erro_ids = ids_vetor_auxiliar - Ids_ref;
            % 
            % % Penalização para erro integral absoluto
            % penalizacao_erro_integral = sum(abs(erro_ids))/1000;

            %% Penalização para zeros nos ganhos
            cost_KP_w = KP_w/10000;
            cost_KI_w = KI_w/1000;
            cost_KP_id = KP_id/1000;
            cost_KI_id = KI_id/1000;
            cost_KP_iq = KP_iq/1000;
            cost_KI_iq = KI_iq/1000;
            custo_ganhos = cost_KP_w + cost_KI_w + cost_KP_id + cost_KI_id + cost_KP_iq + cost_KI_iq;
            
            if (~KP_w && ~KI_w) || (~KP_id && ~KI_id) || (~KP_iq && ~KI_iq)
                cost = inf;
                return;
            end

            cost = 2*custo_erros + custo_ganhos;

        end

        %% Plotagens
        % Plota Motor
        function obj = plotMotor(obj)
            %% Parâmetros da Simulação

            f = 10000;                                     % Frequencia de amostragem do sinal
            Tsc = 1/f;                                     % Periodo de amostragem do sinal
            p = 10;                                        % Numero de partes que o intervalo discreto e dividido
            h = Tsc/p;                                     % Passo de amostragem continuo
            Tsimu = 2;                                    % Tempo de Simulação
            Np = Tsimu/Tsc;                                % Número de Pontos (vetores)

            %% Parâmetros do Motor

            Polos = 8;                                   % Numero de polos
            frequencia = 60;                             % Frequência elétrica nominal

            Lr = obj.Lm + obj.Lls;                               % Indutância dos enrolamentos do rotor

            J =  0.1633;                                 % Inércia do motor
            K = 0.12;                                    % Coeficiente de atrito do eixo do motor
            weles = 2*pi*frequencia;                     % Velocidade sincrona em radianos elétricos por segundo
            wr = 0;                                      % Velocidade inicial
            P = 4*736;                                   % Potência do motor

            %% Reatâncias para espaço de estados

            Xm = weles*obj.Lm;                               % Reatância de magnetização
            Xls = weles*obj.Lls;                             % Reatância de dispersão do estator
            Xlr = weles*obj.Lls;                             % Reatância de dispersão do rotor
            Xml = 1/(1/Xm + 1/Xls + 1/Xlr);

            %% Constantes para solução mecânica do motor

            A =  - K/J;
            B1 =   Polos/(2*J);
            B2 = - Polos/(2*J);

            %% Ganhos Controladores

            KP_w = obj.bestSolution(1);
            KI_w = obj.bestSolution(2);

            KP_id = obj.bestSolution(3);
            KI_id = obj.bestSolution(4);

            KP_iq = obj.bestSolution(5);
            KI_iq = obj.bestSolution(6);

            %% Inicialização das variaveis

            Fqs = 0;
            Fds = 0;
            Fqr = 0;
            Fdr = 0;
            Ids = 0;
            Ids_ant = 0;
            Iqs = 0;
            lambda_dr_est = 0;
            theta = 0;
            UI_w = 0;
            UI_id = 0;
            UI_iq = 0;
            t = 0:Tsc:Np*Tsc-Tsc;
            custo_erros = 0;
            Iqs_atrasado = 0;
            Ids_atrasado = 0;

            %% buffer de atraso
            % Inicializando buffers para armazenar os valores anteriores
            buffer_size = 10;  % Tamanho do atraso (número de amostras)
            buffer_index = 1;
            Iqs_buffer = zeros(1, buffer_size);
            Ids_buffer = zeros(1, buffer_size);

            %% Torque de Carga

            Tn = P*Polos/(2*weles);

            Tl =  Tn*0.75*20* ((t-1.3).*(t >= 1.3) - (t-1.35).*(t >= 1.35));

            %% Corrente Id de referência

            lambda_nonminal = 127/(2*pi*frequencia)/(obj.Lm);

            Ids_ref = 5*lambda_nonminal*((t-0).*(t >= 0) - (t-0.2).*(t >= 0.2));

            %% Velocidade de Referência

            w_nom_ref = 2*pi*60;

            w_ref = w_nom_ref * ((t-0.2).*(t >= 0.2) - (t-1.2).*(t >= 1.2));

            %% Constantes
            cTl = h*B2*Tl;
            cA = h*A;
            cB1 = h*B1;
            cXls = 1/Xls;
            cTe = 3/2*Polos/2*1/weles;
            c1 = Xml/Xls;
            c2 = Xml/Xlr;
            c3 = (2*Lr-Tsc*obj.Rr)/(2*Lr+Tsc*obj.Rr);
            c4 = (obj.Lm*obj.Rr*Tsc)/(2*Lr+obj.Rr*Tsc);
            c5 = (obj.Lm*obj.Rr*Tsc)/(2*Lr+obj.Rr*Tsc);
            c6 = h*weles;
            c7 = c6*obj.Rs/Xls;
            c8 = c6*obj.Rr/Xlr;
            c9 = obj.Lm*obj.Rr/Lr;

            % Definir o desvio padrão do ruído
            std_dev_noise = 1;  % Ajuste este valor conforme necessário

            %% Loop Simulação Motor
            for k = 1:Np
                % Estimador de fluxo rotórico e orientação do sist. de referência

                lambda_dr_est = lambda_dr_est*c3 + Ids_atrasado*c4 + Ids_ant*c5;
                Ids_ant = Ids_atrasado;

                % Define um valor mínimo para lambda_dr_est para evitar divisão por um valor pequeno
                lambda_min = max(lambda_dr_est, 0.1);
                
                % Calcula wsl de forma eficiente
                wsl = c9 * Iqs_atrasado / lambda_min * (lambda_dr_est > 0.1);

                wr_est = wr + wsl;
                w = wr_est;
                theta = theta + Tsc*w;

                %Velocidade
                e_w = w_ref(k) - wr;
                UI_w = UI_w + e_w*Tsc;
                U_w = KP_w*e_w + KI_w*UI_w;
                iqs_ref = U_w;

                %Servos de corrente
                e_id = Ids_ref(k) - Ids_atrasado;
                UI_id = UI_id + e_id*Tsc;
                U_Id = KP_id*e_id + KI_id*UI_id;

                e_iq = iqs_ref - Iqs_atrasado;
                UI_iq = UI_iq*e_iq*Tsc;
                U_Iq = KP_iq*e_iq + KI_iq*UI_iq;

                %% Solucionando a EDO eletrica (euler)
                limit = 127 * sqrt(2);
                
                U_Iq = min(max(U_Iq, -limit), limit);
                U_Id = min(max(U_Id, -limit), limit);
                Vq = U_Iq;
                Vd = U_Id;

                for ksuper=1:p
                    Fqm = c1*Fqs + c2*Fqr;
                    Fdm = c1*Fds + c2*Fdr;

                    Fqs = Fqs + c6*Vq - h*w*Fds - c7*(Fqs-Fqm);
                    Fds = Fds + c6*Vd + h*w*Fqs - c7*(Fds-Fdm);

                    Fqr = Fqr - h*(w-wr)*Fdr + c8*(Fqm-Fqr);
                    Fdr = Fdr - h*(wr-w)*Fqr + c8*(Fdm-Fdr);

                    Iqs = (Fqs-Fqm)*cXls;
                    Ids = (Fds-Fdm)*cXls;

                    % Adicionar ruído às medições de corrente
                    Iqs_measured = Iqs + std_dev_noise * randn();
                    Ids_measured = Ids + std_dev_noise * randn();

                    Te = cTe*(Fds*Iqs - Fqs*Ids);

                    % Solução mecânica

                    wr = wr + cA*wr + cB1*Te + cTl(k);

                end

                % Extração do valor atrasado do buffer usando o índice circular
                Iqs_atrasado = Iqs_buffer(buffer_index);
                Ids_atrasado = Ids_buffer(buffer_index);
            
                % Atualiza os buffers nos índices circulares com os novos valores
                Iqs_buffer(buffer_index) = Iqs_measured;
                Ids_buffer(buffer_index) = Ids_measured;
            
                % Atualiza o índice circular
                buffer_index = mod(buffer_index, buffer_size) + 1;

                obj.iqs_vetor(k) = Iqs;
                iqs_atrasado_vetor(k) = Iqs_atrasado;
                ids_atrasado_vetor(k) = Ids_atrasado;
                obj.ids_vetor(k) = Ids;
                obj.te_vetor(k) = Te;
                obj.wr_vetor(k) = wr;
            end

            %% graficos
            figure
            % Plotando os dados
            plot(t,obj.wr_vetor, t, w_ref); % tom de cinza escuro
            legend('Wr','Ref')

            figure
            % Plotando os dados
            plot(t,obj.ids_vetor,t,ids_atrasado_vetor, t, Ids_ref); % tom de cinza escuro
            legend('Ids', 'Ids atrasado', 'Referencia')

            figure
            % Plotando os dados
            plot(t,obj.iqs_vetor, t ,iqs_ref, t, iqs_atrasado_vetor); % tom de cinza escuro
            legend('Iqs', 'Ref', 'Iqs Atrasado')

            figure
            % Plotando os dados
            plot(t,obj.te_vetor,t, Tl); % tom de cinza escuro
            legend('Te', 'TL')
        end

        function plotGA(obj)
            % Plotar a melhor aptidão ao longo das gerações
            figure;
            plot(1:obj.numGenerations, obj.bestFitnessHistory, 'b-', 'LineWidth', 2);
            title('Melhor Aptidão por Geração');
            xlabel('Geração');
            ylabel('Melhor Aptidão');
            grid on;

            % Plotar a aptidão média ao longo das gerações
            figure;
            plot(1:obj.numGenerations, obj.averageFitnessHistory, 'r-', 'LineWidth', 2);
            title('Aptidão Média por Geração');
            xlabel('Geração');
            ylabel('Aptidão Média');
            grid on;

            % Plotando a evolução de cada gene
            figure;
            for i = 1:4
                subplot(4, 1, i);
                plot(1:obj.numGenerations, obj.bestIndividuals(:, i), 'LineWidth', 2);
                title(['Evolução do Gene ', num2str(i)]);
                xlabel('Geração');
                ylabel(['Gene ', num2str(i)]);
                grid on;
            end
        end
    end
end