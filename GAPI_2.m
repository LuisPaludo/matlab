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
        wrpm_vetor
        iqs_vetor
        ids_vetor
        t_vetor
        executionTime
        simp_wr
        simp_iqs
        simp_ids
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
            obj.simp_wr = 891.11;
            obj.simp_iqs = 4.85;
            obj.simp_ids = 9.72;

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
                     = selectionArgs{2};
            end
            if strcmp(selection, 'Boltzmann')
                % disp('Digite a temperatura Inicial')
                tempInicial = selectionArgs{2}(1);
                % disp('Digite a temperatura Final')
                tempFinal = selectionArgs{2}(2);
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

                    if strcmp(selection, 'Rank')
                        obj.parent1 = obj.rankSelection();
                        obj.parent2 = obj.rankSelection();
                    end

                    if strcmp(selection, 'Estoc')
                        parents = obj.stochasticUniversalSelection(estocSize);
                        obj.parent1 = parents(1,:);
                        obj.parent2 = parents(2,:);
                    end

                    if strcmp(selection, 'Boltzmann')
                        obj.parent1 = obj.boltzmannSelection(tempInicial, tempFinal);
                        obj.parent2 = obj.boltzmannSelection(tempInicial, tempFinal);
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
                        if strcmp(crossover, 'Arit')
                            offspring = obj.arithmeticCrossover();
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
                    if strcmp(mutation, 'Inv')
                        offspring = obj.inversionMutate(offspring);
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

        % Seleção por Rank
        function selected = rankSelection(obj)
            % Objetivo: Os indivíduos são classificados com base em sua aptidão,
            % e a probabilidade de seleção é determinada pela sua classificação, não pela sua aptidão absoluta.
            % pop: População atual
            % fitnessScores: Vetor de aptidões dos indivíduos

            % Invertendo as aptidões
            % A inversão é feita para que menores custos resultem em maiores valores de aptidão.
            invertedFitnessScores = 1 ./ obj.fitnessScores;
            % Tratando valores NaN
            % Substituindo NaNs por um valor pequeno para evitar problemas na seleção
            nanIndices = isnan(invertedFitnessScores);
            invertedFitnessScores(nanIndices) = 1e-10;
            % Classificar os indivíduos
            [~, sortedIndices] = sort(invertedFitnessScores);

            % Atribuir probabilidades baseadas na classificação
            numIndividuals = length(invertedFitnessScores);
            rankProbabilities = ((1:numIndividuals) / sum(1:numIndividuals))';

            % Probabilidades cumulativas
            cumulativeProbabilities = cumsum(rankProbabilities);

            % Selecionando um indivíduo
            r = rand;
            search = find(cumulativeProbabilities >= r, 1, 'first');
            selectedIdx = sortedIndices(search);
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

        % Seleção de Boltzmann
        function selected = boltzmannSelection(obj, tempInicial, tempFinal)
            % Objetivo: Alterar a escala de aptidão dos indivíduos com base em uma "temperatura" que diminui ao longo do tempo
            % pop: População atual
            % fitnessScores: Vetor de aptidões dos indivíduos
            % gen: Geração atual
            % maxGen: Número máximo de gerações

            % Parâmetros da seleção de Boltzmann
            % Os valores de T_initial e T_final podem ser ajustados de acordo com as necessidades específicas do seu problema.
            % A "temperatura" começa alta e diminui a cada geração, alterando a escala de aptidão. Temperaturas mais altas
            % favorecem a exploração (menos pressão seletiva), e temperaturas mais baixas favorecem a explotação (maior pressão seletiva).
            T_initial = tempInicial; % Temperatura inicial
            T_final = tempFinal; % Temperatura final
            % Calculando a temperatura atual com base na geração
            T = T_initial * (T_final / T_initial) ^ (obj.generation / obj.numGenerations);

            % Ajustando as aptidões com base na distribuição de Boltzmann
            % As aptidões são ajustadas usando a função exponencial negativa baseada na temperatura atual.
            % Isso altera as probabilidades de seleção dos indivíduos.
            adjustedFitnessScores = exp(-(obj.fitnessScores/10) / T);
            totalAdjustedFitness = sum(adjustedFitnessScores);

            % Probabilidades de seleção
            selectionProbabilities = adjustedFitnessScores / totalAdjustedFitness;

            % Probabilidades cumulativas
            cumulativeProbabilities = cumsum(selectionProbabilities);

            % Selecionando um indivíduo
            r = rand;
            selectedIdx = find(cumulativeProbabilities >= r, 1, 'first');
            selected = obj.pop(selectedIdx, :);
            if(isempty(selected))
                disp('vai dar erro')
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

        % Crossover Aritmético
        function offspring = arithmeticCrossover(obj)
            % Objetivo: Produzir um novo indivíduo (descendente) a partir de dois
            % indivíduos existentes (pais) usando crossover aritmético.
            % parent1, parent2: Dois indivíduos pais

            % Número de genes (parâmetros)
            numGenes = length(obj.parent1);

            % Inicializar descendente com zeros
            offspring = zeros(1, numGenes);

            % Crossover aritmético para cada gene
            for i = 1:numGenes
                % Coeficiente aleatório para cada gene
                alpha = rand();

                % Valor do gene do descendente como média ponderada
                offspring(i) = alpha * obj.parent1(i) + (1 - alpha) * obj.parent2(i);
            end

            % Este método é eficaz para problemas com variáveis contínuas.
            % Ele pode criar descendentes dentro do intervalo definido pelos pais,
            % mas também pode limitar a diversidade se os pais forem muito semelhantes.
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

        % Mutação por Inversão
        function mutatedIndividual = inversionMutate(obj, individual)
            % Objetivo: Aplicar a mutação por inversão em um indivíduo.
            % individual: Um indivíduo (solução) que será potencialmente mutado.
            % mutationRate: Probabilidade de mutação de cada gene

            % Número de genes (parâmetros)
            numGenes = length(individual);

            % Copiar o indivíduo para a mutação
            mutatedIndividual = individual;

            % Verificar se a mutação ocorrerá
            if rand <= obj.mutationRate
                % Escolhendo dois pontos aleatórios para a inversão
                % Dois pontos são selecionados aleatoriamente dentro do cromossomo. Estes pontos definem o segmento que será invertido.
                points = randperm(numGenes, 2);
                lowerPoint = min(points);
                upperPoint = max(points);

                % Invertendo o segmento entre os dois pontos
                % Os genes entre os dois pontos são invertidos em ordem. Isso pode ser útil em problemas onde a ordem dos genes é importante.
                mutatedIndividual(lowerPoint:upperPoint) = mutatedIndividual(upperPoint:-1:lowerPoint);
            end
        end

        %% Funções de Custo

        function [cost] = f_motor(obj, x)
            %% Parâmetros da Simulação

            f = 10000;                                     % Frequencia de amostragem do sinal
            Tsc = 1/f;                                     % Periodo de amostragem do sinal
            p = 10;                                        % Numero de partes que o intervalo discreto e dividido
            h = Tsc/p;                                     % Passo de amostragem continuo
            Tsimu = 1;                                    % Tempo de Simulação
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
            Iqs_buffer = zeros(1, buffer_size);
            Ids_buffer = zeros(1, buffer_size);

            %% Torque de Carga

            Tn = P*Polos/(2*weles);

            Tl =  Tn*0.75*10* ((t-0.7).*(t >= 0.7) - (t-0.8).*(t >= 0.8));

            %% Corrente Id de referência

            lambda_nonminal = 127/(2*pi*frequencia)/(obj.Lm);

            Ids_ref = 5*lambda_nonminal*((t-0).*(t >= 0) - (t-0.2).*(t >= 0.2));

            %% Velocidade de Referência

            w_nom_ref = 2*pi*60;

            w_ref = 2*w_nom_ref * ((t-0.2).*(t >= 0.2) - (t-0.7).*(t >= 0.7));

            %% Loop Simulação Motor
            for k = 1:Np
                % Estimador de fluxo rotórico e orientação do sist. de referência

                lambda_dr_est = lambda_dr_est*((2*Lr-Tsc*obj.Rr)/(2*Lr+Tsc*obj.Rr)) + Ids_atrasado*((obj.Lm*obj.Rr*Tsc)/(2*Lr+obj.Rr*Tsc)) + Ids_ant*((obj.Lm*obj.Rr*Tsc)/(2*Lr+obj.Rr*Tsc));
                Ids_ant = Ids_atrasado;

                if(lambda_dr_est > 0.1)
                    wsl = (obj.Lm*obj.Rr*Iqs_atrasado)/(Lr*lambda_dr_est);
                else
                    wsl = 0;
                end

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

                if(U_Id >= 127*sqrt(2))
                    U_Id = 127*sqrt(2);
                end
                if(U_Id <= -127*sqrt(2))
                    U_Id = -127*sqrt(2);
                end

                e_iq = iqs_ref - Iqs_atrasado;
                UI_iq = UI_iq*e_iq*Tsc;
                U_Iq = KP_iq*e_iq + KI_iq*UI_iq;

                if(U_Iq >= 127*sqrt(2))
                    U_Iq = 127*sqrt(2);
                end
                if(U_Iq <= -127*sqrt(2))
                    U_Iq = -127*sqrt(2);
                end

                %% Solucionando a EDO eletrica (euler)
                Vq = U_Iq;
                Vd = U_Id;

                %% Calculando as Tensões Va Vb e Vc

                Valfa = Vq*cos(theta) + Vd*sin(theta);
                Vbeta = -Vq*sin(theta) + Vd*cos(theta);

                %Transf. inversa clarke
                Va = Valfa;
                Vb = -0.5*Valfa - sqrt(3)/2*Vbeta;
                Vc = -0.5*Valfa + sqrt(3)/2*Vbeta;

                Vmax = 127*sqrt(2);

                if abs(Va) > Vmax || abs(Vb) > Vmax || abs(Vc) > Vmax
                    % Calcule o fator de redução baseado na maior tensão
                    scalingFactor = Vmax / max(abs([Va, Vb, Vc]));

                    Vd = Vd * scalingFactor;
                    Vq = Vq * scalingFactor;
                end


                for ksuper=1:p

                    Fqm = Xml/Xls*Fqs + Xml/Xlr*Fqr;
                    Fdm = Xml/Xls*Fds + Xml/Xlr*Fdr;

                    Fqs = Fqs + h*weles*(Vq - w/weles*Fds - obj.Rs/Xls*(Fqs-Fqm));
                    Fds = Fds + h*weles*(Vd + w/weles*Fqs - obj.Rs/Xls*(Fds-Fdm));

                    Fqr = Fqr - h*weles*((w-wr)*Fdr/weles + obj.Rr/Xlr*(Fqr-Fqm));
                    Fdr = Fdr - h*weles*((wr-w)*Fqr/weles + obj.Rr/Xlr*(Fdr-Fdm));

                    Iqs = (Fqs-Fqm)/Xls;
                    Ids = (Fds-Fdm)/Xls;

                    Te = 3/2*Polos/2*1/weles*(Fds*Iqs - Fqs*Ids);

                    % Solução mecânica

                    wr = wr + h*(A*wr + B1*Te + B2*Tl(k));

                end

                % Desloca o buffer para a direita (atualiza o atraso)
                Iqs_buffer = [Iqs, Iqs_buffer(1:end-1)];
                Ids_buffer = [Ids, Ids_buffer(1:end-1)];
                
                % Extrai os valores atrasados do buffer
                Iqs_atrasado = Iqs_buffer(end);  % Valor com atraso de buffer_size amostras
                Ids_atrasado = Ids_buffer(end);  % Valor com atraso de buffer_size amostras

                %% Penalização para zeros nos ganhos
                cost_KP_w = KP_w/10000;
                cost_KI_w = KI_w/1000;
                cost_KP_id = KP_id/1000;
                cost_KI_id = KI_id/1000;
                cost_KP_iq = KP_iq/1000;
                cost_KI_iq = KI_iq/1000;
                custo_ganhos = cost_KP_w + cost_KI_w + cost_KP_id + cost_KI_id + cost_KP_iq + cost_KI_iq;

                % Aplicando os pesos no cálculo do erro
                custo_erros = custo_erros + Tsc*sqrt(e_w^2 + e_id^2);

                cost = 2*custo_erros + custo_ganhos;

            end
        end

        %% Plotagens
        % Plota Motor
        function obj = plotMotor(obj)
            %% Simulação do Motor com os parâmetros obtidos

            % Parâmetros de Simulação
            f = 10000;                                     % Frequencia de amostragem do sinal
            Tsc = 1/f;                                     % Periodo de amostragem do sinal
            p = 10;                                        % Numero de partes que o intervalo discreto e dividido
            h = Tsc/p;                                     % Passo de amostragem continuo
            Tsimu = 1;
            Np = Tsimu/Tsc;

            % Parâmetros do motor
            Polos = 8;                                   % Numero de polos
            frequencia = 60;                             % Frequência elétrica nominal
            Rs = obj.bestSolution(1);                        % Resistência Estatórica
            Lls = obj.bestSolution(4);
            Lm = obj.bestSolution(3);
            Llr = obj.bestSolution(4);
            Rr = obj.bestSolution(2);
            J =  0.1633;                                 % Inércia do motor
            K = 0.12;                                    % Coeficiente de atrito do eixo do motor
            weles = 2*pi*frequencia;                     % Velocidade sincrona em radianos elétricos por segundo
            w = weles;                                   % Velocidade do sistema síncrono
            Tl = 0;                                      % Torque de carga
            wr = 0;

            Vq = 150;
            Vd = 5.7;

            % Reatâncias para espaço de estados

            Xm = weles*Lm;                     % Reatância de magnetização
            Xls = weles*Lls;                   % Reatância de dispersão do estator
            Xlr = weles*Llr;                   % Reatância de dispersão do rotor
            Xml = 1/(1/Xm + 1/Xls + 1/Xlr);

            % Constantes para solução mecânica do motor

            A =  - K/J;
            B1 =   Polos/(2*J);
            B2 = - Polos/(2*J);

            % Inicialização das variaveis
            Fqs = 0;
            Fds = 0;
            Fqr = 0;
            Fdr = 0;

            % Inicialização de vetores
            obj.wrpm_vetor = zeros(1, Np);
            obj.iqs_vetor = zeros(1, Np);
            obj.ids_vetor = zeros(1, Np);
            obj.t_vetor = linspace(0, Tsimu, Np); % Cria um vetor de tempo linearmente espaçado

            % Loop Simulação Motor
            for k = 1:Np

                for ksuper=1:p

                    Fqm = Xml/Xls*Fqs + Xml/Xlr*Fqr;
                    Fdm = Xml/Xls*Fds + Xml/Xlr*Fdr;

                    Fqs = Fqs + h*weles*(Vq - w/weles*Fds - Rs/Xls*(Fqs-Fqm));
                    Fds = Fds + h*weles*(Vd + w/weles*Fqs - Rs/Xls*(Fds-Fdm));

                    Fqr = Fqr - h*weles*((w-wr)*Fdr/weles + Rr/Xlr*(Fqr-Fqm));
                    Fdr = Fdr - h*weles*((wr-w)*Fqr/weles + Rr/Xlr*(Fdr-Fdm));

                    Iqs = (Fqs-Fqm)/Xls;
                    Ids = (Fds-Fdm)/Xls;

                    Te = 3/2*Polos/2*1/weles*(Fds*Iqs - Fqs*Ids);

                    % Solução mecânica

                    wr = wr + h*(A*wr + B1*Te + B2*Tl);
                    wrpm = wr*2/Polos*60/(2*pi);

                end

                obj.wrpm_vetor(k) = wrpm;
                obj.iqs_vetor(k) = Iqs;
                obj.ids_vetor(k) = Ids;
                obj.t_vetor(k) = Tsc*(k-1);
            end


            figure; % Cria uma nova janela de figura
            plot(obj.t_vetor, obj.wrpm_vetor,obj.t_vetor, obj.real_wr_vetor, 'LineWidth', 2); % Plota wrpm_vetor
            title('Velocidade (wrpm) x Tempo');
            legend('Estimada','Real')
            xlabel('Tempo (s)');
            ylabel('Velocidade (wrpm)');
            grid on; % Adiciona uma grade para melhor visualização

            figure; % Cria outra nova janela de figura
            plot(obj.t_vetor, obj.iqs_vetor,obj.t_vetor, obj.real_iqs_vetor, 'LineWidth', 2); % Plota iqs_vetor
            title('Corrente iqs x Tempo');
            legend('Estimada','Real')
            xlabel('Tempo (s)');
            ylabel('Corrente iqs (A)');
            grid on;

            figure; % Cria mais uma nova janela de figura
            plot(obj.t_vetor, obj.ids_vetor,obj.t_vetor, obj.real_ids_vetor, 'LineWidth', 2); % Plota ids_vetor
            title('Corrente ids x Tempo');
            legend('Estimada','Real')
            xlabel('Tempo (s)');
            ylabel('Corrente ids (A)');
            grid on;
        end

        function obj = plotMotorComplement(obj)
            %% Simulação do Motor com os parâmetros obtidos

            % Parâmetros de Simulação
            f = 10000;                                     % Frequencia de amostragem do sinal
            Tsc = 1/f;                                     % Periodo de amostragem do sinal
            p = 10;                                        % Numero de partes que o intervalo discreto e dividido
            h = Tsc/p;                                     % Passo de amostragem continuo
            Tsimu = 1;
            Np = Tsimu/Tsc;

            % Parâmetros do motor
            Polos = 8;                                   % Numero de polos
            frequencia = 60;                             % Frequência elétrica nominal
            Rs = obj.Rs_solution;                        % Resistência Estatórica
            Lls = obj.bestSolution(1);
            Lm = obj.Lm_solution;
            Llr = obj.bestSolution(1);
            Rr = obj.Rr_solution;
            J =  0.1633;                                 % Inércia do motor
            K = 0.12;                                    % Coeficiente de atrito do eixo do motor
            weles = 2*pi*frequencia;                     % Velocidade sincrona em radianos elétricos por segundo
            w = weles;                                   % Velocidade do sistema síncrono
            Tl = 0;                                      % Torque de carga
            wr = 0;

            Vq = 150;
            Vd = 5.7;

            % Reatâncias para espaço de estados

            Xm = weles*Lm;                     % Reatância de magnetização
            Xls = weles*Lls;                   % Reatância de dispersão do estator
            Xlr = weles*Llr;                   % Reatância de dispersão do rotor
            Xml = 1/(1/Xm + 1/Xls + 1/Xlr);

            % Constantes para solução mecânica do motor

            A =  - K/J;
            B1 =   Polos/(2*J);
            B2 = - Polos/(2*J);

            % Inicialização das variaveis
            Fqs = 0;
            Fds = 0;
            Fqr = 0;
            Fdr = 0;

            % Inicialização de vetores
            obj.wrpm_vetor = zeros(1, Np);
            obj.iqs_vetor = zeros(1, Np);
            obj.ids_vetor = zeros(1, Np);
            obj.t_vetor = linspace(0, Tsimu, Np); % Cria um vetor de tempo linearmente espaçado

            % Loop Simulação Motor
            for k = 1:Np

                for ksuper=1:p

                    Fqm = Xml/Xls*Fqs + Xml/Xlr*Fqr;
                    Fdm = Xml/Xls*Fds + Xml/Xlr*Fdr;

                    Fqs = Fqs + h*weles*(Vq - w/weles*Fds - Rs/Xls*(Fqs-Fqm));
                    Fds = Fds + h*weles*(Vd + w/weles*Fqs - Rs/Xls*(Fds-Fdm));

                    Fqr = Fqr - h*weles*((w-wr)*Fdr/weles + Rr/Xlr*(Fqr-Fqm));
                    Fdr = Fdr - h*weles*((wr-w)*Fqr/weles + Rr/Xlr*(Fdr-Fdm));

                    Iqs = (Fqs-Fqm)/Xls;
                    Ids = (Fds-Fdm)/Xls;

                    Te = 3/2*Polos/2*1/weles*(Fds*Iqs - Fqs*Ids);

                    % Solução mecânica

                    wr = wr + h*(A*wr + B1*Te + B2*Tl);
                    wrpm = wr*2/Polos*60/(2*pi);

                end

                obj.wrpm_vetor(k) = wrpm;
                obj.iqs_vetor(k) = Iqs;
                obj.ids_vetor(k) = Ids;
                obj.t_vetor(k) = Tsc*(k-1);
            end


            figure; % Cria uma nova janela de figura
            plot(obj.t_vetor, obj.wrpm_vetor,obj.t_vetor, obj.real_wr_vetor, 'LineWidth', 2); % Plota wrpm_vetor
            title('Velocidade (wrpm) x Tempo');
            legend('Estimada','Real')
            xlabel('Tempo (s)');
            ylabel('Velocidade (wrpm)');
            grid on; % Adiciona uma grade para melhor visualização

            figure; % Cria outra nova janela de figura
            plot(obj.t_vetor, obj.iqs_vetor,obj.t_vetor, obj.real_iqs_vetor, 'LineWidth', 2); % Plota iqs_vetor
            title('Corrente iqs x Tempo');
            legend('Estimada','Real')
            xlabel('Tempo (s)');
            ylabel('Corrente iqs (A)');
            grid on;

            figure; % Cria mais uma nova janela de figura
            plot(obj.t_vetor, obj.ids_vetor,obj.t_vetor, obj.real_ids_vetor, 'LineWidth', 2); % Plota ids_vetor
            title('Corrente ids x Tempo');
            legend('Estimada','Real')
            xlabel('Tempo (s)');
            ylabel('Corrente ids (A)');
            grid on;
        end

        % Plota Motor
        function obj = plotMotorSimplified(obj)
            %% Simulação do Motor com os parâmetros obtidos

            % Parâmetros de Simulação
            f = 10000;                                     % Frequencia de amostragem do sinal
            Tsc = 1/f;                                     % Periodo de amostragem do sinal
            p = 10;                                        % Numero de partes que o intervalo discreto e dividido
            h = Tsc/p;                                     % Passo de amostragem continuo
            Tsimu = 1;
            Np = Tsimu/Tsc;

            % Parâmetros do motor
            Polos = 8;                                   % Numero de polos
            frequencia = 60;                             % Frequência elétrica nominal
            Rs = obj.bestSolution(1);                        % Resistência Estatórica
            Lls = obj.bestSolution(4);
            Lm = obj.bestSolution(3);
            Llr = obj.bestSolution(4);
            Rr = obj.bestSolution(2);
            J =  0.1633;                                 % Inércia do motor
            K = 0.12;                                    % Coeficiente de atrito do eixo do motor
            weles = 2*pi*frequencia;                     % Velocidade sincrona em radianos elétricos por segundo
            w = weles;                                   % Velocidade do sistema síncrono
            Tl = 0;                                      % Torque de carga
            wr = 2*pi*frequencia;

            Vq = 150;
            Vd = 5.7;

            % Reatâncias para espaço de estados

            Xm = weles*Lm;                     % Reatância de magnetização
            Xls = weles*Lls;                   % Reatância de dispersão do estator
            Xlr = weles*Llr;                   % Reatância de dispersão do rotor
            Xml = 1/(1/Xm + 1/Xls + 1/Xlr);

            % Constantes para solução mecânica do motor

            A =  - K/J;
            B1 =   Polos/(2*J);
            B2 = - Polos/(2*J);

            % Inicialização das variaveis
            Fqs = 0;
            Fds = 0;
            Fqr = 0;
            Fdr = 0;

            % Inicialização de vetores
            obj.wrpm_vetor = zeros(1, Np);
            obj.iqs_vetor = zeros(1, Np);
            obj.ids_vetor = zeros(1, Np);
            obj.t_vetor = linspace(0, Tsimu, Np); % Cria um vetor de tempo linearmente espaçado

            % Loop Simulação Motor
            for k = 1:Np

                for ksuper=1:p

                    Fqm = Xml/Xls*Fqs + Xml/Xlr*Fqr;
                    Fdm = Xml/Xls*Fds + Xml/Xlr*Fdr;

                    Fqs = Fqs + h*weles*(Vq - w/weles*Fds - Rs/Xls*(Fqs-Fqm));
                    Fds = Fds + h*weles*(Vd + w/weles*Fqs - Rs/Xls*(Fds-Fdm));

                    Fqr = Fqr - h*weles*((w-wr)*Fdr/weles + Rr/Xlr*(Fqr-Fqm));
                    Fdr = Fdr - h*weles*((wr-w)*Fqr/weles + Rr/Xlr*(Fdr-Fdm));

                    Iqs = (Fqs-Fqm)/Xls;
                    Ids = (Fds-Fdm)/Xls;

                    Te = 3/2*Polos/2*1/weles*(Fds*Iqs - Fqs*Ids);

                    % Solução mecânica

                    wr = wr + h*(A*wr + B1*Te + B2*Tl);
                    wrpm = wr*2/Polos*60/(2*pi);

                end

                obj.wrpm_vetor(k) = wrpm;
                obj.iqs_vetor(k) = Iqs;
                obj.ids_vetor(k) = Ids;
                obj.t_vetor(k) = Tsc*(k-1);
            end
            %
            % load('Dados/Iqs.mat');
            % load('Dados/Ids.mat');
            % load('Dados/wr.mat');

            constante_wr_vetor = 891.11 * ones(size(obj.t_vetor));
            constante_iqs_vetor = 4.85 * ones(size(obj.t_vetor));
            constante_ids_vetor = 9.72 * ones(size(obj.t_vetor));

            figure; % Cria uma nova janela de figura
            plot(obj.t_vetor, obj.wrpm_vetor,obj.t_vetor, constante_wr_vetor, 'LineWidth', 2); % Plota wrpm_vetor
            title('Velocidade (wrpm) x Tempo');
            legend('Estimada','Real')
            xlabel('Tempo (s)');
            ylabel('Velocidade (wrpm)');
            grid on; % Adiciona uma grade para melhor visualização
            xlim([0.8, 1])
            ylim([0.98*891.11, 1.02*891.11])

            figure; % Cria outra nova janela de figura
            plot(obj.t_vetor, obj.iqs_vetor,obj.t_vetor, constante_iqs_vetor, 'LineWidth', 2); % Plota iqs_vetor
            title('Corrente iqs x Tempo');
            legend('Estimada','Real')
            xlabel('Tempo (s)');
            ylabel('Corrente iqs (A)');
            grid on;
            xlim([0.8, 1])
            % ylim([0.98*4.85, 1.02*4.85])

            figure; % Cria mais uma nova janela de figura
            plot(obj.t_vetor, obj.ids_vetor,obj.t_vetor, constante_ids_vetor, 'LineWidth', 2); % Plota ids_vetor
            title('Corrente ids x Tempo');
            legend('Estimada','Real')
            xlabel('Tempo (s)');
            ylabel('Corrente ids (A)');
            grid on;
            xlim([0.8, 1])
            ylim([0.98*9.72, 1.02*9.72])
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