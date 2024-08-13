classdef GAPI
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
        real_iqs_vetor
        real_ids_vetor
        real_wr_vetor
        simp_wr
        simp_iqs
        simp_ids
        Rs {mustBePositive, mustBeFloat}
        Rr {mustBePositive, mustBeFloat}
        Lm {mustBePositive, mustBeFloat}
        Lls {mustBePositive, mustBeFloat}
        Llr {mustBePositive, mustBeFloat}
        cfc1
        cfc2
        cfc3
        cfc4
        cfc5
        buffer
    end

    methods
        % Construtor da classe
        function obj = GAPI(popSize, numGenerations, crossoverRate, mutationRate, paramBounds, Rs, Rr, Lm, Lls, buffer)
            obj.popSize = popSize;
            obj.numGenerations = numGenerations;
            obj.crossoverRate = crossoverRate;
            obj.mutationRate = mutationRate;
            obj.paramNum = 6;
            obj.pop = zeros(popSize, 6);
            obj.paramBounds = paramBounds;
            obj.lowerBounds = repmat(paramBounds(:,1)', popSize, 1);
            obj.upperBounds = repmat(paramBounds(:,2)', popSize, 1);
            obj.fitnessScores = zeros(popSize, 1);
            obj.bestFitnessHistory = zeros(numGenerations, 1);
            obj.averageFitnessHistory = zeros(numGenerations, 1);
            obj.bestIndividuals = zeros(numGenerations, length(paramBounds));
            obj.generation = 0;

            obj.Rs = Rs;
            obj.Rr = Rr;
            obj.Lm = Lm;
            obj.Lls = Lls;
            obj.Llr = Lls;

            obj.buffer = buffer;

            % var_filter = obj.Initilize_filters();
            var_filter = [0, 0, 0, 0, 0];
            obj.cfc1 = var_filter(1);
            obj.cfc2 = var_filter(2);
            obj.cfc3 = var_filter(3);
            obj.cfc4 = var_filter(4);
            obj.cfc5 = var_filter(5);

            obj.CostFunction = @obj.f_motor;
        end

        function var_filter = Initilize_filters(obj)
            %% Parâmetros da Simulação
            % Frequencia de amostragem do sinal (Hz)
            fs = 10000;
            % Periodo de amostragem do sinal (s)
            Ts = 1/fs;
            % Numero de partes que o intervalo discreto e dividido
            p = 1000;
            % Passo de amostragem continuo (s)
            h = Ts/p;
            %% Parâmetros do filtro de corrente (continuo)
            % ordem -> 2, discretização por Backward, fc = 600Hz, ganho unitário
            % Frequência de Corte (Hz)
            fc = 600;
            % Frequência natural
            wn = 2*pi*fc;
            % Coeficiente de amortecimento
            zeta = 0.707;

            % Função de Transferência do filtro (Domínio da frequência)
            %
            %    Y(s)               wn^2
            %   ______  =  ________________________
            %
            %    X(s)       s^2 + zeta*wn*s + wn^2
            %
            Ft_filtro = tf(wn^2,[1 2*zeta*wn wn^2]);
            % Função de transferência Discretizada
            %
            %    Y(Z)        a0 + a1*z^(-1) + a2*z^(-2)
            %   ______  =  ______________________________
            %
            %    X(Z)        b0 + b1*z^(-1) + b2*z^(-2)
            %
            Ft_filtro_discretizada = c2d(Ft_filtro,h,'tustin');
            % Coleta dos coeficientes
            %
            %    Y(K) = (a0/b0)*X(K) + (a1/b0)*X(K-1) + (a2/b0)*X(K-2) - (b1/b0)*Y(K-1) - (b2/b0)*Y(K-2)
            %
            [Ft_num,Ft_den] = tfdata(Ft_filtro_discretizada,'v');

            a0 = Ft_num(1);
            a1 = Ft_num(2);
            a2 = Ft_num(3);

            b0 = Ft_den(1);
            b1 = Ft_den(2);
            b2 = Ft_den(3);
            % Reescrevendo a equação
            %
            %    Y(K) = cfc1*X(K) + cfc2*X(K-1) + cfc3*X(K-2) + cfc4*Y(K-1) + cfc5*Y(K-2)
            %
            cf1 = a0/b0;
            cf2 = a1/b0;
            cf3 = a2/b0;
            cf4 = -b1/b0;
            cf5 = -b2/b0;

            var_filter = [cf1, cf2, cf3, cf4, cf5];
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
                estocSize = selectionArgs{2};
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

        function penalty = inf_if_zero(obj, value)
            if value == 0
                penalty = inf;
            else
                penalty = value;
            end
        end
        %% Funções de Custo
        function [cost] = f_motor_2(obj, x)
            %% Função auxiliar para penalização
            %% Parâmetros da Simulação
            % Frequencia de amostragem do sinal (Hz)
            fs = 10000;
            % Periodo de amostragem do sinal (s)
            Ts = 1/fs;
            % Numero de partes que o intervalo discreto e dividido
            p = 10;
            % Passo de amostragem continuo (s)
            h = Ts/p;
            % Tempo de Simulação (s)
            T_simulacao = 0.4;
            % Número de pontos gerados em uma simulação completa (Discreto)
            Np = T_simulacao/Ts;

            %% Parâmetros da Rede
            % Frequência da rede (Hz)
            frequencia = 60;

            %% Parâmetros do Motor

            % Numero de polos
            Polos = 8;
            % Indutância dos enrolamentos do rotor
            Lr = obj.Lm + obj.Llr;
            % Coefieciente de Inércia do motor
            J =  0.1633;
            % Coeficiente de atrito do eixo do motor
            K = 0.12;
            % Potência do motor (W)
            P = 4*736;

            %% Parâmetros do Modelo Matemático
            % Velocidade sincrona em radianos elétricos por segundo
            w_rede = 2*pi*frequencia;

            % Reatância de magnetização
            Xm = w_rede*obj.Lm;
            % Reatância de dispersão do estator
            Xls = w_rede*obj.Lls;
            % Reatância de dispersão do rotor
            Xlr = w_rede*obj.Llr;
            Xml = 1/(1/Xm + 1/Xls + 1/Xlr);

            % Constantes do modelo mecânico
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

            %% Condições iniciais
            % Velocidade inicial
            wr = 0;

            %% Inicialização das variaveis
            Fqs = 0;
            Fds = 0;
            Fqr = 0;
            Fdr = 0;
            Ids = 0;
            Ids_ant = 0;
            Iqs = 0;
            lambda_dr_est = 0;
            theta_park = 0;
            UI_w = 0;
            UI_id = 0;
            UI_iq = 0;
            t = 0:Ts:Np*Ts-Ts;
            t_continuo = 0;
            cost = 0;
            custo_erros = 0;
            e_w_anterior = 0;
            e_iq_anterior = 0;
            e_id_anterior = 0;

            Ia_k1 = 0;
            Ia_k0 = 0;

            Ia_filtrado_k1 = 0;
            Ia_filtrado = 0;

            Ids_filtrado_atrasado = 0;
            Iqs_filtrado_atrasado = 0;


            Ib_k1 = 0;
            Ib_k0 = 0;

            Ib_filtrado_k1 = 0;
            Ib_filtrado = 0;


            Ic_k1 = 0;
            Ic_k0 = 0;

            Ic_filtrado_k1 = 0;
            Ic_filtrado = 0;


            % Constantes Auxiliares (redução de cálculos)
            pi23 = pi*2/3;
            c23 = 2/3;

            %% Torque de Carga
            Tn = P*Polos/(2*w_rede);
            Tl =  Tn/2 *10*0* ((t-0.6).*(t >= 0.6) - (t-0.7).*(t >= 0.7));

            %% Corrente Id de referência
            lambda_nonminal = 127/(2*pi*frequencia)/(obj.Lm);
            Ids_ref = 20*lambda_nonminal * ((t-0).*(t >= 0) - (t-0.05).*(t >= 0.05));

            %% Velocidade de Referência
            w_nom_ref = 2*pi*60;
            w_ref = 4*w_nom_ref * ((t-0.05).*(t >= 0.05) - (t-0.3).*(t >= 0.3));

            %% buffer de atraso
            % Inicializando buffers para armazenar os valores anteriores
            buffer_size = obj.buffer;  % Tamanho do atraso (número de amostras)
            Iqs_buffer = zeros(1, buffer_size);
            Ids_buffer = zeros(1, buffer_size);

            %% Loop Simulação Motor
            for k = 1:Np
                %% Estimador de fluxo rotórico e orientação do sist. de referência

                lambda_dr_est = lambda_dr_est*( ...
                    (2*Lr-Ts*obj.Rr)/(2*Lr+Ts*obj.Rr)) + Ids_filtrado_atrasado*( ...
                    (obj.Lm*obj.Rr*Ts)/(2*Lr+obj.Rr*Ts)) + Ids_ant*( ...
                    (obj.Lm*obj.Rr*Ts)/(2*Lr+obj.Rr*Ts) ...
                    );
                Ids_ant = Ids_filtrado_atrasado;

                if(lambda_dr_est > 0.1)
                    wsl = (obj.Lm*obj.Rr*Iqs_filtrado_atrasado)/(Lr*lambda_dr_est);
                else
                    wsl = 0;
                end

                wr_est = wr + wsl;
                w = wr_est;
                theta_park = theta_park + Ts*w;

                %% Calculando as correntes Ia Ib e Ic

                %Velocidade
                e_w = w_ref(k) - wr;
                UI_w = UI_w + e_w*Ts;
                U_w = KP_w*e_w + KI_w*UI_w;
                iqs_ref = U_w;

                %Servos de corrente
                e_id = Ids_ref(k) - Ids_filtrado_atrasado*(1 + 0.1*(2 * rand() - 1));
                UI_id = UI_id + e_id*Ts;
                U_Id = KP_id*e_id + KI_id*UI_id;

                if(U_Id >= 127*sqrt(2))
                    U_Id = 127*sqrt(2);
                end
                if(U_Id <= -127*sqrt(2))
                    U_Id = -127*sqrt(2);
                end

                e_iq = iqs_ref - Iqs_filtrado_atrasado*(1 + 0.1*(2 * rand() - 1));
                UI_iq = UI_iq*e_iq*Ts;
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

                for ksuper=1:p
                    Fqm = Xml/Xls*Fqs + Xml/Xlr*Fqr;
                    Fdm = Xml/Xls*Fds + Xml/Xlr*Fdr;

                    Fqs = Fqs + h*w_rede*(Vq - w/w_rede*Fds - obj.Rs/Xls*(Fqs-Fqm));
                    Fds = Fds + h*w_rede*(Vd + w/w_rede*Fqs - obj.Rs/Xls*(Fds-Fdm));

                    Fqr = Fqr - h*w_rede*((w-wr)*Fdr/w_rede + obj.Rr/Xlr*(Fqr-Fqm));
                    Fdr = Fdr - h*w_rede*((wr-w)*Fqr/w_rede + obj.Rr/Xlr*(Fdr-Fdm));

                    Iqs = (Fqs-Fqm)/Xls;
                    Ids = (Fds-Fdm)/Xls;

                    % Torque Eletromagnético
                    Te = 3/2*Polos/2*1/w_rede*(Fds*Iqs - Fqs*Ids);

                    % Solução mecânica
                    wr = wr + h*(A*wr + B1*Te + B2*Tl(k));
                    %
                    % % Transformada inversa de park para as correntes dq0 -> abc
                    % % Eixo q alinhado com a fase a
                    % Ia = Ids*sin(theta_park) + Iqs*cos(theta_park);
                    % Ib = Ids*sin(theta_park - pi23) + Iqs*cos(theta_park - pi23);
                    % Ic = Ids*sin(theta_park + pi23) + Iqs*cos(theta_park + pi23);
                    %
                    % amplitude_ruido = 0.2;
                    %
                    % % Frequências para os ruídos de baixa e alta frequência
                    % frequencia_baixa = 50;   % 50 Hz para simular ruído de baixa frequência
                    % frequencia_alta = 1000;  % 1000 Hz para simular interferência de alta frequência
                    %
                    % % Ruído branco gaussiano
                    % ruido_branco_Ia = amplitude_ruido * 0.5 * randn(size(Ia));
                    % ruido_branco_Ib = amplitude_ruido * 0.5 * randn(size(Ib));
                    % ruido_branco_Ic = amplitude_ruido * 0.5 * randn(size(Ic));
                    %
                    % % Ruído de baixa frequência (50 Hz)
                    % ruido_baixa_frequencia_Ia = amplitude_ruido * 0.3 * sin(2*pi*frequencia_baixa*t_continuo);
                    % ruido_baixa_frequencia_Ib = amplitude_ruido * 0.3 * sin(2*pi*frequencia_baixa*t_continuo);
                    % ruido_baixa_frequencia_Ic = amplitude_ruido * 0.3 * sin(2*pi*frequencia_baixa*t_continuo);
                    %
                    % % Ruído de alta frequência (1000 Hz)
                    % ruido_alta_frequencia_Ia = amplitude_ruido * 0.2 * sin(2*pi*frequencia_alta*t_continuo);
                    % ruido_alta_frequencia_Ib = amplitude_ruido * 0.2 * sin(2*pi*frequencia_alta*t_continuo);
                    % ruido_alta_frequencia_Ic = amplitude_ruido * 0.2 * sin(2*pi*frequencia_alta*t_continuo);
                    %
                    % % Combinação dos ruídos para os sinais medidos
                    % Ia_medido = Ia + ruido_branco_Ia + ruido_baixa_frequencia_Ia + ruido_alta_frequencia_Ia;
                    % Ib_medido = Ib + ruido_branco_Ib + ruido_baixa_frequencia_Ib + ruido_alta_frequencia_Ib;
                    % Ic_medido = Ic + ruido_branco_Ic + ruido_baixa_frequencia_Ic + ruido_alta_frequencia_Ic;
                    %
                    % %% Filtro Analógico de sinais
                    % % Correntes ABC
                    % % 1º Etapa -> Armazenamento das variáveis anteriores (abordagem sem o
                    % % uso de vetores)
                    % % IA
                    % Ia_k2 = Ia_k1;
                    % Ia_k1 = Ia_k0;
                    % Ia_k0 = Ia_medido;
                    % Ia_filtrado_k2 = Ia_filtrado_k1;
                    % Ia_filtrado_k1 = Ia_filtrado;
                    % % IB
                    % Ib_k2 = Ib_k1;
                    % Ib_k1 = Ib_k0;
                    % Ib_k0 = Ib_medido;
                    % Ib_filtrado_k2 = Ib_filtrado_k1;
                    % Ib_filtrado_k1 = Ib_filtrado;
                    % % IC
                    % Ic_k2 = Ic_k1;
                    % Ic_k1 = Ic_k0;
                    % Ic_k0 = Ic_medido;
                    % Ic_filtrado_k2 = Ic_filtrado_k1;
                    % Ic_filtrado_k1 = Ic_filtrado;
                    %
                    % % Equações de filtro
                    % Ia_filtrado = obj.cfc1*Ia_medido + obj.cfc2*Ia_k1 + obj.cfc3*Ia_k2 + obj.cfc4*Ia_filtrado_k1 + obj.cfc5*Ia_filtrado_k2;
                    % Ib_filtrado = obj.cfc1*Ib_medido + obj.cfc2*Ib_k1 + obj.cfc3*Ib_k2 + obj.cfc4*Ib_filtrado_k1 + obj.cfc5*Ib_filtrado_k2;
                    % Ic_filtrado = obj.cfc1*Ic_medido + obj.cfc2*Ic_k1 + obj.cfc3*Ic_k2 + obj.cfc4*Ic_filtrado_k1 + obj.cfc5*Ic_filtrado_k2;
                end

                % Ids_filtrado = c23*(Ia_filtrado*sin(theta_park) + ...
                %     Ib_filtrado*sin(theta_park - pi23) + ...
                %     Ic_filtrado*sin(theta_park + pi23));
                % Iqs_filtrado = c23*(Ia_filtrado*cos(theta_park) + ...
                %     Ib_filtrado*cos(theta_park - pi23) + ...
                %     Ic_filtrado*cos(theta_park + pi23));
                % %
                Iqs_filtrado = Iqs;
                Ids_filtrado = Ids;


                % Desloca o buffer para a direita (atualiza o atraso)
                Iqs_buffer = [Iqs_filtrado, Iqs_buffer(1:end-1)];
                Ids_buffer = [Ids_filtrado, Ids_buffer(1:end-1)];

                % Extrai os valores atrasados do buffer
                Iqs_filtrado_atrasado = Iqs_buffer(end);  % Valor com atraso de buffer_size amostras
                Ids_filtrado_atrasado = Ids_buffer(end);  % Valor com atraso de buffer_size amostras

                %% Penalização para zeros nos ganhos
                cost_KP_w = obj.inf_if_zero(KP_w)/1000;
                cost_KI_w = obj.inf_if_zero(KI_w)/1000;
                cost_KP_id = obj.inf_if_zero(KP_id)/1000;
                cost_KI_id = obj.inf_if_zero(KI_id)/1000;
                cost_KP_iq = obj.inf_if_zero(KP_iq)/1000;
                cost_KI_iq = obj.inf_if_zero(KI_iq)/1000;

                % % Cálculo das derivadas dos erros
                % de_w = (e_w - e_w_anterior) / Ts;
                % de_iq = (e_iq - e_iq_anterior) / Ts;
                % de_id = (e_id - e_id_anterior) / Ts;
                %
                % e_w_anterior = e_w;
                % e_iq_anterior = e_iq;
                % e_id_anterior = e_id;
                %
                % % Penalização pela taxa de variação do erro
                % K_derivada = 1;
                % penalizacao_derivada = K_derivada * (abs(de_w) + abs(de_iq) + abs(de_id));  % K_derivada é um ganho
                %
                % % Penalização por erro persistente (indicando longo tempo de acomodação)
                % K_tempo = 1;
                % penalizacao_tempo_acomodacao = K_tempo * sum(abs(e_w) + abs(e_iq) + abs(e_id));  % K_tempo é um ganho
                %
                % Adicionar à função de custo
                custo_erros = custo_erros + Ts * (abs(e_w) + abs(e_iq) + abs(e_id)) ;
                %
                % custo_ganhos = cost_KP_w + cost_KI_w + cost_KP_id + cost_KI_id + cost_KP_iq + cost_KI_iq;

                cost = custo_erros;

            end

        end

        function [cost] = f_motor(obj, x)
            %% Parâmetros da Simulação

            f = 10000;                                     % Frequencia de amostragem do sinal
            Tsc = 1/f;                                     % Periodo de amostragem do sinal
            p = 10;                                        % Numero de partes que o intervalo discreto e dividido
            h = Tsc/p;                                     % Passo de amostragem continuo
            Tsimu = 0.5;                                    % Tempo de Simulação
            Np = Tsimu/Tsc;                                % Número de Pontos (vetores)

            %% Parâmetros do Motor

            Polos = 8;                                   % Numero de polos
            frequencia = 60;                             % Frequência elétrica nominal
            Lr = obj.Lm + obj.Llr;                               % Indutância dos enrolamentos do rotor
            J =  0.1633;                                 % Inércia do motor
            K = 0.12;                                    % Coeficiente de atrito do eixo do motor
            weles = 2*pi*frequencia;                     % Velocidade sincrona em radianos elétricos por segundo
            wr = 0;                                      % Velocidade inicial
            P = 4*736;                                   % Potência do motor

            %% Reatâncias para espaço de estados

            Xm = weles*obj.Lm;                               % Reatância de magnetização
            Xls = weles*obj.Lls;                             % Reatância de dispersão do estator
            Xlr = weles*obj.Llr;                             % Reatância de dispersão do rotor
            Xml = 1/(1/Xm + 1/Xls + 1/Xlr);

            %% Constantes para solução mecânica do motor

            A =  - K/J;
            B1 =   Polos/(2*J);
            B2 = - Polos/(2*J);

            %% Ganhos Controladores

            KP_w = x(1);
            KI_w =  x(2);

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
            cost = 0;
            t = 0:Tsc:Np*Tsc-Tsc;

            %% Torque de Carga

            Tn = P*Polos/(2*weles);

            Tl =  Tn*0.75 *20* ((t-0.3).*(t >= 0.3) - (t-0.5).*(t >= 0.5));

            %% Corrente Id de referência

            lambda_nonminal = 127/(2*pi*frequencia)/(obj.Lm);

            Ids_ref = 20*lambda_nonminal * ((t-0).*(t >= 0) - (t-0.05).*(t >= 0.05));

            %% Velocidade de Referência

            w_nom_ref = 2*pi*60;

            w_ref = 4*w_nom_ref * ((t-0.05).*(t >= 0.05) - (t-0.3).*(t >= 0.3));

            %% Loop Simulação Motor
            for k = 1:Np
                %% Estimador de fluxo rotórico e orientação do sist. de referência

                lambda_dr_est = lambda_dr_est*((2*Lr-Tsc*obj.Rr)/(2*Lr+Tsc*obj.Rr)) + Ids*((obj.Lm*obj.Rr*Tsc)/(2*Lr+obj.Rr*Tsc)) + Ids_ant*((obj.Lm*obj.Rr*Tsc)/(2*Lr+obj.Rr*Tsc));
                Ids_ant = Ids;

                if(lambda_dr_est > 0.1)
                    wsl = (obj.Lm*obj.Rr*Iqs)/(Lr*lambda_dr_est);
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
                e_id = Ids_ref(k) - Ids;
                UI_id = UI_id + e_id*Tsc;
                U_Id = KP_id*e_id + KI_id*UI_id;

                if(U_Id >= 127*sqrt(2))
                    U_Id = 127*sqrt(2);
                end
                if(U_Id <= -127*sqrt(2))
                    U_Id = -127*sqrt(2);
                end

                e_iq = iqs_ref - Iqs;
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

                cost = cost + Tsc*(abs(e_w) + abs(e_iq) + abs(e_id));

            end
            % Penalização para zeros nos ganhos
            cost_KP_w = obj.inf_if_zero(KP_w)/1000;
            cost_KI_w = obj.inf_if_zero(KI_w)/1000;
            cost_KP_id = obj.inf_if_zero(KP_id)/1000;
            cost_KI_id = obj.inf_if_zero(KI_id)/1000;
            cost_KP_iq = obj.inf_if_zero(KP_iq)/1000;
            cost_KI_iq = obj.inf_if_zero(KI_iq)/1000;
            custo_ganhos = cost_KP_w + cost_KI_w + cost_KP_id + cost_KI_id + cost_KP_iq + cost_KI_iq;
            
            cost = cost + custo_ganhos;
        end
    end
end