obj.Lm = 0.0367;
obj.Llr = 0.0028;
obj.Rs = 0.6759;
obj.Rr = 0.2759;
obj.Lls = 0.0028;

x = [33.2169413, 270.54687,  9809.72592 ,  36107.8129];


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
T_simulacao = 10;                                     
% Número de pontos gerados em uma simulação completa (Discreto)
Np = T_simulacao/Ts;                                
% Número de pontos gerados em uma simulação completa (Contínuo)
Npc = p*(Np + 1);

%% Parâmetros da Rede
% Frequência da rede (Hz)
frequencia = 60;
% Tensão de pico de fase (V)
V_pico_fase = 180;
% Tensão Rms Fase-Fase (V)
V_rms_fase_fase = 220;
% Tensão do Barramento CC (Vdc)
Vdc = V_rms_fase_fase*sqrt(2);

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
% Velocidade do sistema síncrono
w_park = 0;

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
cfc1 = a0/b0;
cfc2 = a1/b0;
cfc3 = a2/b0;
cfc4 = -b1/b0;
cfc5 = -b2/b0;

%% Ganhos Controladores
KP_w = x(1);
KI_w = x(2);

KP_id = x(3);
KI_id = x(4);

KP_iq = x(3);
KI_iq = x(4);

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

Ia_k2 = 0;
Ia_k1 = 0;
Ia_k0 = 0;
Ia_filtrado_k2 = 0;
Ia_filtrado_k1 = 0;
Ia_filtrado = 0;

Ids_filtrado_atrasado = 0;
Iqs_filtrado_atrasado = 0;

Ib_k2 = 0;
Ib_k1 = 0;
Ib_k0 = 0;
Ib_filtrado_k2 = 0;
Ib_filtrado_k1 = 0;
Ib_filtrado = 0;

Ic_k2 = 0;
Ic_k1 = 0;
Ic_k0 = 0;
Ic_filtrado_k2 = 0;
Ic_filtrado_k1 = 0;
Ic_filtrado = 0;

Ia_d_k2 = 0;
Ia_d_k1 = 0;
Ia_d_k0 = 0;
Ia_d_filtrado_k2 = 0;
Ia_d_filtrado_k1 = 0;
Ia_d_filtrado = 0;

Ib_d_k2 = 0;
Ib_d_k1 = 0;
Ib_d_k0 = 0;
Ib_d_filtrado_k2 = 0;
Ib_d_filtrado_k1 = 0;
Ib_d_filtrado = 0;

Ic_d_k2 = 0;
Ic_d_k1 = 0;
Ic_d_k0 = 0;
Ic_d_filtrado_k2 = 0;
Ic_d_filtrado_k1 = 0;
Ic_d_filtrado = 0;

Iqs_filtrado = 0;
Ids_filtrado = 0;

Vd_continuo_k2 = 0;
Vd_continuo_k1 = 0;
Vd_continuo_k0 = 0;
Vd_continuo_filtrado_k2 = 0;
Vd_continuo_filtrado_k1 = 0;
Vd_continuo_filtrado = 0;

Vq_continuo_k2 = 0;
Vq_continuo_k1 = 0;
Vq_continuo_k0 = 0;
Vq_continuo_filtrado_k2 = 0;
Vq_continuo_filtrado_k1 = 0;
Vq_continuo_filtrado = 0;

Vd_k2 = 0;
Vd_k1 = 0;
Vd_k0 = 0;
Vd_filtrado_k2 = 0;
Vd_filtrado_k1 = 0;
Vd_filtrado = 0;

Vq_k2 = 0;
Vq_k1 = 0;
Vq_k0 = 0;
Vq_filtrado_k2 = 0;
Vq_filtrado_k1 = 0;
Vq_filtrado = 0;

% Constantes Auxiliares (redução de cálculos)
pi23 = pi*2/3;
c23 = 2/3;
p2 = p/2;
pi180 = 180/pi;
Vdc23 = 2/3*Vdc;
r23 = sqrt(2/3);
Ts2 = Ts/2;
P34 = 3*Polos/4;
rpm = 2/Polos*60/(2*pi);

%% Torque de Carga
Tn = P*Polos/(2*w_rede);
Tl =  Tn/2 *10*0* ((t-0.6).*(t >= 0.6) - (t-0.7).*(t >= 0.7));

%% Corrente Id de referência
lambda_nonminal = 127/(2*pi*frequencia)/(obj.Lm);
Ids_ref = 10*lambda_nonminal * ((t-0).*(t >= 0) - (t-0.1).*(t >= 0.1));

%% Velocidade de Referência
w_nom_ref = 2*pi*60;
w_ref = 4*w_nom_ref * ((t-0.1).*(t >= 0.1) - (t-0.35).*(t >= 0.35));

%% buffer de atraso
% Inicializando buffers para armazenar os valores anteriores
buffer_size = 10;  % Tamanho do atraso (número de amostras)
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

    iqs_ref = U_w;  % K_predict é um ganho preditivo

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

    %% Calculando as Tensões Va Vb e Vc
    Valfa = Vq*cos(theta_park) + Vd*sin(theta_park);
    Vbeta = -Vq*sin(theta_park) + Vd*cos(theta_park);

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
        % Ia_filtrado = cfc1*Ia_medido + cfc2*Ia_k1 + cfc3*Ia_k2 + cfc4*Ia_filtrado_k1 + cfc5*Ia_filtrado_k2;
        % Ib_filtrado = cfc1*Ib_medido + cfc2*Ib_k1 + cfc3*Ib_k2 + cfc4*Ib_filtrado_k1 + cfc5*Ib_filtrado_k2;
        % Ic_filtrado = cfc1*Ic_medido + cfc2*Ic_k1 + cfc3*Ic_k2 + cfc4*Ic_filtrado_k1 + cfc5*Ic_filtrado_k2;
        % 
        % 
        % k_aux = p*k + ksuper;
        % 
        % t_continuo = t_continuo + h;
        % t_continuo_vetor(k_aux) = t_continuo;
    end
   
    % % Equações de filtro
    % % IA
    % Ia_d_k2 = Ia_d_k1;
    % Ia_d_k1 = Ia_d_k0;
    % Ia_d_k0 = Ia_filtrado;
    % Ia_d_filtrado_k2 = Ia_d_filtrado_k1;
    % Ia_d_filtrado_k1 = Ia_d_filtrado;
    % % IB
    % Ib_d_k2 = Ib_d_k1;
    % Ib_d_k1 = Ib_d_k0;
    % Ib_d_k0 = Ib_filtrado;
    % Ib_d_filtrado_k2 = Ib_d_filtrado_k1;
    % Ib_d_filtrado_k1 = Ib_d_filtrado;
    % % IC
    % Ic_d_k2 = Ic_d_k1;
    % Ic_d_k1 = Ic_d_k0;
    % Ic_d_k0 = Ic_filtrado;
    % Ic_d_filtrado_k2 = Ic_d_filtrado_k1;
    % Ic_d_filtrado_k1 = Ic_d_filtrado;
    % 
    % Ia_d_filtrado = cf1*Ia_filtrado + cf2*Ia_d_k1 + cf3*Ia_d_k2 + cf4*Ia_d_filtrado_k1 + cf5*Ia_d_filtrado_k2;
    % Ib_d_filtrado = cf1*Ib_filtrado + cf2*Ib_d_k1 + cf3*Ib_d_k2 + cf4*Ib_d_filtrado_k1 + cf5*Ib_d_filtrado_k2;
    % Ic_d_filtrado = cf1*Ic_filtrado + cf2*Ic_d_k1 + cf3*Ic_d_k2 + cf4*Ic_d_filtrado_k1 + cf5*Ic_d_filtrado_k2;
    

    % Ids_filtrado = c23*(Ia_filtrado*sin(theta_park) + ...
    %     Ib_filtrado*sin(theta_park - pi23) + ...
    %     Ic_filtrado*sin(theta_park + pi23));
    % Iqs_filtrado = c23*(Ia_filtrado*cos(theta_park) + ...
    %     Ib_filtrado*cos(theta_park - pi23) + ...
    %     Ic_filtrado*cos(theta_park + pi23));

    Ids_filtrado = Ids;
    Iqs_filtrado = Iqs;

    
    % Desloca o buffer para a direita (atualiza o atraso)
    Iqs_buffer = [Iqs_filtrado, Iqs_buffer(1:end-1)];
    Ids_buffer = [Ids_filtrado, Ids_buffer(1:end-1)];

    % Extrai os valores atrasados do buffer
    Iqs_filtrado_atrasado = Iqs_buffer(end);  % Valor com atraso de buffer_size amostras
    Ids_filtrado_atrasado = Ids_buffer(end);  % Valor com atraso de buffer_size amostras

    cost = cost + Ts * (2*abs(e_w) + abs(e_iq) + abs(e_id));


    wr_vetor(k) = wr;
    ids_vetor(k) = Ids;
    Te_vetor(k) = Te;

    Va_vetor(k) = Va;
    Vb_vetor(k) = Vb;
    Vc_vetor(k) = Vc;

    Valpha_vetor = Valfa;
    Vbeta_vetor = Vbeta;

    Vd_vetor(k) = Vd;
    Vq_vetor(k) = Vq;

    Iqs_vetor(k) = Iqs;
    Iqs_filtrado_vetor(k) = Iqs_filtrado;
    Iqs_atrasado_vetor(k) = Iqs_filtrado_atrasado;
    Ids_vetor(k) = Ids;
    Ids_filtrado_vetor(k) = Ids_filtrado;
    Ids_atrasado_vetor(k) = Ids_filtrado_atrasado;
end

t = 0:Ts:Np*Ts-Ts;

figure
% Plotando os dados
plot(t,wr_vetor,t,w_ref); % tom de cinza escuro
legend('Wr','Ref')
% 
% figure
% % Plotando os dados
% plot(t,ids_vetor); % tom de cinza escuro
% legend('ID')
% 
% figure
% % Plotando os dados
% plot(t,Te_vetor,t, Tl); % tom de cinza escuro
% legend('Torque')
% 
% figure
% % Plotando os dados
% plot(t,Va_vetor,t, Vb_vetor, t, Vc_vetor); % tom de cinza escuro
% legend('Va', 'Vb', 'Vc')
% 
% figure
% % Plotando os dados
% plot(t,Ia_vetor,t, Ib_vetor, t, Ic_vetor); % tom de cinza escuro
% legend('Ia', 'Ib', 'Ic')
% 
% figure
% % Plotando os dados
% plot(t,Ia_medido_vetor,t, Ib_medido_vetor, t, Ic_medido_vetor); % tom de cinza escuro
% legend('Ia_medido', 'Ib_medido', 'Ic_medido')
% 
% figure
% % Plotando os dados
% plot(t,Ia_filtrado_vetor,t, Ib_filtrado_vetor, t, Ic_filtrado_vetor); % tom de cinza escuro
% legend('Ia_filtrado', 'Ib_filtrado', 'Ic_filtrado')
% 
% 
figure
% Plotando os dados
plot(t,Iqs_vetor,t, Iqs_filtrado_vetor); % tom de cinza escuro
legend('Iqs', 'Iqs filtrado')
% 
% figure
% % Plotando os dados
% plot(t,Ids_vetor,t, Ids_filtrado_vetor); % tom de cinza escuro
% legend('Ids', 'Ids filtrado')

figure
% Plotando os dados
plot(t,Iqs_filtrado_vetor,t, Iqs_atrasado_vetor); % tom de cinza escuro
legend('Iqs', 'Iqs filtrado')

% figure
% % Plotando os dados
% plot(t,Ids_filtrado_vetor,t, Ids_atrasado_vetor); % tom de cinza escuro
% legend('Ids filtrado', 'Ids atrasado')