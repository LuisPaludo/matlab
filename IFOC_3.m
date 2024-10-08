% close all;
clear all;
close all;
clc;

%% Parâmetros da Simulação

f = 10000;                                     % Frequencia de amostragem do sinal
Tsc = 1/f;                                     % Periodo de amostragem do sinal
p = 10;                                        % Numero de partes que o intervalo discreto e dividido
h = Tsc/p;                                     % Passo de amostragem continuo
Tsimu = 1;                                    % Tempo de Simulação
Np = Tsimu/Tsc;                                % Número de Pontos (vetores)

%% Parâmetros do Motor

Polos = 8;                                   % Numero de polos
Vf = 127;
frequencia = 60;                             % Frequência elétrica nominal
Rs = 0.6759;                                 % Resistência do estator
Lls = 0.00280;                               % Indutância de dispersão dos enrolamentos do estator
Lm = 0.0387;                                 % Indutância mútua dos enrolamentos - estator/estator - estator/rotor - rotor/rotor
Llr = 0.00280;                               % Indutância de dispersão dos enrolamentos do rotor
Lr = Lm + Llr;                               % Indutância dos enrolamentos do rotor
Ls = Lm + Lls;                               % Indutância dos enrolamentos do estator
Rr = 0.2615;                                 % Resistencia do rotor
J =  0.1633;                                 % Inércia do motor
K = 0.12;                                    % Coeficiente de atrito do eixo do motor
weles = 2*pi*frequencia;                     % Velocidade sincrona em radianos elétricos por segundo
wr = 0;                                      % Velocidade inicial
P = 4*736;                                   % Potência do motor

%% Reatâncias para espaço de estados

Xm = weles*Lm;                               % Reatância de magnetização
Xls = weles*Lls;                             % Reatância de dispersão do estator
Xlr = weles*Llr;                             % Reatância de dispersão do rotor
Xml = 1/(1/Xm + 1/Xls + 1/Xlr);

%% Constantes para solução mecânica do motor

A =  - K/J;
B1 =   Polos/(2*J);
B2 = - Polos/(2*J);

%% Ganhos Controladores
% 1696.8194              0      1.3668319      259.29596       13.96912              0
% Geração: 100 Melhor Valor: 1.3562 Melhor Individuo: 1140.7216              0      2.7360757      285.60373     0.40253971              0
% Geração: 100 Melhor Valor: 7.8133 Melhor Individuo: 8364.777      31.419262       9.022495      289.80748      85.198494       87.75967
% Geração: 100 Melhor Valor: 22.1557 Melhor Individuo: 1043.4093       1.312594      71.357879      1052.7372       49.78579      80.287684
KP_w = 347.0288     ;
KI_w =    0  ;

KP_id = 86.85888 ;
KI_id =    0;

KP_iq = 12.4916   ;
KI_iq = 0;

%% Inicialização das variaveis

Fqs = 0;
Fds = 0;
Fqr = 0;
Fdr = 0;
Ids = 0;
Ids_ant = 0;
Iqs = 0;
Vq = 0;
Vd = 0;
lambda_dr_est = 0;
theta = 0;
UI_w = 0;
UI_id = 0;
UI_iq = 0;
ialpha_vetor = zeros(1,Np);
ibeta_vetor = zeros(1,Np);
iqs_vetor = zeros(1,Np);
ids_vetor = zeros(1,Np);
Ids_ref_vetor = zeros(1,Np);
Te_vetor = zeros(1,Np);
wrpm_vetor = zeros(1,Np);
wr_vetor = zeros(1,Np);
w_ref_vetor = zeros(1,Np);
t = 0:Tsc:Np*Tsc-Tsc;
Va_vetor = zeros(1,Np);
Vb_vetor = zeros(1,Np);
Vc_vetor = zeros(1,Np);
Vd_vetor = zeros(1,Np);
Vd_est_vetor = zeros(1,Np);
Vq_vetor = zeros(1,Np);
Vq_est_vetor = zeros(1,Np);
lambda_qr_vetor = zeros(1,Np);
lambda_dr_vetor = zeros(1,Np);
lambda_dr_est_vetor = zeros(1,Np);
Ia_vetor = zeros(1,Np);
Ib_vetor = zeros(1,Np);
Ic_vetor = zeros(1,Np);
U_w_vetor = zeros(1,Np);
U_iq_vetor = zeros(1,Np);
U_id_vetor = zeros(1,Np);
e_w_vetor = zeros(1,Np);
e_iq_vetor = zeros(1,Np);
e_id_vetor = zeros(1,Np);
iqs_atrasado_vetor = zeros(1,Np);
ids_atrasado_vetor = zeros(1,Np);

%% Torque de Carga

Tn = P*Polos/(2*weles);

Tl =  Tn*0.75 *100* ((t-0.6).*(t >= 0.6) - (t-0.61).*(t >= 0.61));

%% Corrente Id de referência

lambda_nonminal = 127/(2*pi*frequencia)/(Lm);

Ids_ref = 10*lambda_nonminal * ((t-0).*(t >= 0) - (t-0.1).*(t >= 0.1));

%% Velocidade de Referência

w_nom_ref = 2*pi*60;

w_ref = 2*w_nom_ref * ((t-0.1).*(t >= 0.1) - (t-0.6).*(t >= 0.6));

%% buffer de atraso
% Inicializando buffers para armazenar os valores anteriores
buffer_size = 10;  % Tamanho do atraso (número de amostras)
buffer_index = 1;
Iqs_buffer = zeros(1, buffer_size);
Ids_buffer = zeros(1, buffer_size);

%% Loop Simulação Motor
for k = 1:Np
    %% Estimador de fluxo rotórico e orientação do sist. de referência

    lambda_dr_est = lambda_dr_est*((2*Lr-Tsc*Rr)/(2*Lr+Tsc*Rr)) + Ids*((Lm*Rr*Tsc)/(2*Lr+Rr*Tsc)) + Ids_ant*((Lm*Rr*Tsc)/(2*Lr+Rr*Tsc));
    Ids_ant = Ids;

    if(lambda_dr_est > 0.1)
        wsl = (Lm*Rr*Iqs)/(Lr*lambda_dr_est);
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

        % Recalcule Va, Vb e Vc com os novos Vd e Vq se desejar
        Valfa = Vq*cos(theta) + Vd*sin(theta);
        Vbeta = -Vq*sin(theta) + Vd*cos(theta);

        Va = Valfa;
        Vb = -0.5*Valfa - sqrt(3)/2*Vbeta;
        Vc = -0.5*Valfa + sqrt(3)/2*Vbeta;
    end


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

        wr = wr + h*(A*wr + B1*Te + B2*Tl(k));
        wrpm = wr*2/Polos*60/(2*pi);

    end

    % Extração do valor atrasado do buffer usando o índice circular
    Iqs_atrasado = Iqs_buffer(buffer_index);
    Ids_atrasado = Ids_buffer(buffer_index);

    % Atualiza os buffers nos índices circulares com os novos valores
    Iqs_buffer(buffer_index) = Iqs;
    Ids_buffer(buffer_index) = Ids;

    % Atualiza o índice circular
    buffer_index = mod(buffer_index, buffer_size) + 1;

    iqs_vetor(k) = Iqs;
    iqs_atrasado_vetor(k) = Iqs_atrasado;
    ids_atrasado_vetor(k) = Ids_atrasado;
    ids_vetor(k) = Ids;
    Te_vetor(k) = Te;
    wrpm_vetor(k) = wrpm;
    wr_vetor(k) = wr;
    Vd_vetor(k) = Vd;
    Vq_vetor(k) = Vq;
    lambda_qr_vetor(k) = Fqr/weles;
    lambda_dr_vetor(k) = Fdr/weles;
    lambda_dr_est_vetor(k) = lambda_dr_est;
    U_w_vetor(k) = U_w;
    U_iq_vetor(k) = UI_iq;
    U_id_vetor(k) = UI_id;
    e_w_vetor(k) = e_w;
    e_iq_vetor(k) = e_iq;
    e_id_vetor(k) = e_id;
end
%% graficos
figure
% Plotando os dados
plot(t,wr_vetor,t,w_ref); % tom de cinza escuro
legend('Wr','Ref')

figure
% Plotando os dados
plot(t,ids_vetor,t,ids_atrasado_vetor); % tom de cinza escuro
legend('Ids', 'Ids atrasado')

figure
% Plotando os dados
plot(t,iqs_vetor, t ,iqs_ref, t, iqs_atrasado_vetor); % tom de cinza escuro
legend('Iqs', 'Ref', 'Iqs Atrasado')

figure
% Plotando os dados
plot(t,Te_vetor,t, Tl); % tom de cinza escuro
legend('Te', 'TL')