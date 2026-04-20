%% Leitura de dados
%Primeiramente temos de fazer os gráficos das frequências cardíacas ao
%longo do tempo. Os dados fornecidos contêm os registos dos intervalos RR
%dos individuos durante uma hora
table_signal = readtable("/Users/mirandaaa/Desktop/Fct- UNI/4º Ano/2º Semestre/ASB/Projeto/ECG/Colo_rectal/Dados da Sofia Silvestre/Intervalos_RR_1h/Doente01/01_Dia1.txt");
signal_raw = table2array(table_signal);

%% RECUPERAÇÃO DE DADOS TEMPORAIS
samples = length(signal_raw);
%Para os posicionar no tempo, temos de fazer o seguinte ciclo while
T_base=[0];
for i=[1:samples-1] %Samples-1 porque o último valor não tem de ser somado para obter o tempo do próximo batimento, visto que não há mais nenhum registo a seguir a esse
    ponto_seguinte = T_base(i) + signal_raw(i);
    T_base = [T_base ponto_seguinte];
end

figure
plot(T_base, signal_raw, 'b', 'LineWidth', 1.5); 
title('Intervalos RR')
xlabel('Tempo (s)')
ylabel('Intervalo RR (s)')
grid on

%% INTERPOLAÇÂO
%A interpolação é necessária porque não temos uma frequência de amostragem
%fixa, logo não conseguimos aplicar nem transformada de fourier nem
%Transformada de Wavelets descreta
%O objetivo será criar pontos para que o sinal fique com intervalos
%constantes de sampling
%new_T = 0:300:sum(signal); se fizesse isto teria que normalizar para poder comparar os sinais, não
%daria, complicava bastante
new_T = 0:500:3600000; %aqui posso atribuir o valor que eu quiser para que todos os sinais tenham o mesmo tamanho em x, convém ser um valoro
rr_intervals = interp1(T_base, signal_raw, new_T, "spline"); %o spline consegue dar uma continuação mais orgânica ao batimento cardíaco, evitando a presneça de muitos outliers
figure, subplot(1,2,1), plot(T_base, signal_raw), title("Sem Interpolação")
subplot(1,2,2), plot(new_T, rr_intervals),title("Com Interpolação")

%% BPM
BPM_s = 60000./signal_raw;
new_BPM = 60000./rr_intervals;
figure, subplot(1,2,1), plot(T_base, BPM_s), title("Sem Interpolação")
subplot(1,2,2), plot(new_T, new_BPM),title("Com Interpolação")

%% SUAVIZAÇÃO
%tem de ser depois da interpolação, porque se for antes da interpolação, os
%intervalos RR ficam contaminados, produzindo dados incorretos a partir
%daqui
%A suavização neste caso só serve para facilitar a visualização do gráfico
%e não para alterar dados antes de estes serem usados em futuras operações
%O objetivo será tentar sempre manter os dados o maiss RAW possíveis
BPM_mean = movmean(BPM_s, 50); %1000 suaviza demaisado o sinal
new_BPM_mean = movmean(new_BPM, 50);
figure,  subplot(1,2,1), plot(T_base, BPM_mean),title("Sem Interpolação"),
subplot(1,2,2), plot(new_T, new_BPM_mean), title("Com Interpolação")

%% BPM MÉDIO
BPM = mean(new_BPM) % que é igual a = mean(BPM), ou seja, o impacto da suavização e interpolação nos dados foi bastante baixo, garantindo a integridade inicial dos dados

%% ESPETRO DE POTÊNCIA
%Permite relacionar com o sistema nervoso (tenho as relações escritas numa
%folha)
fs = 1/0.5; %agora temos uma frequência
HRV = pwelch(rr_intervals, [],[],[], fs);
figure, plot(HRV), title("HRV")

%% HRV
% RMSSD
rmssd = sqrt(mean((diff(intervalo_RR*1000)).^2));
% SDNN
sdnn=sqrt(mean((intervalos_RR-mean(intervalos_RR*1000)).^2));
% pNN50


%% CÁLCULO DA DIFERENÇA (Delta RR) E SEPARAÇÃO (SNS vs PNS)
% Calcular a diferença de tempo entre o batimento atual e o anterior
% diff() encurta o vetor em 1 por isso colocamos um 0 no início para manter o tamanho
delta_rr = [0; diff(signal_raw)];

%FC aumenta, o t entre batimentos diminui (delta_rr negativo)
idx_sns = find(delta_rr < 0);
T_sns = T_base(idx_sns);       % Eixo temporal só com os SNS
rr_sns = signal_raw(idx_sns);    % Valores de RR associados ao SNS

%FC diminui, o t entre batimentos aumenta (delta_rr positivo)
idx_pns = find(delta_rr > 0);
T_pns = T_base(idx_pns);       %eixo temporal só com instantes PNS
rr_pns = signal_raw(idx_pns);    %valores de RR associados ao PNS

figure;
plot(T_base, signal_raw, 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5); hold on;

%plota apenas os do SNS a vermelho
plot(T_sns, rr_sns, 'r.', 'MarkerSize', 12); 
plot(T_pns, rr_pns, 'b.', 'MarkerSize', 12); 
title('Separação SNS vs PNS')
xlabel('Tempo (s)')
ylabel('Intervalo RR (s)')
legend('Sinal Original', 'SNS', 'PNS')

%% INTERPOLAÇÃO SPLINE 
% ao extrair os pontos, os dados T_sns e T_pns estão cheios de missing
% data basicamnete (necessário reamostragem)
%necessario pois pwelch necessita que a distância temporal entre todas as amostras seja igual

fs_interp = 4; % 4 é o standard que encontrei que se utiliza HRV (HRV max = 0.4Hz)
T_uniforme = T_base(1) : (1/fs_interp) : T_base(end); %ENCONTREI ESTA FORMA MAIS SIMPLES DE CRIAR O EIXO TEMPORAL (uma amostra a cada 0.25 segundos)
%começa no primeiro batimento (T_base(1))) e vai até ao ultimo

%onde os pontos estavam(T_sns), qual era o valor deles(rr_sns), e para onde vão os pontos(T_uniforme)
sns_interp = spline(T_sns, rr_sns, T_uniforme);
pns_interp = spline(T_pns, rr_pns, T_uniforme);

%Centrar o sinal no zero (Remover a componente DC (fazer a média)).
%correspnde ao pico de 0Hz no pwelch (assim vamos conseguir ver as HF e as
%LF melhor pois retiramos o pico grande)
%sns_interp = sns_interp - mean(sns_interp);
%pns_interp = pns_interp - mean(pns_interp);


%% VISUALIZAÇÃO 
figure
subplot(2,1,1)
plot(T_uniforme, sns_interp, 'r');
title('Série SNS')
ylabel('Amplitude (s)')
grid on
subplot(2,1,2)
plot(T_uniforme, pns_interp, 'b');
title('Série PNS (Marcadores de Desaceleração) - Pronta para PSD')
xlabel('Tempo (s)')
ylabel('Amplitude (s)')
grid on

%% ANÁLISE ESPETRAL (PWelch)


%% VISUALIZAÇÃO


%% 10. EXTRAÇÃO DAS Frequencias? (LF & HF)
