%% Leitura de dados
%Primeiramente temos de fazer os gráficos das frequências cardíacas ao
%longo do tempo. Os dados fornecidos contêm os registos dos intervalos RR
%dos individuos durante uma hora
table_signal = readtable("/Users/mirandaaa/Desktop/Fct- UNI/4º Ano/2º Semestre/ASB/Projeto/ECG/Colo_rectal/Dados da Sofia Silvestre/Intervalos_RR_1h/Doente01/01_Dia-1.txt");
signal_raw = table2array(table_signal);

%% ELIMINAÇÃO DE OUTLIERS
%Verificar vários pacientes para perceber se são mesmo outliers
%for i = 1:length(signal_raw)
    %if abs(x(i+1)-x(i)) >= 700 
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

%% CÁLCULO DA DIFERENÇA (Delta RR) E SEPARAÇÃO (SNS vs SNP)
% Calcular a diferença de tempo entre o batimento atual e o anterior
% diff() encurta o vetor em 1 por isso colocamos um 0 no início para manter o tamanho
delta_rr = [0; diff(signal_raw)];

%FC aumenta, o t entre batimentos diminui (delta_rr negativo)
idx_sns = find(delta_rr < 0);
T_sns = T_base(idx_sns);       % Eixo temporal só com os SNS
rr_sns = signal_raw(idx_sns);    % Valores de RR associados ao SNS

%FC diminui, o t entre batimentos aumenta (delta_rr positivo)
idx_pns = find(delta_rr > 0);
T_pns = T_base(idx_pns);       %eixo temporal só com instantes SNP
rr_pns = signal_raw(idx_pns);    %valores de RR associados ao SNP

figure;
plot(T_base, signal_raw, 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5); hold on;

%plota apenas os do SNS a vermelho
plot(T_sns, rr_sns, 'r.', 'MarkerSize', 12); 
plot(T_pns, rr_pns, 'b.', 'MarkerSize', 12); 
title('Separação SNS vs PNS')
xlabel('Tempo (s)')
ylabel('Intervalo RR (s)')
legend('Sinal Original', 'SNS', 'PNS')

%% CÁLCULO LINEAR
% O SDNN é simplesmente o desvio padrão da série de intervalos RR normais
SDNN_valor = std(signal_raw); 
fprintf('Valor do SDNN (Domínio do Tempo): %f ms\n\n', SDNN_valor);
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
%Usamos o detrend para remover o pico de 0Hz e tendências muito lentas
sns_interp = detrend(sns_interp);
pns_interp = detrend(pns_interp);


%% VISUALIZAÇÃO 
figure
subplot(2,1,1)
plot(T_uniforme, sns_interp, 'r');
title('SNS')
ylabel('Amplitude')
grid on
subplot(2,1,2)
plot(T_uniforme, pns_interp, 'b');
title('SNP')
xlabel('Tempo')
ylabel('Amplitude')
grid on


%% PWELCH E POTÊNCIAS
window_length = 300 * fs_interp; %janelas de 5 minutos (300 seg)
overlap = round(window_length / 2);

[pxx_sns, f_sns] = pwelch(sns_interp, window_length, overlap, [], fs_interp);
[pxx_pns, f_pns] = pwelch(pns_interp, window_length, overlap, [], fs_interp);

%% VISUALIZAÇÃO FINAL
figure;
%subplot do SNS
subplot(2,1,1);
plot(f_sns, pxx_sns, 'r', 'LineWidth', 1.5);
xlim([0.03 0.5]); % Zoom para focar no LF/HF
title('Potências SNS (Aceleração)');
grid on;

%subplot do SNP
subplot(2,1,2);
plot(f_pns, pxx_pns, 'b', 'LineWidth', 1.5);
xlim([0.03 0.5]);
title('Potências SNP (Desaceleração)');
xlabel('Frequência (Hz)');
grid on;
%% EXTRAÇÃO DE LF E HF 
LF_range = [0.04 0.15];
HF_range = [0.15 0.40];

LF_SNS = bandpower(pxx_sns, f_sns, LF_range, 'psd');
HF_SNS = bandpower(pxx_sns, f_sns, HF_range, 'psd');

LF_PNS = bandpower(pxx_pns, f_pns, LF_range, 'psd');
HF_PNS = bandpower(pxx_pns, f_pns, HF_range, 'psd');

fprintf('SNS:\n LF: %f | HF: %f\n', LF_SNS, HF_SNS);
fprintf('SNP:\n LF: %f | HF: %f\n', LF_PNS, HF_PNS);