%% Leitura de dados
%Primeiramente temos de fazer os gráficos das frequências cardíacas ao
%longo do tempo. Os dados fornecidos contêm os registos dos intervalos RR
%dos individuos durante uma hora
table_signal = readtable("C:\Users\guilh\OneDrive\Documentos\MATLAB\ECG\ECG\Colo_rectal\Dados da Sofia Silvestre\Intervalos_RR_1h\Doente01\01_Dia1.txt");
signal_raw = table2array(table_signal);
signal_raw = signal_raw./1000; %(passar para segundos)

%% RECUPERAÇÃO DE DADOS TEMPORAIS
samples = length(signal_raw);
%Para os posicionar no tempo, temos de fazer o seguinte ciclo while
T_base=[0];
for i=[1:samples-1] %Samples-1 porque o último valor não tem de ser somado para obter o tempo do próximo batimento, visto que não há mais nenhum registo a seguir a esse
    ponto_seguinte = T_base(i) + signal_raw(i);
    T_base = [T_base ponto_seguinte];
end

%% ELIMINAÇÃO DE OUTLIERS
%Verificar vários pacientes para perceber se são mesmo outliers
signal_out = filloutliers(signal_raw,'pchip', 'movmean',35);
signal_out = movmean(signal_out, 30)
%% Visualização
figure
subplot(2,1,1)
plot(T_base, signal_raw, 'b', 'LineWidth', 1.5); 
title('Intervalos RR - Raw')
xlabel('Tempo (s)')
ylabel('Intervalo RR (s)')
grid on
subplot(2,1,2)
plot(T_base, signal_out, 'b', 'LineWidth', 1.5); 
title('Intervalos RR - Suavizado')
xlabel('Tempo (s)')
ylabel('Intervalo RR (s)')
grid on


%% 10. EXTRAÇÃO DAS Frequencias? (LF & HF)
%STFT
fs_interp = 4; % 4Hz
T_uniforme = T_base(1) : (1/fs_interp) : T_base(end);
signal_interp = interp1(T_base, signal_out, T_uniforme, "pchip");
window = hamming(256, "periodic") %evita cortes abrubtos nas frequências e criação de artefactos
[s, f, t] = stft(signal_interp, fs_interp,'Window', window,'OverlapLength', 64,'FFTLength', 512);
PSD = abs(s).^2;

% LF/HF indices
LF_idx = (f >= 0.04) & (f < 0.15);
HF_idx = (f >= 0.15) & (f < 0.40);

% Band powers over time
LF_t = trapz(f(LF_idx), PSD(LF_idx,:),1);
HF_t = trapz(f(HF_idx), PSD(HF_idx,:),1);

%Ativações
SNS_t = LF_t ./ HF_t; %Mede a ativação do SNS em relação ao SNP
PNS_t = HF_t;

%Normalização
SNS_norm = SNS_t/mean(SNS_t);
%SNS_norm = SNS_norm-mean(SNS_norm);
PNS_norm = PNS_t/mean(PNS_t);
%PNS_norm = PNS_norm-mean(PNS_norm);
figure, subplot(2,1,1), plot(t, PNS_norm), subplot(2,1,2), plot(t,SNS_norm)
%% Poincaré
%Este método basicamente é um espaço de fases a duas dimensões com um
%formato eliptico amplo e central associado a um comportamento saudável e
%um comportamento eliptico curto e na periieria associado a doentes

RR = signal_out(1:end-1);
RR_1 = signal_out(2:end);
center_y = mean(RR);
center_x = mean(RR_1);

%Parâmeteros do gráfico com todos os dados
CCD = sqrt(center_x^2+center_y^2)%baseline cardiac cycle duration (CCD) comprimento do vetor que aponta para o centro
SD1 = sqrt(1/2*std(diff(RR))^2);
SD2 = sqrt(abs(2*std(RR)^2-1/2*std(diff(RR))^2));

%Poincaré
figure
scatter(signal_raw(1:end-1), signal_out(2:end), 'filled')
xlabel("Intervalos RR")
ylabel("Intervalos RR_+_1")
title("Gráfico Poincaré")
hold on
plot(center_x,center_y, 'r.', "MarkerSize",10)
text(center_x, center_y, '  Centro', 'FontWeight','bold');

% Variação dos Parâmetros com o tempo
CCD_t = [];
SD1_t = []
SD2_t = []
%Precisamos de usar o T_base, que representa o tempo do nosso sinal
%Vamos calculando os mesmos parâmetros para intervalos de 15 segundos (lembrar que a
%data está em milisegundos) e juntar num array para obtermos a variação no
%tempo
w = 15 %Intervalo
T_w = 0:w:max(T_base)+w %este +w existe para garantir que todos os dados são incluídos
for i = 1:length(T_w)-1
    if T_w(i) <= max(T_base) %Não interessa se o valor de T_w(i+1) é maior do que o max(T_base) porque o sistema assume só os valores menores que isso que não é conjunto vazio
        idx = find(T_base >= T_w(i)  & T_base <= T_w(i+1)); %vai buscar os indexs
    else
        break %Se T_w(i) for maior do que T_base, ent dá conjunto vazio, então acaba aqui o ciclo
    end
    rr_15 = signal_out(idx); %encontra os valores do sinal com os indexes respetivos
    %calcula cada parâmetro e adciona ao array
    sd1 = sqrt(1/2*std(diff(rr_15))^2);
    sd2 = sqrt(abs(2*std(rr_15)^2-1/2*std(diff(rr_15))^2));
    ccd = sqrt(mean(rr_15)^2+mean(rr_15(2:end))^2);
    SD1_t = [SD1_t, sd1];
    SD2_t = [SD2_t, sd2];
    CCD_t = [CCD_t, ccd];
end
T_SD = w:w:w*length(SD1_t); %atualiza os dados temporais

%Plotar SD1 e SD2 (estes gráfico servem só para visualizar os dados e ter uma ideia deles, mas não
%permitem chegar a nenhuma conclusão)
figure
subplot(3,1,1)
plot(T_SD, SD1_t)
title("SD1 Over Time")
subplot(3,1,2)
plot(T_SD,SD2_t)
title("SD2 Over Time")
subplot(3,1,3)
plot(T_SD,CCD_t)
title("CCD Over Time") %É normal ser muito parecido ao sinal original visto que ele mede a distância ao centro de todos os pontos do sinal, que varia de acordo com o sinal

%Cálculos do CSI e CPI
CCD_det = detrend(CCD_t);
SD1_det = detrend(SD1_t);
SD2_det = detrend(SD2_t);
D = CCD + CCD_det;
D_flipped = 2.*mean(D)-D;

%Supostamente existe uma constante k que define o peso de cada SD para o
%CSI e CPI
%Como não nos foram sugeridos os valores dos k's, vamos obter a partir da
%normalização:
%- Fazer k a partir das constantes de cada parâmetro é irrealista porque
%não considera a variação dos valores no tempo
%- Fazer k ponto a ponto pode gerar algum ruído devido à alta sensibilidade
%a variações
%- Vamos fazer k aplicado em janelas e considerar esse o k para todos os
%valores da janela
ks = get_k(SD2_det, D_flipped, 5);
kp = get_k(SD1_det, D, 5);


CPI = kp.*(SD1 + SD1_det)+D; %O simpático contribui positivamente para o CCD, por isso é que não está flipped
CSI = ks.*(SD2 + SD2_det)+D_flipped; %O parassimpático contribui negativamente para o CCD, por isso é que está flipped

%Normalização
CPI_norm = CPI/mean(CPI);
CSI_norm = CSI/mean(CSI);

%Detrend
CPI_det = detrend(CPI)
CSI_det = detrend(CSI)


%Plotar CSI e CPI com um filtro a suavizar
figure
subplot(2,1,1)
plot(T_SD, movmean(CPI_det,3), 'b')
title("CPI")

subplot(2,1,2)
plot(T_SD, sgolayfilt(CSI_det,3,11), 'red')
title("CSI")
%%
figure
boxplot([CSI_norm(:), CPI_norm(:)])
xticklabels({'SNS','PNS'})
ylabel('Normalized Power')
title('Distribuição dos índices autonómicos')