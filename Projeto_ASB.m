%% LER DADOS
%Primeiramente temos de fazer os gráficos das frequências cardíacas ao
%longo do tempo. Os dados fornecidos contêm os registos dos intervalos RR
%dos individuos durante uma hora
table_signal = readtable("C:\Users\guilh\OneDrive\Documentos\MATLAB\ECG\ECG\Colo_rectal\Dados da Sofia Silvestre\Intervalos_RR_1h\Doente01\01_Dia1.txt");
signal = table2array(table_signal);

%% RECUPERAÇÃO DE DADOS TEMPORAIS
samples = length(signal);
%Para os posicionar no tempo, temos de fazer o seguinte ciclo while
T_base=[0];
for i=[1:samples-1] %Samples-1 porque o último valor não tem de ser somado para obter o tempo do próximo batimento, visto que não há mais nenhum registo a seguir a esse
    ponto_seguinte = T_base(i) + signal(i);
    T_base = [T_base ponto_seguinte];
end

%% INTERPOLAÇÂO
%A interpolação é necessária porque não temos uma frequência de amostragem
%fixa, logo não conseguimos aplicar nem transformada de fourier nem
%Transformada de Wavelets descreta
%O objetivo será criar pontos para que o sinal fique com intervalos
%constantes de sampling
%new_T = 0:300:sum(signal); se fizesse isto teria que normalizar para poder comparar os sinais, não
%daria, complicava bastante
new_T = 0:200:3600000; %aqui posso atribuir o valor que eu quiser para que todos os sinais tenham o mesmo tamanho em x
rr_intervals = interp1(T_base, signal, new_T, "spline"); %o spline consegue dar uma continuação mais orgânica ao batimento cardíaco, evitando a presneça de muitos outliers
figure, subplot(1,2,1), plot(T_base, signal), title("Sem Interpolação")
subplot(1,2,2), plot(new_T, rr_intervals),title("Com Interpolação")


%% BPM
BPM_s = 60000./signal;
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
BPM_mean = movmean(BPM_s, 10); %1000 suaviza demaisado o sinal
new_BPM_mean = movmean(new_BPM, 10);
figure, subplot(1,2,1), plot(T_base, BPM_mean),title("Sem Interpolação"),
subplot(1,2,2), plot(new_T, new_BPM_mean), title("Com Interpolação")

%% BPM MÉDIO
BPM = mean(new_BPM) % que é igual a = mean(BPM), ou seja, o impacto da suavização e interpolação nos dados foi bastante baixo, garantindo a integridade inicial dos dados

%% HRV A TRANSFORMADA DE FOURIER ESTÁ A FUNCIONAR MAL!!!!! TALVEZ FAZER A SHORT TIME FOURIER TRANSFORM?????
%Permite relacionar com o sistema nervoso (tenho as relações escritas numa
%folha)
fs = 1/0.2; %agora temos uma frequência
HRV = pwelch(rr_intervals, [],[],[], fs);
figure, plot(HRV), title("HRV")


