function k = get_k(SD, D, w)
    k=[];
    n = length(SD)
    for i = 1:w:n %vai andar de w em w // A cada ciclo i vai ser um valor de 5 em 5
        if i == n %se acontecer i ser igual a n, vai ser igual ao k anterior
            k = [k,k_one];
        elseif i+w-1 <= n
            SD_std = std(SD(i:i+w-1)); %+w-1 para garantir que temos w valores
            D_std = std(D(i:i+w-1));
            k_one = SD_std/D_std;
            k = [k,repmat(k_one,1,w)]; %adiciona o valor de k w vezes ao vetor para acertar os eixos
        else
            SD_std = std(SD(i:n));
            D_std = std(D(i:n));
            k_one = SD_std/D_std;
            k = [k,repmat(k_one,1,n-i+1)];
        end
    end
end
