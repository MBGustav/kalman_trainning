function [H] = general_H(U, list_W)
  % U       - valores de entrada da rede neural
  % list_W  - lista de matriz de pesos
  % k_total: total de camadas a serem dessenvolvidas
  
  %calcula H1:
  #W1 = list_W{1};
  #v1 = W1 * U
  #H = {sigmoid_f(v1)};
  
  k_total = size(list_W, 2);
  for k=1:k_total
    Wk = list_W{k};
    if k==1  
      y = sigmoid_f(Wk * U);
    else
      y = sigmoid_f(Wk * prev_y);
    end
    H{k} = y;
    
    prev_y = y;
  endfor
  
end
