///////////////////////////// EXERCÍCIO 1 /////////////////////////////

// Função do algoritmo iterativo de Jacobi
   // Variáveis de entrada:
      // A: matriz invertível A do sistema Ax = b
      // b: vetor b do sistema Ax = b
      // x0: vetor x inicial
      // E: tolerância para a diferença entre as iterações
      // M: número máximo de iterações
      // norm_type: norma utilizada
   // Variáveis de saída:
      // xk: vetor solução encontrado pelo algoritmo
      // final_delta: variação entre a penúltima e a última iteração
      // k: número de iterações realizadas
      // residue: erro em relação à solução exigida
      
function [xk, final_delta, k, residue] = Jacobi(A, b, x0, E, M, norm_type)
    // Dimensão de A
    n = size(A, 1)
    // Matriz triangular inferior de A
    L = tril(A, -1)
    // Gerando a matriz diagonal de A
    D = diag(diag(A))
    // Matriz triangular superior de A
    U = triu(A, 1)
    // Inicializando o número de iterações
    k = 0
    
    // Inicializando o vetor solução como o vetor inicial
    // e o vetor anterior a ele como um vetor grande
    // (para garantir que entrará no loop while)
    xk = x0
    x_previous = ones(n,1) * 1e30
    
    // Enquanto a variação entre as iterações for menor que a tolerância
    // e o número de iterações for menor que o máximo pedido...
    while norm(xk - x_previous, norm_type) > E && k < M
        // O x anterior passa a ser o x atual
        x_previous = xk
        // Calcula Dx
        Dx = -(L+U)*xk + b
        // Calcula o novo x
        xk = Dx ./ diag(D)
        // Aumenta em 1 o número de iterações
        k = k + 1
    end
    
    // Calculando a variação da última iteração
    final_delta = norm(xk - x_previous, norm_type)
    // Calculando o erro em relação à solução requerida
    residue = norm(b - A * xk, norm_type)
    
endfunction


///////////////////////////// EXERCÍCIO 2 /////////////////////////////


// Função do algoritmo iterativo de Gauss-Seidel com inversa
   // Variáveis de entrada:
      // A: matriz invertível A do sistema Ax = b
      // b: vetor b do sistema Ax = b
      // x0: vetor x inicial
      // E: tolerância para a diferença entre as iterações
      // M: número máximo de iterações
      // norm_type: norma utilizada
   // Variáveis de saída:
      // xk: vetor solução encontrado pelo algoritmo
      // final_delta: variação entre a penúltima e a última iteração
      // k: número de iterações realizadas
      // residue: erro em relação à solução exigida
      
function [xk, final_delta, k, residue] = Gauss_Seidel_inv(A, b, x0, E, M, norm_type)
    // Dimensão de A
    n = size(A, 1)
    // Matriz triangular inferior de A
    L = tril(A, -1)
    // Gerando a matriz diagonal de A
    D = diag(diag(A))
    // Matriz triangular superior de A
    U = triu(A, 1)
    // Inicializando o número de iterações
    k = 0
    
    // Inicializando o vetor solução como o vetor inicial
    // e o vetor anterior a ele como um vetor grande
    // (para garantir que entrará no loop while)
    xk = x0
    x_previous = ones(n,1) * 1e30
    
    // Enquanto a variação entre as iterações for menor que a tolerância
    // e o número de iterações for menor que o máximo pedido...
    while norm(xk - x_previous, norm_type) > E && k < M
        // O x anterior passa a ser o x atual
        x_previous = xk
        // Calculando o novo xk utilizando inversas
        xk = -inv(L+D)*U*xk + inv(L+D)*b
        // Aumenta em 1 o número de iterações
        k = k + 1
    end
    
    // Calculando a variação da última iteração
    final_delta = norm(xk - x_previous, norm_type)
    // Calculando o erro em relação à solução requerida
    residue = norm(b - A * xk, norm_type)
    
endfunction


// Função para resolver um sistema linear triangular inferior
   // Variáveis de entrada:
      // L: a matriz triangular inferior do sistema Lx = b
      // b: vetor b do sistema Lx = b
   // Variáveis de saída:
      // x: solução do sistema
      
function [x] = Solve_triangular_system(L, b)
    // Dimensão da L
    n = size(L,1)
    // Inicializando o vetor x
    x = zeros(n,1)
    // Calculando x por substituição
    x(1) = b(1)/L(1,1)
    for row = 2:n
        sub_factor = L(row,1:(row-1))*x(1:(row-1))
        x(row) = (b(row) - sub_factor)/L(row,row)
    end
endfunction


// Função do algoritmo iterativo de Gauss-Seidel sem inversa
   // Variáveis de entrada:
      // A: matriz invertível A do sistema Ax = b
      // b: vetor b do sistema Ax = b
      // x0: vetor x inicial
      // E: tolerância para a diferença entre as iterações
      // M: número máximo de iterações
      // norm_type: norma utilizada
   // Variáveis de saída:
      // xk: vetor solução encontrado pelo algoritmo
      // final_delta: variação entre a penúltima e a última iteração
      // k: número de iterações realizadas
      // residue: erro em relação à solução exigida
      
function [xk, final_delta, k, residue] = Gauss_Seidel(A, b, x0, E, M, norm_type)
    // Dimensão de A
    n = size(A, 1)
    // Matriz triangular inferior de A
    L = tril(A, -1)
    // Gerando a matriz diagonal de A
    D = diag(diag(A))
    // Matriz triangular superior de A
    U = triu(A, 1)
    
    // Inicializando o número de iterações
    k = 0
    
    // Inicializando o vetor solução como o vetor inicial
    // e o vetor anterior a ele como um vetor grande
    // (para garantir que entrará no loop while)
    xk = x0
    x_previous = ones(n,1) * 1e30
    
    // Enquanto a variação entre as iterações for menor que a tolerância
    // e o número de iterações for menor que o máximo pedido...
    while norm(xk - x_previous, norm_type) > E && k < M
        // O x anterior passa a ser o x atual
        x_previous = xk
        // Calculando o novo x utilizando a função de 
        // resolução de sistemas triangulares
        xk = Solve_triangular_system(L+D, -U*xk+b)
        // Aumenta em 1 o número de iterações
        k = k + 1
    end
    
    // Calculando a variação da última iteração
    final_delta = norm(xk - x_previous, norm_type)
    // Calculando o erro em relação à solução requerida
    residue = norm(b - A * xk, norm_type)
    
endfunction


///////////////////////////// EXERCÍCIO 3 /////////////////////////////

disp("---------- Exercício 3 ----------")

disp("----- Teste 1 -----")

// Sistema e vetor inicial
A1 = [1 -4 2; 0 2 4; 6 -1 -2]
b1 = [2; 1; 1]
x0 = [0; 0; 0]

// Testando com a função de Jacobi
[x1_jacobi, fd1_jacobi, k1_jacobi, r1_jacobi] = Jacobi(A1, b1, x0, 10e-4, 50, 2)

disp("Solução com Jacobi:", x1_jacobi)
disp("Variação final com Jacobi:", fd1_jacobi)
disp("Número de iterações com Jacobi:", k1_jacobi)
disp("Resíduo com Jacobi:", r1_jacobi)

// Testando com a função de Gauss-Seidel com inversa
[x1_gsi, fd1_gsi, k1_gsi, r1_gsi] = Gauss_Seidel_inv(A1, b1, x0, 10e-4, 50, 2)

disp("Solução com Gauss-Seidel com inversa:", x1_gsi)
disp("Variação final com Gauss-Seidel com inversa:", fd1_gsi)
disp("Número de iterações com Gauss-Seidel com inversa:", k1_gsi)
disp("Resíduo com Gauss-Seidel com inversa:", r1_gsi)

// Testando com a função de Gauss-Seidel
[x1_gs, fd1_gs, k1_gs, r1_gs] = Gauss_Seidel(A1, b1, x0, 10e-4, 50, 2)

disp("Solução com Gauss-Seidel:", x1_gs)
disp("Variação final com Gauss-Seidel:", fd1_gs)
disp("Número de iterações com Gauss-Seidel:", k1_gs)
disp("Resíduo com Gauss-Seidel:", r1_gs)

disp("----- Teste 2 -----")

// Sistema com as linhas permutadas
A2 = [6 -1 -2; 1 -4 2; 0 2 4]
b2 = [1; 2; 1]

// Testando com a função de Jacobi
[x2_jacobi, fd2_jacobi, k2_jacobi, r2_jacobi] = Jacobi(A2, b2, x0, 10e-4, 50, 2)

disp("Solução com Jacobi:", x2_jacobi)
disp("Variação final com Jacobi:", fd2_jacobi)

disp("Número de iterações com Jacobi:", k2_jacobi)
disp("Resíduo com Jacobi:", r2_jacobi)

// Testando com a função de Gauss-Seidel com inversa
[x2_gsi, fd2_gsi, k2_gsi, r2_gsi] = Gauss_Seidel_inv(A2, b2, x0, 10e-4, 50, 2)

disp("Solução com Gauss-Seidel com inversa:", x2_gsi)
disp("Variação final com Gauss-Seidel com inversa:", fd2_gsi)
disp("Número de iterações com Gauss-Seidel com inversa:", k2_gsi)
disp("Resíduo com Gauss-Seidel com inversa:", r2_gsi)

// Testando com a função de Gauss-Seidel
[x2_gs, fd2_gs, k2_gs, r2_gs] = Gauss_Seidel(A2, b2, x0, 10e-4, 50, 2)

disp("Solução com Gauss-Seidel:", x2_gs)
disp("Variação final com Gauss-Seidel:", fd2_gs)
disp("Número de iterações com Gauss-Seidel:", k2_gs)
disp("Resíduo com Gauss-Seidel:", r2_gs)

///////////////////////////// EXERCÍCIO 4 /////////////////////////////

disp("---------- Exercício 4a ----------")

// Matriz e vetor do sistema
A3 = [2 -1 1; 2 2 2; -1 -1 2]
b3 = [-1; 4; -5]

// Testando com a função de Jacobi
[x3_jacobi, fd3_jacobi, k3_jacobi, r3_jacobi] = Jacobi(A3, b3, x0, 10e-4, 25, 2)

disp("Solução:", x3_jacobi)
disp("Variação final:", fd3_jacobi)
disp("Número de iterações:", k3_jacobi)
disp("Resíduo:", r3_jacobi)

disp("---------- Exercício 4b ----------")

// Testando com a função de Gauss-Seidel com a norma-infinito
[x3_gs, fd3_gs, k3_gs, r3_gs] = Gauss_Seidel(A3, b3, x0, 10e-5, 25, 'inf')

disp("Solução:", x3_gs)
disp("Variação final:", fd3_gs)
disp("Número de iterações:", k3_gs)
disp("Resíduo:", r3_gs)

///////////////////////////// EXERCÍCIO 5 /////////////////////////////

disp("---------- Exercício 5a ----------")

A4 = [1 0 -1; -1/2 1 -1/4; 1 -1/2 1]
b4 = [0.2; -1.425; 2]

// Testando com a função de Gauss-Seidel
[x4_gs, fd4_gs, k4_gs, r4_gs] = Gauss_Seidel(A4, b4, x0, 10e-2, 300, 2)

disp("Solução:", x4_gs)
disp("Variação final:", fd4_gs)
disp("Número de iterações:", k4_gs)
disp("Resíduo:", r4_gs)

disp("---------- Exercício 5b ----------")

A5 = [1 0 -2; -1/2 1 -1/4; 1 -1/2 1]
b5 = [0.2; -1.425; 2]

// Testando com a função de Gauss-Seidel
[x5_gs, fd5_gs, k5_gs, r5_gs] = Gauss_Seidel(A5, b5, x0, 10e-2, 300, 2)

disp("Solução:", x5_gs)
disp("Variação final:", fd5_gs)
disp("Número de iterações:", k5_gs)
disp("Resíduo:", r5_gs)

///////////////////////////// EXERCÍCIO 6 /////////////////////////////

disp("---------- Exercício 6 ----------")

// Função para gerar uma matriz estritamente diagonal dominante
   // Variáveis de entrada:
      // n: dimensão da matriz
   // Variáveis de saída:
      // A: uma matriz nxn estritamente diagonal dominante
      
function [A] = SDD_matrix(n)
    // Gera uma matriz nxn aleatória
    A = rand(n,n)
    // Para cada linha...
    for row = 1:n
        // Calcula a soma de todos os termos fora da diagonal
        row_sum = sum(abs(A(row,:))) - abs(A(row,row))
        // Adiciona essa soma ao elemento da diagonal para
        // garantir que ele será maior que essa soma
        A(row,row) = abs(A(row,row)) + row_sum
    end
endfunction

// Dimensões das matrizes a serem testadas
n = [10, 100, 1000, 2000]
// Matriz com os tempos de execução
timers = zeros(4,2)

// Para cada dimensão...
for i = 1:4
    // Gera uma matriz estritamente diagonal dominante
    A = SDD_matrix(n(i))
    // Gera um vetor b de dimensão coerente
    b = rand(n(i),1)
    // Gera um vetor inicial de zeros de dimensão coerente
    x0 = zeros(n(i),1)
    // Salva o tempo de execução da Gauss-Seidel sem inversa
    tic()
    x = Gauss_Seidel(A, b, x0, 10e-4, 50, 2)
    timers(i,1) = toc()
    // Salva o tempo de execução da Gauss-Seidel com inversa
    tic()
    x = Gauss_Seidel_inv(A, b, x0, 10e-4, 50, 2)
    timers(i,2) = toc()
end

disp(" Sem inv     Com inv ")
disp(timers)















