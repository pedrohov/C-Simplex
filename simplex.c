/*
 * Pesquisa Operacional - 2o sem/2018
 * IFMG Campus Formiga
 * Pedro Henrique Oliveira Veloso (0002346)
 * Savio Silva Simoes (0011942)
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "simplex.h"

// Estrutura:
struct Tmodel {
    Matrix A;           // Matriz A.
    Matrix b;           // Vetor de igualdade das restricoes.
    Matrix custo;       // Vetor de custo (min).
    Matrix base;        // Vetor com as variaveis na base.
    Matrix naoBase;     // Vetor com as variaveis nao-base.
    Matrix artificiais; // Vetor com as variaveis artificiais.
    Matrix B;           // Matriz das colunas de cada variavel base.
    Matrix solucao;     // Vetor solucao das variaveis presentes na base.
};

// Sub-rotinas:
Model carregaModelo(char arquivo[]) {
     // Abre o arquivo para leitura:
    FILE *parq;
    parq = fopen(arquivo, "r");

    // Caso o arquivo nao exista, retorna NULL:
    if(parq == NULL)
    {
        printf("Arquivo nao encontrado.\n");
        return NULL;
    } 

    // Determinar quantidade de linhas e colunas;
    int i, j, k;
    double valor, bigM;
    fscanf(parq, "%d", &i);
    fscanf(parq, "%d", &j);

    // Cria vetor b:
    Matrix b = matCria(i, 1);
    for(k = 0; k < i; k++) {
        fscanf(parq, "%lf", &valor);
        matPut(b, k, 0, valor);
    }

    // Cria vetor de custo (Add 'i' variaveis artificiais):
    Matrix c = matCria(j + i, 1);
    for(k = 0; k < j; k++) {
        fscanf(parq, "%lf", &valor);
        matPut(c, k, 0, valor);
    }

    // Cria vetor base, nao-base e variaveis artificiais:
    Matrix base = matCria(i, 1);
    Matrix naoBase = matCria(j, 1);
    Matrix artificiais = matCria(i, 1);
    bigM = 99999;
    for(k = 0; k < (i + j); k++) {
        if(k < j)
            matPut(naoBase, k, 0, k);
        else {
            matPut(base, (i + j) - k - 1, 0, k);
            matPut(artificiais, (i + j) - k - 1, 0, k);
            matPut(c, k, 0, bigM);
        }
    }

    // Cria matriz A:
    Matrix A;
    A = matCria(i, (j + i));

    // Le e salva cada elemento na matriz:
    for(i = 0; i < matNlinhas(A); i++) {
        for(k = 0; k < (matNcolunas(A) - matNlinhas(A)); k++) {
            fscanf(parq, "%lf", &valor);
            matPut(A, i, k, valor);
        }
        // Forma colunas extras para variaveis artificiais:
        matPut(A, i, (i + j), 1);
    }

    // Empacota os dados em um novo modelo:
    Model modelo = criaModelo(A, b, c, base, naoBase, artificiais);

    fclose(parq);
    return modelo;
}

Model criaModelo(Matrix A, Matrix b, Matrix c, Matrix base, Matrix naoBase, Matrix artificiais) {
    // Prerequisitos (matrizes e vetores validos):
    if((A == NULL) || (b == NULL) || (c == NULL) || (base == NULL) || (naoBase == NULL))
        return NULL;

    // Aloca memoria para o modelo:
    Model novo = (Model) malloc(sizeof(struct Tmodel));

    // Atribui valores:
    novo -> A = A;
    novo -> b = b;
    novo -> custo = c;
    novo -> base = base;
    novo -> naoBase = naoBase;
    novo -> artificiais = artificiais;
    novo -> solucao = NULL;

    // Cria a matriz basica:
    Matrix B = matCria(matNlinhas(A), matNlinhas(A));
    novo -> B = B;

    return novo;
}

double calcObjetivo(Model modelo) {
    // Calcula o valor da funcao objetivo do modelo:
    double objetivo = 0;
    int i;
    for(i = 0; i < matNlinhas(modelo -> base); i++) {
        int baseI = matGet(modelo -> base, i, 0);
        objetivo += matGet(modelo -> custo, baseI, 0) * matGet(modelo -> solucao, i, 0);
    }

    return objetivo;
}

int simplex(Model modelo) {
    /* 
     * Executa o simplex e atualiza o campo 'solucao' do modelo.
     * Retorna -1 se o modelo nao possuir solucao;
     * Retorna -2 se a solucao for ilimitada;
     * Retorna -3 se houverem multiplas solucoes;
     * Retorna >= 0 se houver uma unica solucao otima. Retorna numero de iteracoes.
     * [DEBUG] matImprime(Matrix m) mostra qualquer matriz/vetor no console.
     * [DEBUG] imprimeModelo(Model m) mostra todos os dados do simplex na iteracao atual.
     */

    // Prerequisitos (modelo nao nulo):
    if(modelo == NULL)
        return -1;

    int m = matNlinhas(modelo -> A);
    int n = matNcolunas(modelo -> A);
    int iteracoes = 0;

    while(1) {
        // Calcula a B^-1:
        geraMatrizBase(modelo);
        Matrix invB = matInversa(modelo -> B);

        // (REQUISITO 03) Calcular a SBF inicial (invB x b):
        modelo -> solucao = matProdutoMatricial(invB, modelo -> b);

        // Calcula o valor da funcao objetivo:
        double objetivo = calcObjetivo(modelo);

        // Monta o vetor de custo basico:
        Matrix custoBase = matCria(matNlinhas(modelo -> base), 1);
        int i;
        for(i = 0; i < m; i++) {
            int indice = matGet(modelo -> base, i, 0);
            double valor = matGet(modelo -> custo, indice, 0);
            matPut(custoBase, i, 0, valor);
        }

        // Determina qual variavel entra na base:
        Matrix custoBaseT = matTransposta(custoBase);
        Matrix custoReduzido = matCria(1, n);
        double menorCusto = INFINITY;
        int indiceMenor = -1;
        for(i = 0; i < matNlinhas(modelo -> naoBase); i++) {

            int indice = matGet(modelo -> naoBase, i, 0);

            // Custo do indice na funcao objetivo:
            double custo = matGet(modelo -> custo, indice, 0);

            // Coluna do indice nao basico da matriz A:
            Matrix Aj = matGetColuna(modelo -> A, indice);

            // (REQUISITO 05) Direcao:
            Matrix direcao = matProdutoMatricial(invB, Aj);

            // resOp1 = custoBase x invB:
            Matrix resOp1 = matProdutoMatricial(custoBaseT, invB);

            // resOp2 = custoBase x invB x Aj:
            Matrix resOp2 = matProdutoMatricial(resOp1, Aj);
            
            // (REQUISITO 04) Custo reduzido:
            custo = custo - matGet(resOp2, 0, 0);
            matPut(custoReduzido, 0, indice, custo);
            
            // Trata valores de custo muito pequenos e negativos (i.e: -1e-16).
            // (REQUISITO 09) Degeneracao - Garantir que custos proximos de 0 sejam +0.
            // Arredonda o valor para 5 casas decimais.
            double fac = pow(10, 5);
            custo = round(custo * fac) / fac;

            // (REQUISITO 09) Degeneracao - Custos reduzidos iguais a zero nao entram na base.
            // Se o custo for negativo:
            if(custo < 0) {
                if(custo < menorCusto) {
                    // Atualiza variavel candidata para entrar na base:
                    indiceMenor = indice;
                    menorCusto = custo;
                }
            }
            
            // [DEBUG] (REQUISITO 05) Mostra direcao encontrada:
            // printf("Direcao factivel %d, custo reduzido %lf\n", indice, custo);
            // matImprime(direcao);

            // Libera matrizes utilizadas nas operacoes:
            matLibera(Aj);
            matLibera(resOp1);
            matLibera(resOp2);
            matLibera(direcao);
        }

        // [DEBUG] Mostra o custo reduzido encontrado:
        // printf("Custo reduzido:\n");
        // matImprime(custoReduzido);

        // Nao houve nenhum indice com custo reduzido negativo.
        // Solucao otima encontrada:
        if(indiceMenor == -1) {
            objetivo = calcObjetivo(modelo);

            // (REQUISITO 09) - Solucao inexistente (Variaveis artificiais presentes na solucao):
            if(artificialNaBase(modelo)){
                printf("[Simplex] Solucao inexistente.\n");
                return -1;
            }
            // (REQUISITO 09) - Multiplos Otimos:
            else if(existeNaoBase0(modelo, custoReduzido)){
                printf("[Simplex] Multiplos Otimos.\nIteracoes = %d\nCusto otimo = %lf\n", iteracoes, objetivo);
                return -3;
            }
            // (REQUISITO 09) - Solucao Otima Unica:
            else {
                printf("[Simplex] Solucao Otima Unica.\nIteracoes = %d\nCusto otimo = %lf\n", iteracoes, objetivo);
                return iteracoes;
            }
        }

        matLibera(custoBaseT);
        matLibera(custoReduzido);

        // [DEBUG] Indice da variavel a entrar na base:
        // printf("Entra na base %d, custo = %lf\n", indiceMenor, menorCusto);

        // (REQUISITO 09) - Solucao Ilimitada:
        Matrix Aj = matGetColuna(modelo -> A, indiceMenor);
        Matrix u = matProdutoMatricial(invB, Aj);
        if(!existePositivo(u)) {
            printf("[Simplex] Solucao Ilimitada.\nCusto otimo = -Infinito\n");
            return -2;
        }
        matLibera(Aj);
        matLibera(u);

        // (REQUISITO 06) Determina o valor de theta:
        double theta = INFINITY;
        int indice = -1;

        for(i = 0; i < m; i++) {
            if(matGet(u, i, 0) > 0) {

                double razao = matGet(modelo -> solucao, i, 0) / matGet(u, i, 0);

                if(razao < theta) {
                    theta = razao;
                    indice = matGet(modelo -> base, i, 0);
                }
            }
        }

        // [DEBUG] Indice da variavel a sair da base:
        // printf("Variavel sai da base: %d, theta = %lf\n", indice, theta);

        // (REQUISITO 07) Atualiza a base:
        for(i = 0; i < m; i++) {
            if(matGet(modelo -> base, i, 0) == indice) {
                matPut(modelo -> solucao, i, 0, theta);
                matPut(modelo -> base, i, 0, indiceMenor);
            }
        }

        // (REQUISITO 07) Atualiza o conjunto de variaveis nao-basicas:
        for(i = 0; i < (n - m); i++) {
            if(matGet(modelo -> naoBase, i, 0) == indiceMenor)
                matPut(modelo -> naoBase, i, 0, indice);
        }

        // Libera a matriz inversa utilizada para a iteracao atual:
        matLibera(invB);

        iteracoes++;
    }

    return -1;
}

int artificialNaBase(Model modelo) {
    // (REQUISITO 09) Checa se existe uma variavel artificial na base:
    if((modelo == NULL) || (modelo -> artificiais == NULL))
        return 0;

    int i, j;
    for(i = 0; i < matNlinhas(modelo -> artificiais); i++) {
        int artificial = matGet(modelo -> artificiais, i, 0);
        for(j = 0; j < matNlinhas(modelo -> base); j++) {
            if(artificial == matGet(modelo -> base, j, 0))
                return 1;
        }
    }

    return 0;
}

int existePositivo(Matrix u) {
    // (REQUISITO 09) Checa se existe um valor positivo na linha 'u':
    if(u == NULL)
        return 0;

    int i;
    for(i = 0; i < matNlinhas(u); i++) {
        if(matGet(u, i, 0) > 0)
            return 1;
    }
    return 0;
}

int existeNaoBase0(Model modelo, Matrix custoReduzido) {
    // (REQUISITO 09) Checa se existe uma variavel nao base com custo reduzido zero:
    if((modelo == NULL) || (custoReduzido == NULL))
        return 0;

    int i, j;
    for(i = 0; i < matNcolunas(custoReduzido); i++) {
        double custo = matGet(custoReduzido, 0, i);

        // Trata valores muito pequenos de custo (i.e: 1e-16).
        // Arredonda o valor para 5 casas decimais:
        double fac = pow(10, 5);
        custo = round(custo * fac) / fac;

        // Encontrou um custo reduzido igual a zero:
        if(custo == 0) {
            // Olha se o indice 'i' existe na base:
            int existe = 0;
            for(j = 0; j < matNlinhas(modelo -> base); j++){
                if(matGet(modelo -> base, j, 0) == i) {
                    existe = 1;
                    break;
                }
            }
            if(!existe)
                return 1;
        }
    }

    return 0;
}

void imprimeModelo(Model modelo) {
    // [DEBUG] Exibe  o estado atual do modelo:
    printf("\nMatriz A:\n");
    matImprime(modelo -> A);
    printf("\nVetor b:\n");
    matImprime(modelo -> b);
    printf("\nVetor custo:\n");
    matImprime(modelo -> custo);
    printf("\nVetor base (indices contam a partir de 0):\n");
    matImprime(modelo -> base);
    printf("\nVetor nao-base (indices contam a partir de 0):\n");
    matImprime(modelo -> naoBase);
    printf("\nMatriz basica:\n");
    matImprime(modelo -> B);

    if(modelo -> solucao != NULL) {
        printf("\nSolucao:\n");
        matImprime(modelo -> solucao);
    } else
        printf("\n");

    return;
}

void geraMatrizBase(Model modelo) {
    // Copia as colunas da base 'modelo -> base' para a matriz B:
    if(modelo == NULL)
        return;

    int i, j;
    double valor;
    for(i = 0; i < matNlinhas(modelo -> A); i++) {
        int base = matGet(modelo -> base, i, 0);

        for(j = 0; j < matNlinhas(modelo -> A); j++) {
            valor = matGet(modelo -> A, j, base);
            matPut(modelo -> B, j, i, valor);
        }
    }
    return;
}

void liberaModelo(Model modelo) {
    // Libera a memoria utilizada pelo modelo:
    if(modelo == NULL)
        return;

    matLibera(modelo -> A);
    matLibera(modelo -> b);
    matLibera(modelo -> custo);
    matLibera(modelo -> base);
    matLibera(modelo -> naoBase);
    matLibera(modelo -> B);
    matLibera(modelo -> solucao);

    return;
}

void outputModelo(Model modelo, int iteracoes, char entrada[], char saida[]) {
    if((modelo == NULL) || (modelo -> solucao == NULL))
        return;

    FILE *arq = fopen(saida, "w");

    if(iteracoes == -1)
        fprintf(arq, "Modelo: %s\nSolucao inexistente\n", entrada);
    else if(iteracoes == -2)
        fprintf(arq, "Modelo: %s\nSolucao ilimitada\n", entrada);
    else if(iteracoes == -3)
        fprintf(arq, "Modelo: %s\nMultiplos otimos\n", entrada);
    else
        fprintf(arq, "Modelo: %s\nSolucao Unica Otima encontrada\n\nIteracoes: %d", entrada, iteracoes);

    fprintf(arq, "\nOtimo: %.2lf\n\n", calcObjetivo(modelo));

    // Exibe o valor das variaveis:
    int i, j, k;
    for(i = 0; i < (matNlinhas(modelo -> custo) - matNlinhas(modelo -> A)); i++) {
        double valor = 0;

        // Procura se a variavel 'i' esta na base:
        for(j = 0; j < matNlinhas(modelo -> base); j++) {
            int base = matGet(modelo -> base, j, 0);
            if(i == base) {
                // Pega o valor da base na solucao:
                valor = matGet(modelo -> solucao, j, 0);
            }
        }

        fprintf(arq, "x[%d] = %.2lf\n", (i + 1), valor);
    }

    fclose(arq);
    return;
}