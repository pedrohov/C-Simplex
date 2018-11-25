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
    Matrix A;
    Matrix b;
    Matrix custo;
    Matrix base;
    Matrix naoBase;
    Matrix B;
    Matrix solucao;
};

// Sub-rotinas:
int simplex(Model modelo) {
    /* 
     * Executa o simplex e atualiza o campo 'solucao' do modelo.
     * Retorna -1 se o modelo nao possuir solucao;
     * Retorna -2 se a solucao for ilimitada;
     * Retorna -3 se houverem multiplas solucoes;
     * Retorna  1 se houver uma unica solucao otima.
     * [DEBUG] matImprime(Matrix m) printa qualquer matriz/vetor no console.
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

            // (REQUISITO 09) - Multiplos Otimos:
            if(existeNaoBase0(modelo, custoReduzido)){
                printf("[Simplex] Multiplos Otimos.\nIteracoes = %d\nCusto otimo = %lf\n", iteracoes, objetivo);
                return -3;
            }
            // (REQUISITO 09) - Solucao Otima Unica:
            else {
                printf("[Simplex] Solucao Otima Unica.\nIteracoes = %d\nCusto otimo = %lf\n", iteracoes, objetivo);
                return 1;
            }
        }

        matLibera(custoBaseT);
        matLibera(custoReduzido);

        // [DEBUG] Indice da variavel a entrar na base:
        // printf("Entra na base %d\n", indiceMenor);

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
    double valor;
    fscanf(parq, "%d", &i);
    fscanf(parq, "%d", &j);

    // Cria vetor b:
    Matrix b = matCria(i, 1);
    for(k = 0; k < i; k++) {
        fscanf(parq, "%lf", &valor);
        matPut(b, k, 0, valor);
    }

    // Cria vetor de custo:
    Matrix c = matCria(j, 1);
    for(k = 0; k < j; k++) {
        fscanf(parq, "%lf", &valor);
        matPut(c, k, 0, valor);
    }

    // Cria vetor base:
    Matrix base = matCria(i, 1);
    for(k = 0; k < i; k++) {
        fscanf(parq, "%lf", &valor);
        matPut(base, k, 0, valor - 1);
    }

    // Cria vetor nao-base:
    Matrix naoBase = matCria((j - i), 1);
    for(k = 0; k < (j - i); k++) {
        fscanf(parq, "%lf", &valor);
        matPut(naoBase, k, 0, valor - 1);
    }

    // Cria matriz A:
    Matrix A;
    A = matCria(i, j);

    // Le e salva cada elemento na matriz:
    for(i = 0; i < matNlinhas(A); i++) {
        for(j = 0; j < matNcolunas(A); j++) {
            fscanf(parq, "%lf", &valor);
            matPut(A, i, j, valor);
        }
    }

    // Empacota os dados em um novo modelo:
    Model modelo = criaModelo(A, b, c, base, naoBase);

    fclose(parq);
    return modelo;
}

Model criaModelo(Matrix A, Matrix b, Matrix c, Matrix base, Matrix naoBase) {
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
    novo -> solucao = NULL;

    // Cria a matriz basica:
    Matrix B = matCria(matNlinhas(A), matNlinhas(A));
    novo -> B = B;

    return novo;
}