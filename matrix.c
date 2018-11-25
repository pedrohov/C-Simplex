/*
 * Pesquisa Operacional - 2o sem/2018
 * IFMG Campus Formiga
 * Pedro Henrique Oliveira Veloso (0002346)
 * Savio Silva Simoes (0011942)
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"

// Estrutura:
struct Tmatrix {
    int lin;    // Quantidade de linhas.
    int col;    // Quantidade de colunas.
    double **p;  // Vetor de ponteiros. Apontam para o inicio de cada linha.
};

// Sub-rotinas:
Matrix matCria(int lin, int col)
{
    // Prerequisitos (lin e col devem ser numeros positivos nao nulos):
    if((lin <= 0) || (col <= 0))
        return NULL;

    // Aloca memoria para uma matriz e define seu tamanho:
    Matrix nova = (Matrix)malloc(sizeof(struct Tmatrix));
    nova -> col = col;
    nova -> lin = lin;
    nova -> p = (double**)malloc(sizeof(double*)*lin);

    // Aloca o numero de colunas necessario p/ cada linha:
    int i, j;
    for(i = 0; i < lin; i++)
        nova -> p[i] =(double*)malloc(sizeof(double)*col);

    // Inicializa a matriz com 0's:
    for(i = 0; i < lin; i++)
        for(j = 0; j < col; j++)
            nova -> p[i][j] = 0;

    return nova;
}

Matrix matCopia(Matrix mat)
{
    // Prerequisitos (Matriz nao nula):
    if(mat == NULL)
        return NULL;

    // Cria uma nova matriz:
    Matrix copia = matCria(mat -> lin, mat -> col);
    int i, j;
    copia -> lin = copia -> lin;
    copia -> col = copia -> col;

    // Copia cada elemento de 'mat':
    for(i = 0; i < mat -> lin; i++)
        for(j = 0; j < mat -> col; j++)
            copia -> p[i][j] = mat -> p[i][j];

    return copia;
}

Matrix matCarrega(char nomeArquivo[])
{
    // Abre o arquivo para leitura:
    FILE *parq;
    parq = fopen(nomeArquivo, "r");

    // Caso o arquivo nao exista, retorna NULL:
    if(parq == NULL)
    {
        printf("Arquivo nao encontrado.\n");
        return NULL;
    } 

    // Determinar quantidade de linhas e colunas;
    int i, j;
    fscanf(parq, "%d", &i);
    fscanf(parq, "%d", &j);

    // Cria matriz:
    Matrix nova;
    nova = matCria(i, j);

    // Le e salva cada elemento na matriz:
    for(i = 0; i < nova -> lin; i++)
        for(j = 0; j < nova -> col; j++)
            fscanf(parq, "%lf", &nova->p[i][j]);

    return nova;
}

Matrix matIdentidade (int n)
{
    // Prerequisitos (n > 0):
    if(n <= 0)
        return NULL;

    // Cria nova matriz:
    Matrix identidade = matCria(n, n);

    // Coloca 1's na linha principal:
    int i;
    for(i = 0; i < n; i++)
        identidade -> p[i][i] = 1;

    return identidade;
}

Matrix matVetMedia (Matrix mat)
{
    // Prerequisitos (Matriz nao nula):
    if(mat == NULL)
        return NULL;

    // Cria uma matriz (n, 1):
    Matrix media = matCria(mat -> col, 1);

    // Armazena a media de cada variavel nas linhas da matriz:
    int i, j;
    for(i = 0; i < mat -> col; i++)
    {
        double soma = 0;
        for(j = 0; j < mat -> lin; j++)
            soma += mat -> p[j][i];

        soma = soma / mat -> lin;
        media -> p[i][0] = soma;
    }
    
    return media; 
}

Matrix matCovariancia (Matrix mat)
{
    // Prerequisitos (Matriz nao nula):
    if(mat == NULL)
        return NULL;

    // Cria a matriz covariancia:
    Matrix cov = matCria(mat -> col, mat -> col);

    // Cria um vetor com as medias de cada coluna:
    double *media = (double*)malloc(sizeof(double) * mat -> col);
    int i, j, k;
    for(j = 0; j < mat -> col; j++)
    {
        // Soma todas as linhas da coluna:
        media[j] = 0;
        for(i = 0; i < mat -> lin; i++)
            media[j] += mat -> p[i][j];

        // Divide pelo numero de linhas:
        media[j] = media[j] / mat -> lin;
    }

    // Calcula a variancia e armazena na diagonal principal:
    for(i = 0; i < mat -> col; i++)
    {
        // Somatorio de (x1i - x1~)(x2i - x2~):
        for(j = 0; j < mat -> col; j++)
        {
            double covar = 0;
            for(k = 0; k < mat -> lin; k++)
                covar = covar + (mat -> p[k][i] - media[i]) * (mat -> p[k][j] - media[j]);

            covar = covar / (mat -> lin - 1);

            // Armazena a variancia na matriz 'cov': 
            cov -> p[i][j] = covar;
        }
    }

    // Libera vetor media:
    free(media);
    return cov;     
}

Matrix matTransposta(Matrix mat)
{
    // Prerequisitos (Matriz nao nula):
    if(mat == NULL)
        return NULL;

    // Cria matriz:
    Matrix transposta = matCria(mat -> col, mat -> lin);

    // Inverte linhas e colunas da original:
    int i, j;
    for(i = 0; i < transposta -> lin; i++)
        for(j = 0; j < transposta -> col; j++)
            transposta -> p[i][j] = mat -> p[j][i];

    return transposta;
}

Matrix matOposta (Matrix mat)
{
    // Prerequisitos (Matriz nao nula):
    if(mat == NULL)
        return NULL;

    // Cria uma copia de 'mat':
    Matrix oposta = matCopia(mat);

    // Multiplica cada elemento por -1:
    int i, j;
    for(i = 0; i < mat -> lin; i++)
        for(j = 0; j < mat -> col; j++)
            oposta -> p[i][j] = -oposta -> p[i][j];
        
    return oposta;
}

Matrix matSoma (Matrix mat1, Matrix mat2)
{
    // Prerequisitos (Matrizes nao nulas, numero de linhas e colunas igual):
    if((mat1 == NULL) || (mat2 == NULL) || (mat1 -> lin != mat2 -> lin) || (mat1 -> col != mat2 -> col))
        return NULL;

    // Soma cada elemento das matrizes:
    int i, j;
    Matrix soma = matCria(mat1 -> lin, mat1 -> col);
    for(i = 0; i < mat1 -> lin; i++)
        for(j = 0; j < mat1 -> col; j++)
            soma -> p[i][j] = mat1 -> p[i][j] + mat2 -> p[i][j];

    return soma;
}

Matrix matSubtrai (Matrix mat1, Matrix mat2)
{
    // Prerequisitos (Matrizes nao nulas, numero de linhas e colunas igual):
    if((mat1 == NULL) || (mat2 == NULL) || (mat1 -> lin != mat2 -> lin) || (mat1 -> col != mat2 -> col))
        return NULL;

    // Subtrai cada elemento das matrizes:
    int i, j;
    Matrix subtracao = matCria(mat1->lin, mat1->col);
    for(i = 0; i < mat1 -> lin; i++)
        for(j = 0; j< mat1 -> col; j++)
            subtracao -> p[i][j] = mat1 -> p[i][j] - mat2 -> p[i][j];

    return subtracao;
}

Matrix matProdutoMatricial (Matrix mat1, Matrix mat2)
{
    // Prerequisitos (Matrizes nao nulas, n.o de colunas de 'mat1' igual ao n.o de linhas de 'mat2'):
    if((mat1 == NULL) || (mat2 == NULL) || (mat1 -> col != mat2 -> lin))
        return NULL;

    // Cria nova matriz:
    Matrix mult = matCria(mat1 -> lin, mat2 -> col);

    // Faz a multiplicacao:
    int i, j, k;
    for(i = 0; i < mat1 -> lin; i++)
        for(k = 0; k < mat2 -> col; k++)
            for(j = 0; j < mat1 -> col; j++)
                mult -> p[i][k] = mult -> p[i][k] + (mat1 -> p[i][j] * mat2 -> p[j][k]);

    return mult;
}

Matrix matDecomposicaoLU(Matrix upper)
{
    // Prerequisitos (Matriz quadrada nao nula):
    if((upper == NULL) || (upper -> lin != upper -> col))
        return NULL;

    // Cria uma nova matriz para receber os fatores:
    Matrix lower = matCria(upper -> lin, upper -> lin);

    // Decompoe:
    int i, j;
    for(i = 0; i < upper -> lin - 1; i++)
    {
        for(j = i + 1; j < upper -> lin; j++)
        {
            // Determina o multiplicador:
            double mult = upper -> p[j][i] / upper -> p[i][i];

            // Salva o multiplicador em 'lower':
            lower -> p[j][i] = mult;

            // Multiplica a linha:
            matTransformaLinha(upper, j, i, mult);
        }
    }

    // Coloca 1's na diagonal principal de 'lower':
    for(i = 0; i < lower -> lin; i++)
        lower -> p[i][i] = 1;

    // Retorna a triangular inferior:
    return lower;
}

Matrix matDecomposicaoPivotLU(Matrix upper, Matrix P)
{
    // Prerequisitos (Matriz quadrada nao nula):
    if((upper == NULL) || (P == NULL) || (upper -> lin != upper -> col) || (P -> lin != upper -> lin) || (P -> col != upper -> col))
        return NULL;

    // Cria uma nova matriz para receber os fatores:
    Matrix lower = matCria(upper -> lin, upper -> lin);

    // Decompoe:
    int i, j;
    for(i = 0; i < upper -> lin - 1; i++)
    {
        // Localiza o pivo:
        int pivo = matLocalizaPivo(upper, i, i);

        // Troca de linhas se necessario:
        if(pivo != i)
        {
            matTrocaLinhas(upper, pivo, i);
            matTrocaLinhas(lower, pivo, i);
            matTrocaLinhas(P, pivo, i);
        }

        // Eliminacao de gauss (se elemento na linha atual nao for zero):
        if(fabs(upper -> p[i][i]) != 0)
            for(j = i + 1; j < upper -> lin; j++)
            {
                // Determina o multiplicador:
                double mult = upper -> p[j][i] / upper -> p[i][i];

                // Salva o multiplicador em 'lower':
                lower -> p[j][i] = mult;

                // Multiplica a linha:
                matTransformaLinha(upper, j, i, mult);
            }
    }

    // Coloca 1's na diagonal principal de 'lower':
    for(i = 0; i < lower -> lin; i++)
        lower -> p[i][i] = 1;

    // Retorna a triangular inferior:
    return lower;
}

Matrix matSolucaoPivotLU(Matrix upper, Matrix lower, Matrix P, Matrix b)
{
    // Prerequisitos (Matrizes quadradas nao nulas, vetor do tamanho da matriz):
    if((upper == NULL) || (lower == NULL) || (P == NULL) || (upper -> lin != upper -> col)
        || (lower -> lin != lower -> col) || (lower -> lin != upper -> lin)
        || (P -> lin != P -> col) || (P -> lin != upper -> lin))
        return NULL;

    // Realiza o produto matricial entre a matriz de permutacoes e a matriz solucao:
    Matrix Pb = matProdutoMatricial(P, b);

    // Substituicoes sucessivas:
    Matrix aux = matSubstSucessiva(lower, Pb);

    // Substituicoes retroativas:
    Matrix res = matSubstRetroativa(upper, aux);

    // Libera memoria:
    matLibera(aux);
    matLibera(Pb);

    // Retorna a matriz solucao:
    return res;
}

Matrix matSubstSucessiva(Matrix lower, Matrix b)
{
    // Prerequisitos (Matriz quadrada nao nula, vetor do tamanho da matriz):
    if((lower == NULL) || (lower -> lin != lower -> col))
        return NULL;
    
    // Cria uma matriz (1, n) auxiliar para armazenar 'y' em (Ly = b):
    Matrix aux = matCria(1, lower -> col);

    // Primeira linha e usada como padrao:
    aux -> p[0][0] = b -> p[0][0] / lower -> p[0][0];

    // Itera sobre cada linha (1 variavel por linha):
    int i, j;
    for(i = 1; i < lower -> lin; i++)
    {
        // Soma cada variavel conhecida com seu coeficiente:
        double soma = 0;
        for(j = 0; j < i; j++)
            soma += lower -> p[i][j] * aux -> p[0][j];

        // Resolve a igualdade e armazena:
        aux -> p[0][i] = (b -> p[i][0] - soma) / lower -> p[i][i];
    }

    return aux;
}

Matrix matSubstRetroativa(Matrix upper, Matrix y)
{
    // Prerequisitos (Matriz quadrada nao nula, vetor(1, n)):
    if((upper == NULL) || (upper -> lin != upper -> col) || (y -> lin != 1) || (y -> col != upper -> col))
        return NULL;

    // Cria uma matriz (1, n) para armazenar 'x' em (Ux = y):
    Matrix res = matCria(1, upper -> lin);

    // Ultima linha e usada como padrao:
    res -> p[0][y -> col - 1] = y -> p[0][y -> col - 1] / upper -> p[upper -> col - 1][upper -> col - 1];

    // Itera sobre cada linha, de baixo para cima:
    int i, j;
    for(i = upper -> col - 2; i >= 0; i--)
    {
        // Soma cada variavel conhecida com seu coeficiente:
        double soma = 0;
        for(j = i; j < upper -> col; j++)
            soma += upper -> p[i][j] * res -> p[0][j];

        // Resolve a igualdade e armazena:
        res -> p[0][i] = (y -> p[0][i] - soma) / upper -> p[i][i];
    }

    return res;
}

Matrix matInversa(Matrix mat)
{
    // Prerequisitos (Matriz quadrada nao nula):
    if((mat == NULL) || (mat -> lin != mat -> col))
        return NULL;

    // Cria matriz identidade para ser usada como conjunto solucao:
    Matrix identidade = matIdentidade(mat -> lin);

    // Cria matriz que armazenara a inversa:
    Matrix inversa = matCria(mat -> lin, mat -> col);

    // Determina as matrizes: lower, upper, e permutacao:
    Matrix P = matIdentidade(mat -> lin);
    Matrix upper = matCopia(mat);
    Matrix lower = matDecomposicaoPivotLU(upper, P);

    // Resolve n sistemas, para cada coluna de identidade:
    int i, j;
    Matrix b = matCria(mat -> lin, 1);
    for(i = 0; i < mat -> col; i++)
    {
        // Resolva o sistema para a coluna atual da identidade:
        int k;
        for(k = 0; k < mat -> lin; k++)
            b -> p[k][0] = identidade -> p[k][i];

        Matrix colInversa = matSolucaoPivotLU(upper, lower, P, b);

        // Copia o resultado do sistema linear para a matriz inversa:
        for(j = 0; j < inversa -> lin; j++) {
            double valor = colInversa -> p[0][j];
            // Trata erro de arredondamento:
            if(valor == 0) {
                valor = 0;
            }
            inversa -> p[j][i] = valor;
        }

        matLibera(colInversa);
    }

    matLibera(b);
    matLibera(P);
    matLibera(upper);
    matLibera(lower);
    matLibera(identidade);
    return inversa;
}

Matrix matGetColuna (Matrix mat, int coluna) {
    // Prerequisitos (Matriz nao nula, coluna valida):
    if((mat == NULL) || (coluna < 0) || (coluna > mat -> col))
        return NULL;

    Matrix col = matCria(mat -> lin, 1);
    int i;
    double valor;
    for(i = 0; i < mat -> lin; i++) {
        valor = matGet(mat, i, coluna);
        matPut(col, i, 0, valor);
    }

    return col;
}

double matDeterminante (Matrix mat)
{
    // Prerequisitos (Matriz quadrada nao nula):
    if((mat == NULL) || (mat -> lin != mat -> col))
        return -1;

    // Cria uma matriz:
    Matrix aux = matCopia(mat);

    // Chama triangulacao superior:
    double det = matSuperior(aux);

    // Libera matriz:
    matLibera(aux);

    // Retorna determinante obtido da triangulacao:
    return det;
}

double matSuperior (Matrix mat)
{
    // Prerequisitos (Matriz quadrada nao nula):
    if((mat == NULL) || (mat -> lin != mat -> col))
        return -1;

    int i, j;
    double det = 1;

    // Eliminacao de Gauss:
    for(i = 0; i < mat -> lin - 1; i++)
    {
        for(j = i + 1; j < mat -> lin; j++)
        {
            double mult = mat -> p[j][i] / mat -> p[i][i];
            matTransformaLinha(mat, j, i, mult);
        }
    }
    
    // Calcula o determinante (diagonal principal):
    for(i = 0; i < mat -> lin; i++)
        det = det * mat -> p[i][i];

    return det;
}

// Com pivo:
double matSuperiorPivot (Matrix mat)
{
    // Prerequisitos (Matriz quadrada nao nula):
    if((mat == NULL) || (mat -> lin != mat -> col))
        return -1;

    int i, j, pivo;
    double det = 1;

    for(i = 0; i < mat -> lin - 1; i++)
    {
        pivo = matLocalizaPivo(mat, i + 1, i);

        if(pivo != i)
        {
            matTrocaLinhas(mat, pivo, i);
            det = det * -1;
        }

        if(fabs(mat -> p[i][i]) != 0)
            for(j = i + 1; j < mat -> lin; j++)
            {
                double mult = mat -> p[j][i] / mat -> p[i][i];
                int k;
                for(k = 0; k < mat -> lin; k++)
                    mat -> p[j][k] = mat -> p[j][k] - (mult * mat -> p[i][k]);
            }
    }
    
    // Calcula o determinante (diagonal principal):
    for(i = 0; i < mat -> lin; i++)
        det = det * mat -> p[i][i];

    return det;
}

double matGet(Matrix mat, int lin, int col)
{
    // Prerequisitos (Matriz nao nula):
    if(mat == NULL)
        return -1;

    return mat -> p[lin][col];
}

double matProdutoEscalar (Matrix mat1, Matrix mat2)
{
    // Prerequisitos (Matrizes nao nulas, mat1(m, 1), mat2(1, n), m = n):
    if((mat1 == NULL) || (mat2 == NULL) || (mat1 -> col != 1) || (mat2 -> lin != 1) || (mat1 -> lin != mat2 -> col))
        return -1;

    // Faz a multiplicacao:
    double resultado = 0;
    int i;
    for(i = 0; i < mat1 -> lin; i++)
        resultado += mat1 -> p[i][0] * mat2 -> p[0][i];

    return resultado;
}

double* matReferenciaLinha (Matrix mat, int lin)
{
    // Prerequisitos (Matrize nao nula):
    if(mat == NULL)
        return 0;

    return mat -> p[lin];
}

int matNcolunas (Matrix mat)
{
    // Prerequisitos (Matrize nao nula):
    if(mat == NULL)
        return 0;

    return mat -> col;
}

int matNlinhas (Matrix mat)
{
    // Prerequisitos (Matrize nao nula):
    if(mat == NULL)
        return 0;

    return mat -> lin;
}

int matIgual (Matrix mat1, Matrix mat2)
{
    // Prerequisitos (Matrizes nao nulas, numero de linhas e colunas igual):
    if((mat1 == NULL) || (mat2 == NULL) || (mat1 -> lin != mat2 -> lin) || (mat1 -> col != mat2 -> col))
        return 0;

    // Compara cada elemento das matrizes:
    int i, j;
    for(i = 0; i < mat1 -> lin; i++)
        for(j = 0; j < mat1 -> col; j++)
            if(mat1 -> p[i][j] != mat2 -> p[i][j])
                return 0;

    return 1;
}

int matLocalizaPivo (Matrix mat , int lin, int col)
{
    // Prerequisitos (Matriz nao nula, linha e coluna validos):
    if((mat == NULL) || (lin < 0) || (lin >= mat -> lin) || (col < 0) || (col >= mat -> col))
        return -1;

    int i;
    int pivo = lin; // A linha inicial sera o pivo.
    double maior = fabs(mat -> p[lin][col]);

    // Percorre as demais linhas procurando um valor maior:
    for(i = lin + 1; i < mat -> lin; i++)
        if(fabs(mat -> p[i][col]) >= maior)
        {
            maior = fabs(mat -> p[i][col]); // Atualiza o valor do maior elemento.
            pivo = i; // Atualiza a posicao da linha com o pivo.
        }

    return pivo;
}

void matSetNcolunas (Matrix mat, int n)
{
    // Prerequisitos (Matriz nao nula, n positivo):
    if((mat == NULL) || (n <= 0))
        return;

    mat -> col = n;
    return;
}

void matSetNlinhas (Matrix mat, int n)
{
    // Prerequisitos (Matriz nao nula, n positivo):
    if((mat == NULL) || (n <= 0))
        return;

    mat -> lin = n;
    return;
}

void matLibera(Matrix mat)
{
    // Prerequisitos (Matriz nao nula):
    if(mat == NULL)
        return;

    // Libera a memoria de cada linha:
    int i;
    for(i = 0; i < mat -> lin; i++)
        free(mat -> p[i]);

    // Libera vetor de ponteiros:
    free(mat -> p);

    // Libera a matriz:
    free(mat);

    return;
}

void matLiberaNmatrizes(Matrix* vet, int n)
{
    // Prerequisitos (Vetor nao nulo, n positivo):
    if((vet == NULL) || (n <= 0))
        return;

    // Libera cada matriz:
    int i;
    for(i = 0; i < n; i++)
        matLibera(vet[i]);

    // Libera o vetor:
    free(vet);

    return;
}

void matImprime(Matrix mat)
{
    // Prerequisitos (Matriz nao nula):
    if(mat == NULL)
        return;

    // Percorre a matriz e imprime cada elemento:
    int i,j;
    for(i = 0; i < mat->lin; i++)
    {
        printf("| ");
        for(j = 0; j < mat -> col - 1; j++)
            printf("%lf\t", mat -> p[i][j]);

        // Nova linha para cada nova linha da matriz:
        printf("%lf |\n", mat -> p[i][j]);
    }
    //printf("\n");

    return;
}

void matPut(Matrix mat, int lin, int col, double valor)
{
    // Prerequisitos (Matriz nao nula, posicao valida):
    if((mat == NULL) || (lin < 0) || (lin >= mat -> lin) || (col < 0))
        return;

    mat -> p[lin][col] = valor;
    return;
}

void matSalva(Matrix mat, char nomeArquivo[])
{
    // Prerequisitos (Matriz nao nula):
    if(mat == NULL)
        return;

    // Abre o arquivo para salvar:
    FILE *psaida = fopen(nomeArquivo, "w");
    if(psaida == NULL)
    {
        printf("Espaço em disco insuficiente.\n");
        return;
    }

    // Imprime n.o de linhas e colunas:
    fprintf(psaida, "%d\t%d\n", mat -> lin, mat -> col);

    // Imprime cada elemento da matriz:
    int i, j;
    for(i = 0; i < mat->lin; i++)
    {
        for(j = 0; j < mat -> col; j++)
            fprintf(psaida, "%f\t", mat -> p[i][j]);
        fprintf(psaida, "\n");
    }

    // Fecha o arquivo:
    fclose(psaida);

    return;
}

void matTrocaLinhas (Matrix mat, int lin1, int lin2)
{
    // Prerequisitos (Matriz nao nula, linhas validas):
    if((mat == NULL) || (lin1 < 0) || (lin1 >= mat -> lin) || (lin2 < 0) || (lin2 >= mat -> lin))
        return;

    // Faz a troca:
    double *aux;
    aux = mat -> p[lin1];
    mat -> p[lin1] = mat -> p[lin2];
    mat -> p[lin2] = aux;

    return;
}

void matMultiplicaEscalar(Matrix mat, double escalar)
{
    // Prerequisitos (Matriz nao nula):
    if(mat == NULL)
        return;

    // Multiplica cada elemento da matriz:
    int i, j;
    for(i = 0; i < mat -> lin; i++)
        for(j = 0; j < mat -> col; j++)
            mat -> p[i][j] *= escalar;

    return;
}

void matMultiplicaLinhaEscalar (Matrix mat, int lin, double escalar)
{
    // Prerequisitos (Matriz nao nula, linha valida):
    if((mat == NULL) || (lin < 0) || (lin >= mat -> lin))
        return;

    // Percorre toda a linha e multiplica:
    int i;
    for(i = 0; i < mat -> col; i++)
        mat -> p[lin][i] = mat -> p[lin][i] * escalar;

    return;
}

void matTransformaLinha (Matrix mat , int linAlvo, int lin, double escalar)
{
    // Prerequisitos (Matriz nao nula, linhas validas):
    if((mat == NULL) || (linAlvo < 0) || (linAlvo >= mat -> lin) || (lin < 0) || (lin >= mat -> lin))
        return;

    int i;
    for(i = 0; i < mat -> col; i++)
        mat -> p[linAlvo][i] = mat -> p[linAlvo][i] - (mat -> p[lin][i] * escalar);

    return;
}