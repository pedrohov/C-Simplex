/*
 * Pesquisa Operacional - 2o sem/2018
 * IFMG Campus Formiga
 * Pedro Henrique Oliveira Veloso (0002346)
 * Savio Silva Simoes (0011942)
 */

#ifndef MATRIX_H_INCLUDED
#define MATRIX_H_INCLUDED

// Tipo:
typedef struct Tmatrix *Matrix;

// Sub-rotinas:
Matrix  matCria (int lin, int col);                             // Cria uma matriz m(lin, col).
Matrix  matCopia (Matrix mat);                                  // Retorna uma copia profunda da matriz 'mat'.
Matrix  matCarrega (char nomeArquivo[]);                        // Carrega uma matriz de um arquivo.
Matrix  matIdentidade (int n);                                  // Gera a matriz identidade m(n, n).
Matrix  matVetMedia (Matrix mat);                               // Retorna uma matriz m(n, 1) com a media de cada coluna.    
Matrix  matCovariancia (Matrix mat);                            // Retorna a matriz covariancia de 'mat'.
Matrix  matTransposta (Matrix mat);                             // Retorna a matriz transposta de 'mat'.
Matrix  matOposta (Matrix mat);                                 // Retorna a matriz oposta de 'mat'.
Matrix  matSoma (Matrix mat1, Matrix mat2);                     // Retorna a soma das matrizes ou NULL (impossivel somar).
Matrix  matSubtrai (Matrix mat1, Matrix mat2);                  // Retorna a subtracao de 'mat1' - 'mat2' ou NULL (impossivel subtrair).
Matrix  matProdutoMatricial (Matrix mat1, Matrix mat2);         // Retorna o produto matricial de 'mat1' (x) 'mat2' ou NULL (impossivel multiplicar).
Matrix  matDecomposicaoLU (Matrix upper);                       // Modifica a matriz 'upper' em sua triangular superior. Retorna a triangular inferior.
Matrix  matDecomposicaoPivotLU (Matrix upper, Matrix P);        // Realiza a decomposicao LU com pivoteamento. Recebe a matriz identidade P.
Matrix  matSolucaoPivotLU(Matrix upper, Matrix lower, Matrix P, Matrix b);       // Retorna uma matriz m(1, n) com a solucao do sistema (mat, b). Matriz b(n, 1).
Matrix  matSubstSucessiva (Matrix lower, Matrix b);             // Retorna uma matriz m(1, n) com a solucao do sistema (lower, b) (Usado por SolucaoPivotLU).
Matrix  matSubstRetroativa (Matrix upper, Matrix y);            // Retorna a solucao de upper(n, n) por y(1, n) (Usado por SolucaoPivotLU).
Matrix  matInversa (Matrix mat);                                // Retorna a matriz inversa de 'mat' ou NULL (impossivel calcular).
Matrix  matGetColuna (Matrix mat, int coluna);                  // Retorna a coluna 'coluna' da matriz 'mat'.
double  matDeterminante (Matrix mat);                           // Retorna o determinante de uma matriz quadrada.
double  matSuperior (Matrix mat);                               // Realiza a triangulacao superior da matriz por Eliminacao de Gauss, retorna seu determinante.
double  matSuperiorPivot (Matrix mat);                          // Realiza a triangulacao por Eliminacao de Gauss com pivoteamento, retorna seu determinante.
double  matGet (Matrix mat, int lin, int col);                  // Retorna o valor na posicao m(lin, col).
double  matProdutoEscalar (Matrix mat1, Matrix mat2);           // Retorna o produto escalar entre mat1(m, 1) x mat2(1, n) ou -1 (impossivel multiplicar).
double* matReferenciaLinha (Matrix mat, int lin);               // Retorna um ponteiro para a linha 'lin' da matriz 'mat' que da acesso a todas as colunas da linha.
int     matNcolunas (Matrix mat);                               // Retorna o numero de colunas da matriz 'mat'.
int     matNlinhas (Matrix mat);                                // Retorna o numero de linhas da matriz 'mat'.
int     matIgual (Matrix mat1, Matrix mat2);                    // Retorna 1 se as matrizes forem iguais ou 0 se diferentes.
int     matLocalizaPivo (Matrix mat , int lin, int col);        // Procura pelo maior elemento igual ou abaixo da linha 'lin' na coluna 'col'.
void    matSetNcolunas (Matrix mat, int n);                     // Altera o numero de colunas de 'mat'.
void    matSetNlinhas (Matrix mat, int n);                      // Altera o numero de linhas de 'mat'.
void    matLibera (Matrix mat);                                 // Libera a memoria utilizada pela matriz 'mat'.
void    matLiberaNmatrizes(Matrix* vet, int n);                 // Libera a memoria de um vetor de matrizes.
void    matImprime (Matrix mat);                                // Imprime a matriz no console.
void    matPut (Matrix mat, int lin, int col, double valor);    // Altera o valor da posicao mat[lin][col].
void    matSalva (Matrix mat, char nomeArquivo[]);              // Salva a matriz em um arquivo.
void    matTrocaLinhas (Matrix mat, int lin1, int lin2);        // Troca as linhas 'lin1' e 'lin2' da matriz 'mat'.
void    matMultiplicaEscalar (Matrix mat, double escalar);      // Multiplica a matriz 'mat' por um numero real.
void    matMultiplicaLinhaEscalar (Matrix mat, int lin, double escalar);           // Multiplica a linha 'lin' por um numero real.
void    matTransformaLinha (Matrix mat , int linAlvo, int lin, double escalar);    // Faz [linAlvo - (lin * escalar)] para todos os elementos da linha (auxilia LU).

#endif // MATRIX_H_INCLUDED