/*
 * Pesquisa Operacional - 2o sem/2018
 * IFMG Campus Formiga
 * Pedro Henrique Oliveira Veloso (0002346)
 * Savio Silva Simoes (0011942)
 */

#ifndef SIMPLEX_H_INCLUDED
#define SIMPLEX_H_INCLUDED

#include "matrix.h"

// Tipo:
typedef struct Tmodel *Model;

// Sub-rotinas:
Model  carregaModelo (char arquivo[]); // Carrega modelos na forma canonica.
Model  criaModelo (Matrix A, Matrix b, Matrix c, Matrix base, Matrix naoBase, Matrix artificiais);
double calcObjetivo (Model modelo);    // Calcula a funcao objetivo do modelo para a solucao atual.
int    simplex (Model modelo);         // Procura por uma solucao unica otima utilizando o simplex.
int    artificialNaBase(Model modelo); // Verifica se nao existe solucao.
int    existePositivo (Matrix u);      // Verifica se a solucao e ilimitada.
int    existeNaoBase0 (Model modelo, Matrix custoReduzido); // Verifica se existem multiplos otimos.
void   imprimeModelo (Model modelo);   // Exibe os dados do modelo.
void   geraMatrizBase (Model modelo);  // Pega as colunas de A que pertencem a base e cria B.
void   liberaModelo (Model modelo);    // Libera a memoria utilizada pelo modelo.
void   outputModelo(Model modelo, int iteracoes, char entrada[], char saida[]); // Saida de dados em arquivo.

#endif
