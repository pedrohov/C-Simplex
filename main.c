/*
 * Pesquisa Operacional - 2o sem/2018
 * IFMG Campus Formiga
 * Pedro Henrique Oliveira Veloso (0002346)
 * Savio Silva Simoes (0011942)
 */

#include <stdlib.h>
#include <stdio.h>
#include "matrix.h"
#include "simplex.h"

int main(int argc, char *argv[]) {

    // Recebe parametros da linha de comando.
    // <modelo.mod> <arquivoSaida>
    if(argc <= 2) {
        printf("Especificar arquivo de modelo e local para saida de dados.");
        return 0;
    }

    Model modelo = carregaModelo(argv[1]);
    int solucao = simplex(modelo);
    outputModelo(modelo, solucao, argv[1], argv[2]);

    // Libera memoria utilizada pelo modelo:
    liberaModelo(modelo);

    return 0;
}
