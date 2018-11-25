# C-Simplex

Trabalho Prático da disciplina de Pesquisa Operacional. O trabalho consiste na identificação de um problema do mundo real que pode ser otimizado por meio de programação linear. O problema foi modelado na forma padrão, dividindo-o em variáveis e restrições de igualdade. Para verificar a solução ótima do problema de PL foi implementado um solver Simplex em C.

# Execução
No Windows é necessário ter a ferramenta make e o compilador gcc acessíveis pelo terminal.
Os parâmetros <modelo> e <saída> são obrigatórios.

```
make
simplex <modelo> <saída>
```
O Simplex resolverá o modelo na forma padrão passado pelo argumento e criará o arquivo de saída com informações da solução.
