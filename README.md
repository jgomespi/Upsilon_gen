# Upsilon_gen
# Pythia8 bottomonium event generation
Repositório com os programas e instruções para a geração, utilizando Pythia8, de eventos do decaimento do méson Upsilon em um par de múons, com o objetivo de calcular a aceptância do CMS.

Todas as informações relevantes para a geração dos eventos estão no arquivo gen.C, os demais arquivos apenas viabilizam a compilação e execução. O número de eventos desejado pode ser setado neste arquivo, na linha 54.


Para conseguir rodar na sua máquina, é necessário possuir o Pythia8 e o ROOT6 devidamente instalados.

Para rodar, basta fazer:

source go_gen;

que vai compilar o programa gen. Para executar, basta fazer:

./gen

