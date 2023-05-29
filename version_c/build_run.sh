
gcc-11 -pg -std=c11 -Wall  -o main -O0 -DNDEBUG -I ./include ./src/*.c main.c -lm
./main > /dev/null

gprof ./main gmon.out -b> prof/report.txt

if [ -f "prof/report.txt" ]; then
    echo "[01] - Arquivo gerado em: prof/report.txt"
    echo "[02] - Gerando arquivo de Perfilamento grafico.. "
    gprof2dot prof/report.txt > prof/report.dot
    dot -Tpng -o prof/report.png prof/report.dot
else 
    echo "Erro - arquivo gmon.out (ou report.txt) Inexistente!"
fi;
