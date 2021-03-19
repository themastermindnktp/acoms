g++ main.cpp -o main

mkdir result

./main -i data/mcb/dataset.fa -w 6 -o result/mcb.positions.conf
#./main -i data/CRP/dataset.fa -w 22 -o result/CRP.positions.conf
#./main -i data/ERE_200/dataset.fa -w 13 -o result/ERE_200.positions.conf --recalc-background
#./main -i data/E2F200/dataset.fa -w 13 -o result/E2F200.positions.conf --recalc-background
#./main -i data/creb/dataset.fa -w 13 -o result/creb.positions.conf --recalc-background

