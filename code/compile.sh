g++ -march=core2 -pthread -std=c++11 -DSFMT_MEXP=607 -I /usr/local/include/ -L /usr/local/lib/ -O3 ./SFMT.c subsim.cpp -lnlopt -lm -o radeim
