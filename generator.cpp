#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>

using namespace std;

int main(int argc, char *argv[3]) {

    if (argc != 3) {
        cout << "Need gamma and dimension of matrix." << endl;
        return 0;
    }

    int N = atoi(argv[2]);
    int G = atoi(argv[1]);

    ofstream M("matrix.txt");
    ofstream V("vector.txt");

    M << N << endl;
    V << N << endl;

    time_t theTime;

    srand(time(&theTime));

    //int t = 0;
    int k = int(N/(G-1));
    int c;

    for (int i = 1; i <= N; ++i) {
        for (int j = 1; j <= N; ++j) {

            if (i == j) {
                M << (1 + rand()%100) << "\t";
            } else {

                c = rand()%k;
                if (c != 0) {
                    M << 0 << "\t";
                } else {
                    M << (1 + rand()%100) << "\t";
                }
            }
        }
        M << endl;
        V << (rand()%2)*(rand()%100) << "\t";
    }
    M.close();
    V.close();

    //cout << t;

    return 0;
}

