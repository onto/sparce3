#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>

using namespace std;

int main(int argc, char *argv[4]) {

    if (argc != 4) {
        cout << "Need gamma, dimension and type of matrix." << endl;
        return 0;
    }

    int G = atoi(argv[1]);
    int N = atoi(argv[2]);
    int T = atoi(argv[3]);

    ofstream M("matrix.txt");
    ofstream V("vector.txt");

    M << N << endl;
    V << N << endl;

    time_t theTime;

    srand(time(&theTime));

    int k = int(N/(G-1));
    if (k == 0) k = 1;
    int c;

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {

            if (T == 0) { //certain

                if (i == j) {
                    M << j << " " << (1 + rand()%100) << "\t";
                } else {

                    c = rand()%k;
                    if (c == 0) {
                        M << j << " " << (1 + rand()%100) << "\t";
                    }
                }

            } else { //band

                if ((i-floor(0.5*G) <= j) && (j < i+ceil(0.5*G))) {
                    M << j << " " << (1 + rand()%100) << "\t";
                }
            }
        }
        M << -1 << endl;

        V << (rand()%2)*(rand()%100) << "\t";
    }
    M.close();
    V.close();

    //cout << t;

    return 0;
}

