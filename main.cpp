#include <sparce.h>
#include <time.h>

using namespace std;

const unsigned int MAX_EXACT  = 1;
const unsigned int MAX_MULT   = 2;
const unsigned int MAX_SUM    = 3;
const unsigned int NORM_EXACT = 4;
const unsigned int NORM_MULT  = 5;
const unsigned int NORM_SUM   = 6;


Coordinates exactOptimizedGeneralElement(SparceMatrix& A, int k, double pivTol, bool byNorm) {

    Coordinates C;
    C.col = C.row = -1;

    int D = A.dimension(), q, j;
    const int n = D - k;
    double norm = 0, t;

    int B[n][n], B1[n][n], G[n][n];

    double T[n];

    for (int i = 0; i < n; ++i) {
        T[i] = 0;
        for (int j = 0; j < n; ++j) {
            B[i][j] = 0;
            B1[i][j] = 0;
            G[i][j] = 0;
        }
    }



    for (int i = k; i < D; ++i) {

        q = A.R[i].F;
        while (q != -1) {
            j = A.R[i].C[q];
            if (j >= k) {
                B[i - k][j - k] = 1;
                t = abs(A.R[i].V[q]);
                if (byNorm) {
                    T[i - k] += t;
                } else {
                    if (t > norm) norm = t;
                }
            }
            q = A.R[i].N[q];
        }
        if (T[i - k] > norm)
            norm = T[i - k];
    }

    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            for (int k = 0; k < n; ++k)
                B1[i][j] += B[i][k] * (1 - B[j][k]); //B'[k][j] = 1 - B[j][k]

    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            for (int k = 0; k < n; ++k)
                G[i][j] += B1[i][k] * B[k][j];

    int minc = D * D, newc;
    double value, v;

    for (int i = k; i < D; ++i) {

        q = A.R[i].F;
        while (q != -1) {
            j = A.R[i].C[q];
            if (j >= k) {
                v = abs(A.R[i].V[q]);
                newc = G[i - k][j - k];
                if ((v >= norm * pivTol) && ((newc < minc) || ((newc == minc) && (value < v)))) {
                    C.row = i;
                    C.col = j;
                    value = v;
                    minc = newc;
                }
            }
            q = A.R[i].N[q];
        }
    }

    return C;
}

Coordinates optimizedGeneralElement(SparceMatrix& A, int k, double pivTol, bool byNorm, bool bySum) {

    Coordinates C;
    C.col = C.row = -1;

    int D = A.dimension();
    const int n = D - k;
    double norm = 0, t;

    int I[n], J[n];
    double T[n];

    for (int i = 0; i < n; ++i) {
        I[i] = 0;
        J[i] = 0;
        T[i] = 0;
    }

    int q, j;

    for (int i = k; i < D; ++i) {

        q = A.R[i].F;
        while (q != -1) {
            j = A.R[i].C[q];
            if (j >= k) {
                I[i - k] += 1;
                J[j - k] += 1;
                t = abs(A.R[i].V[q]);
                if (byNorm) {
                    T[i - k] += t;
                } else {
                    if (t > norm) norm = t;
                }
            }
            q = A.R[i].N[q];
        }
        if ((byNorm) && (T[i - k] > norm))
            norm = T[i - k];
    }

    if (byNorm) {
        for (int i = 0; i < n; ++i) {
            if (T[i] > norm)
                norm = T[i];
            if ((I[i] == 0) || (J[i] == 0))
                return C;
        }
    }

    int minc = D * D, newc;
    double value, v;

    for (int i = k; i < D; ++i) {

        q = A.R[i].F;
        while (q != -1) {
            j = A.R[i].C[q];
            if (j >= k) {
                v = abs(A.R[i].V[q]);
                newc = bySum ? (I[i-k]+J[j-k]-1) :(I[i - k] - 1) * (J[j - k] - 1);
                if ((v >= norm * pivTol) && ((newc < minc) || ((newc == minc) && (value < v)))) {
                    C.row = i;
                    C.col = j;
                    value = v;
                    minc = newc;
                }
            }
            q = A.R[i].N[q];
        }
    }

    return C;
}

Coordinates getOptimizedElement(SparceMatrix& A, int i, double pivTol, int genElem) {

    switch (genElem) {
    case NORM_EXACT:
        return exactOptimizedGeneralElement(A,i,pivTol,1);
        break;
    case NORM_MULT:
        return optimizedGeneralElement(A,i,pivTol,1,0);
        break;
    case NORM_SUM:
        return optimizedGeneralElement(A,i,pivTol,1,1);
        break;
    case MAX_EXACT:
        return exactOptimizedGeneralElement(A,i,pivTol,0);
        break;
    case MAX_MULT:
        return optimizedGeneralElement(A,i,pivTol,0,0);
        break;
    case MAX_SUM:
        return optimizedGeneralElement(A,i,pivTol,0,1);
        break;
    }
    return exactOptimizedGeneralElement(A,i,pivTol,1);
}

LUP LU(SparceMatrix A, double pivTol, int genElem) {

    int D = A.dimension();
    double t,t3;
    LUP T;
    T.L = SparceMatrix(D);
    T.Pr = SparceMatrix(D);
    T.Pc = SparceMatrix(D);

    for (int i = 0; i < D; ++i) {
        T.Pr.set(i, i, 1);
        T.Pc.set(i, i, 1);
    }

    for (int i = 0; i < D; ++i) {

        Coordinates g = getOptimizedElement(A,i,pivTol,genElem);

        if (g.col == -1)
            throw "Can't LU correctly";

        A.swapColRow(i, g.row, i, g.col);
        T.L.swapColRow(i, g.row, i, g.col);
        T.Pr.swapRow(i, g.row);
        T.Pc.swapCol(i, g.col);

        t = A.get(i, i);
        for (int k = i + 1; k < D; ++k) {
            t3 = A.get(k, i);
            if (t3 == 0)
                continue;
            A.R[k] = A.R[k] - (A.R[i]*(t3/t));
            T.L.set(k, i, t3 / t);
            A.remove(k, i);
        }
        T.L.set(i, i, 1);
    }

    T.U = A;

    return T;
}

int main() {

    SparceMatrix M, M1;
//    SparceVector V;
    M.fromFile("matrix.txt");
//    cout << M.cellsNum() << endl;

//    V.fromFile("vector.txt");

//    for (double i = 1; i >= 1e-6; i *= 0.1) {

//        try {
//            time_t t = clock();
//            LUP T = LU(M, i, MAX_EXACT);
//            t = clock() - t;
//            T.Pr = T.Pr.transpose();
//            T.Pc = T.Pc.transpose();
//            M1 = T.Pr * T.L * T.U * T.Pc - M;
//            cout << "Norm,Exact \t" << "PivTol: " << i << "\t Error: " << M1.normM() << "\t cells: " << T.L.cellsNum() + T.U.cellsNum() << "\t time: " << double(t)/CLOCKS_PER_SEC << endl;

//        } catch (const char s[]) {
//            cout << "Norm,Exact \t" << "PivTol: " << i << "\t Error: " << s << endl;
//        }

//        cout << endl;

//    }

    double pivTol = 0.01;

    ofstream out("output.txt",ios_base::app);
    for (int i = 1; i <= 6; i++) {
        try {
            LUP T = LU(M, pivTol, i);
            out << T.L.cellsNum() + T.U.cellsNum() - T.L.D << " ";
        } catch (const char s[]) {
            out << 0 << " ";
        }
    }
    out << endl;


    return 0;
}
