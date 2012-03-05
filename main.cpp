#include <sparce.h>
#include <time.h>

using namespace std;

Coordinates exactOptimizedGeneralElement(SparceMatrix& A, int k, double pivTol) {

    Coordinates C;
    C.col = C.row = -1;

    int D = A.dimension(), q, j;
    const int n = D - k;

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

    double norm = 0;

    for (int i = k; i < D; ++i) {

        q = A.R[i].F;
        while (true) {
            if (q == -1)
                break;
            j = A.R[i].C[q];
            if (j >= k) {
                B[i - k][j - k] = 1;
                T[i - k] += abs(A.R[i].V[q]);
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

Coordinates optimizedGeneralElement(SparceMatrix& A, int k, double pivTol) {

    Coordinates C;
    C.col = C.row = -1;

    int D = A.dimension();
    const int n = D - k;

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
        while (true) {
            if (q == -1)
                break;
            j = A.R[i].C[q];
            if (j >= k) {
                I[i - k] += 1;
                J[j - k] += 1;
                T[i - k] += abs(A.R[i].V[q]);
            }
            q = A.R[i].N[q];
        }
    }

    double norm = 0;

    for (int i = 0; i < n; ++i) {
        if (T[i] > norm)
            norm = T[i];
        if ((I[i] == 0) || (J[i] == 0))
            return C;
    }

    int minc = D * D, newc;
    double value, v;

    for (int i = k; i < D; ++i) {

        q = A.R[i].F;
        while (q != -1) {
            j = A.R[i].C[q];
            if (j >= k) {
                v = abs(A.R[i].V[q]);
                newc = (I[i - k] - 1) * (J[j - k] - 1);
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

//SparceVector solveByGauss(SparceMatrix A, SparceVector B, double PivTol, bool exact) {

//    int N = A.dimension();
//    double t, t1, t2, t3;
//    SparceVector X(N);
//    int P[N], k;

//    for (int i = 0; i < N; ++i)
//        P[i] = i;

//    for (int i = 0; i < N - 1; ++i) {

//        Coordinates g = exact ? exactOptimizedGeneralElement(A, i, PivTol) : optimizedGeneralElement(A, i, PivTol);

//        if (g.col == -1)
//            throw "Can't solve by Gauss";

//        A.swapColRow(i, g.row, i, g.col);
//        B.swap(i, g.row);
//        if (i != g.col) {
//            k = P[i];
//            P[i] = P[g.col];
//            P[g.col] = k;
//        }

//        t = A.get(i, i);

//        for (int k = i + 1; k < N; ++k) {
//            t3 = A.get(k, i);
//            if (t3 == 0)
//                continue;
//            for (int j = i + 1; j < N; ++j) {
//                t2 = A.get(i, j);
//                if (t2 == 0)
//                    continue;
//                t1 = A.get(k, j);
//                A.set(k, j, t1 - t2 * t3 / t);
//            }
//            A.remove(k, i);
//            t2 = B.get(i);
//            if (t2 == 0)
//                continue;
//            t1 = B.get(k);
//            B.set(k, t1 - t2 * t3 / t);
//        }
//    }

//    for (int i = N - 1; i >= 0; --i) {
//        t3 = B.get(i);
//        if (t3 == 0)
//            continue;
//        t = A.get(i, i);
//        for (int j = i - 1; j >= 0; --j) {
//            t1 = A.get(j, i);
//            if (t1 == 0)
//                continue;
//            t2 = B.get(j);
//            B.set(j, t2 - t1 * t3 / t);
//            A.remove(j, i);
//        }
//        X.set(P[i], t3 / t);
//    }

//    return X;
//}

LUP LU(SparceMatrix A, double pivTol, bool exact) {

    int D = A.dimension();
    double t,t1,t2, t3;
    LUP T;
    T.L = SparceMatrix(D);
    T.Pr = SparceMatrix(D);
    T.Pc = SparceMatrix(D);

    for (int i = 0; i < D; ++i) {
        T.Pr.set(i, i, 1);
        T.Pc.set(i, i, 1);
    }

    for (int i = 0; i < D; ++i) {

        Coordinates g = exact ? exactOptimizedGeneralElement(A, i, pivTol) : optimizedGeneralElement(A, i, pivTol);

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
    cout << M.cellsNum() << endl;

//    V.fromFile("vector.txt");

    for (double i = 1; i >= 1e-10; i *= 0.1) {

        try {

            time_t t = clock();
            LUP T = LU(M, i, 0);
            t = clock() - t;
            T.Pr = T.Pr.transpose();
            T.Pc = T.Pc.transpose();
            M1 = T.Pr * T.L * T.U * T.Pc - M;
//            cout << T.L <<  endl << T.U << endl << T.Pr << endl << T.Pc << endl;
//            cout << M1 << endl;
//            M1 = M1 - M;
//            cout << M1 << endl;
            cout << "PivTol: " << i << "\t Error: " << M1.normM() << "\t cells: " << T.L.cellsNum() + T.U.cellsNum() << "\t time: " << double(t)/CLOCKS_PER_SEC <<endl;

        } catch (const char s[]) {
            cout << "PivTol: " << i << "\t Error: " << s << endl;
        }

        try {
            time_t t = clock();
            LUP T = LU(M, i, 1);
            t = clock() - t;
            T.Pr = T.Pr.transpose();
            T.Pc = T.Pc.transpose();
            M1 = T.Pr * T.L * T.U * T.Pc - M;
            cout << "PivTol: " << i << "\t Error: " << M1.normM() << "\t cells: " << T.L.cellsNum() + T.U.cellsNum() << "\t time: " << double(t)/CLOCKS_PER_SEC << endl;

        } catch (const char s[]) {
            cout << "PivTol: " << i << "\t Error: " << s << endl;
        }

        cout << endl;

    }

    return 0;
}
