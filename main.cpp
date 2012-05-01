#include "sparce.h"
#include <time.h>
#include <stdlib.h>
#include <omp.h>

using namespace std;

const unsigned int MAX_EXACT  = 0;
const unsigned int MAX_SUM    = 1;
const unsigned int MAX_MULT   = 2;
const unsigned int MAX_MULT2  = 3;
const unsigned int NORM_EXACT = 4;
const unsigned int NORM_SUM   = 5;
const unsigned int NORM_MULT  = 6;
const unsigned int NORM_MULT2 = 7;

const unsigned int NUM_THREADS = 2;

Coordinates exactOptimizedGeneralElement(SparceMatrix& A, int k, double pivTol,
                                         bool byNorm) {

    int D = A.dimension();
    const int n = D - k;

    vector<vector<int> > B(n, vector<int> (n,0)),  B1(n, vector<int> (n,0)),
            G(n, vector<int> (n,0));

    unsigned int newc, minc = D * D, minc1 = D * D;
    int q, j;
    double v, S, norm = 0, value = 0, value1 = 0;

    Coordinates C, C1;
    C.col = C.row = C1.col = C1.row = -1;

#pragma omp parallel for private(q,j,v,S)
    for (int i = k; i < D; ++i) {

        S = 0;
        q = A.R[i].F;
        while (q != -1) {
            j = A.R[i].C[q];
            if (j >= k) {
                B[i-k][j-k] = 1;
                v = abs(A.R[i].V[q]);
                if (byNorm)
                    S += v;
                else if (v > S)
                    S = v;
            }
            q = A.R[i].N[q];
        }

#pragma omp critical
        {
        if (S > norm)
            norm = S;
        }
    }

#pragma omp parallel for
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            for (int k = 0; k < n; ++k)
                B1[i][j] += B[i][k] * (1 - B[j][k]); //B'[k][j] = 1 - B[j][k]

#pragma omp parallel for
    for (int i = 0; i < n; ++i)
        for (int k = 0; k < n; ++k)
            for (int j = 0; j < n; ++j)
                G[i][j] += B1[i][k] * B[k][j];

#pragma omp parallel for private(q,j,v,newc) firstprivate(minc1,value1,C1)
    for (int i = k; i < D; ++i) {

        q = A.R[i].F;
        while (q != -1) {
            j = A.R[i].C[q];
            if (j >= k) {
                v = abs(A.R[i].V[q]);
                newc = G[i-k][j-k];
                if ((v >= norm * pivTol) &&
                        ((newc < minc1) || ((newc == minc1) && (value1 < v)))) {
                    C1.row = i;
                    C1.col = j;
                    value1 = v;
                    minc1 = newc;
                }
            }
            q = A.R[i].N[q];
        }

#pragma omp critical
        {
        if ((minc1 < minc) || ((minc1 == minc) && (value < value1))) {
            C.row = C1.row;
            C.col = C1.col;
            value = value1;
            minc = minc1;
        }
        }
    }

    return C;
}

Coordinates optimizedGeneralElement(SparceMatrix& A, int k, double pivTol,
                                    bool byNorm, bool bySum, bool byMult2) {

    int D = A.dimension();
    const int n = D - k;

    vector<int> I(n,0), J(n,0);

    unsigned int newc, minc = D * D, minc1 = D * D;
    int q, j;
    double v, S, norm = 0, value = 0, value1 = 0;

    Coordinates C, C1;
    C.col = C.row = C1.col = C1.row = -1;

#pragma omp parallel for private(q,j,v,S)
    for (int i = k; i < D; ++i) {

        S = 0;
        q = A.R[i].F;
        while (q != -1) {
            j = A.R[i].C[q];
            if (j >= k) {
                I[i-k]++;
#pragma omp atomic
                J[j-k]++;
                v = abs(A.R[i].V[q]);
                if (byNorm)
                    S += v;
                else if (v > S)
                    S = v;
            }
            q = A.R[i].N[q];
        }

#pragma omp critical
        {
        if (S > norm)
            norm = S;
        }
    }

    for (int i = 0; i < n; ++i) {
        if ((I[i] == 0) || (J[i] == 0))
            return C;
    }

#pragma omp parallel for private(q,j,v,newc) firstprivate(minc1,value1,C1)
    for (int i = k; i < D; ++i) {

        q = A.R[i].F;
        while (q != -1) {
            j = A.R[i].C[q];
            if (j >= k) {
                v = abs(A.R[i].V[q]);
                newc = bySum ? (I[i-k]+J[j-k]) :
                               (byMult2 ? (I[i-k]*J[j-k]) :
                                          (I[i-k]-1)*(J[j-k]-1));
                if ((v >= norm * pivTol) &&
                        ((newc < minc1) || ((newc == minc1) && (value1 < v)))) {
                    C1.row = i;
                    C1.col = j;
                    value1 = v;
                    minc1 = newc;
                }
            }
            q = A.R[i].N[q];
        }

#pragma omp critical
        {
        if ((minc1 < minc) || ((minc1 == minc) && (value < value1))) {
            C.row = C1.row;
            C.col = C1.col;
            value = value1;
            minc = minc1;
        }
        }
    }

    return C;
}

Coordinates getOptimizedElement(SparceMatrix& A, int i, double pivTol,
                                int genElem) {

    switch (genElem) {
    case NORM_EXACT:
        return exactOptimizedGeneralElement(A,i,pivTol,1);
        break;
    case NORM_MULT:
        return optimizedGeneralElement(A,i,pivTol,1,0,0);
        break;
    case NORM_MULT2:
        return optimizedGeneralElement(A,i,pivTol,1,0,1);
        break;
    case NORM_SUM:
        return optimizedGeneralElement(A,i,pivTol,1,1,0);
        break;
    case MAX_EXACT:
        return exactOptimizedGeneralElement(A,i,pivTol,0);
        break;
    case MAX_MULT:
        return optimizedGeneralElement(A,i,pivTol,0,0,0);
        break;
    case MAX_MULT2:
        return optimizedGeneralElement(A,i,pivTol,0,0,1);
        break;
    case MAX_SUM:
        return optimizedGeneralElement(A,i,pivTol,0,1,0);
        break;
    }
    return exactOptimizedGeneralElement(A,i,pivTol,1);
}

LUP LU(SparceMatrix A, double pivTol, int genElem) {

    int D = A.dimension();
    double t,t1;
    LUP T;
    T.L = SparceMatrix(D);
    T.Pr = PermutationMatrix(D);
    T.Pc = PermutationMatrix(D);

    for (int i = 0; i < D; ++i) {

        Coordinates g = getOptimizedElement(A, i, pivTol, genElem);
        //genElem determines the way of choosing general element,
        //pivTol determines limit of value of general element

        if (g.col == -1)
            throw "Can not LU correctly";

        A.swapColRow(i, g.row, i, g.col);
        T.L.swapRow(i, g.row);
        T.Pr.swapRow(i, g.row);
        T.Pc.swapCol(i, g.col);

        t = A.get(i, i);
#pragma omp parallel for private(t1)
        for (int k = i + 1; k < D; ++k) {
            t1 = A.get(k, i);
            if (t1 != 0) {
                A.R[k] = A.R[k] - (A.R[i]*(t1/t)); //A.R[k] - k-th row of matrix A
                T.L.set(k, i, t1/t);
                A.remove(k, i);
            }
        }
        T.L.set(i, i, 1);
    }

    T.U = A;

    return T;
}

int main(int argc, char *argv[]) {

#ifdef _OPENMP
    omp_set_num_threads(NUM_THREADS);
#endif

    SparceMatrix M, M1;
//    SparceVector V, V1, V2;
    M.fromFile("matrix.txt");
//    cout << M.cellsNum() << endl;

//    for (double i = 1; i >= 1e-6; i *= 0.1) {

//        try {
//            time_t t = clock();
//            LUP T = LU(M, i, MAX_MULT);
//            t = clock() - t;
//            //T.Pr = T.Pr.transpose();
//            //T.Pc = T.Pc.transpose();
//            M1 = T.Pr * T.L * T.U * T.Pc - M;
//            cout << "Norm,Exact \t" << "PivTol: " << i << "\t Error: " << M1.normM() << "\t cells: " << T.L.cellsNum() + T.U.cellsNum() << "\t time: " << double(t)/CLOCKS_PER_SEC << endl;

//        } catch (const char s[]) {
//            cout << "Norm,Exact \t" << "PivTol: " << i << "\t Error: " << s << endl;
//        }
//        try {
//            time_t t = clock();
//            LUP T = LU(M, i, MAX_EXACT);
//            t = clock() - t;
//            //T.Pr = T.Pr.transpose();
//            //T.Pc = T.Pc.transpose();
//            M1 = T.Pr * T.L * T.U * T.Pc - M;
//            cout << "Norm,Exact \t" << "PivTol: " << i << "\t Error: " << M1.normM() << "\t cells: " << T.L.cellsNum() + T.U.cellsNum() << "\t time: " << double(t)/CLOCKS_PER_SEC << endl;

//        } catch (const char s[]) {
//            cout << "Norm,Exact \t" << "PivTol: " << i << "\t Error: " << s << endl;
//        }

//        cout << endl;

//    }

//    double pivTol = atof(argv[1]);
//    ofstream out("output.txt",ios_base::app);
//    ofstream tout("toutput.txt",ios_base::app);
//    out << M.cellsNum() << " ";
//    for (int i = 0; i < 8; i++) {
//        if ((i == 0) || (i == 4)) {
//            out << "0 ";
//            tout << "-1 ";
//            continue;
//        }
//        try {
//            time_t t = clock();
//            LUP T = LU(M, pivTol, i);
//            tout << double(clock() - t)/CLOCKS_PER_SEC << " ";
//            out << T.L.cellsNum() + T.U.cellsNum() - T.L.D << " ";
//        } catch (const char s[]) {
//            out << "0 ";
//            tout << "-1 ";
//        }
//    }
//    out << endl;
//    tout << endl;

//    double i = 0.01;
//    try {
//        time_t t = clock();
//        LUP T = LU(M, i, MAX_MULT2);
//        t = clock() - t;
//        //T.Pr = T.Pr.transpose();
//        //T.Pc = T.Pc.transpose();
//        M1 = T.Pr * T.L * T.U * T.Pc - M;
//        cout << "Norm,Exact \t" << "PivTol: " << i << "\t Error: " << M1.normM() << "\t cells: " << T.L.cellsNum() + T.U.cellsNum() << "\t time: " << double(t)/CLOCKS_PER_SEC << endl;

//    } catch (const char s[]) {
//        cout << "Norm,Exact \t" << "PivTol: " << i << "\t Error: " << s << endl;
//    }

    //double pivTol = atof(argv[1]);
//    double pivTol = 0.01;
//    try {
//        time_t t = clock();
//        LUP T = LU(M, pivTol, MAX_MULT2);
//        t = clock() - t;
//        cout << "Norm,Exact \t" << "PivTol: " << pivTol << "\t Error: " << M1.normM() << "\t cells: " << T.L.cellsNum() + T.U.cellsNum() << "\t time: " << double(t)/CLOCKS_PER_SEC << endl;

//    } catch (const char s[]) {
//        cout << "Error: " << s << endl;
//    }

    cout << M << endl;

    return 0;
}
