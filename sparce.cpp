#include "sparce.h"

double eps = 0.0;

//##################//
//                  //
//   SparceMatrix   //
//                  //
//##################//

SparceMatrix::SparceMatrix() {

    D = 0;
}

SparceMatrix::SparceMatrix(int dimension) {

    D = dimension;
    R = new SparceVector[D];
}

SparceMatrix::SparceMatrix(const SparceMatrix& M1) {

    D = M1.D;
    R = M1.R;

}

int SparceMatrix::dimension() {

    return D;
}

int SparceMatrix::cellsNum() const {

    int sum = 0;
    for (int i = 0; i < D; ++i)
        sum += R[i].cellsNum();
    return sum;
}

void SparceMatrix::fromFile(const char file[]) {

    ifstream in(file);

    in >> D;
    R = new SparceVector[D];

    int q, x;

    for (int row = 0; row < D; ++row) {
        R[row].V.clear();
        R[row].C.clear();
        R[row].N.clear();
        R[row].F = 0;
        R[row].D = D;

        q = 0;
        for (int i = 0; i < D; ++i) {
            in >> x;
            if (x != 0) {
                R[row].C.push_back(i);
                R[row].V.push_back(x);
                R[row].N.push_back(q+1);
                q++;
            }

        }
        if (R[row].N.size() > 0)
            R[row].N.at(R[row].N.size()-1) = -1;
    }
    in.close();
}

void SparceMatrix::add(int row, int col, double value) {

    R[row].add(col,value);
}

void SparceMatrix::set(int row, int col, double value) {

    R[row].set(col,value);
}

double SparceMatrix::get(int row, int col) const {

    return R[row].get(col);
}

void SparceMatrix::remove(int row, int col) {

    R[row].remove(col);
}

void SparceMatrix::swapCol(int c1, int c2) {

    for (int i = 0; i < D; ++i)
        R[i].swap(c1,c2);
}

void SparceMatrix::swapRow(int r1, int r2) {

    SparceVector T;

    T = R[r1];
    R[r1] = R[r2];
    R[r2] = T;
}

void SparceMatrix::swapColRow(int r1, int r2, int c1, int c2) {

    swapCol(c1, c2);
    swapRow(r1, r2);
}

SparceMatrix SparceMatrix::transpose() { // TODO: make faster

    SparceMatrix M(D);

    for (int i = 0; i < D; ++i)
        for (int j = 0; j < D; ++j)
            M.add(j, i, get(i, j));

    return M;
}

double SparceMatrix::normM() {

    double norm = 0, sum;

    for (int i = 0; i < D; ++i) {

        sum = 0;
        for (int j = 0; j < R[i].cellsNum(); ++j)
            sum += abs(R[i].V[j]);
        if (sum > norm)
            norm = sum;
    }
    return norm;
}

//double SparceMatrix::normL() {

//}

ostream& operator <<(ostream& os, SparceMatrix& M1) {

    int D = M1.dimension();

    for (int i = 0; i < D; ++i) {
        os << M1.R[i] << endl;
    }
    return os;
}

SparceMatrix SparceMatrix::operator +(const SparceMatrix& M1) {

    SparceMatrix M(D);

    for (int i = 0; i < D; ++i)
        for (int j = 0; j < D; ++j)
            M.R[i] = R[i] + M1.R[i];

    return M;
}

SparceMatrix SparceMatrix::operator -(const SparceMatrix& M1) {

    SparceMatrix M(D);

    for (int i = 0; i < D; ++i)
        for (int j = 0; j < D; ++j)
            M.R[i] = R[i] - M1.R[i];

    return M;
}

SparceMatrix SparceMatrix::operator *(const double x) {

    SparceMatrix M(*this);

    for (int i = 0; i < D; ++i)
        for (int j = 0; j < M.R[i].cellsNum(); ++j)
            M.R[i].V[j] *= x;

    return M;
}

SparceMatrix SparceMatrix::operator *(const SparceMatrix& M1) { // TODO: make faster

    SparceMatrix M(D);
    double S;

    for (int i = 0; i < D; ++i)
        for (int j = 0; j < D; ++j) {
            S = 0;
            for (int k = 0; k < D; ++k) {
                S += get(i, k) * M1.get(k, j);
            }
            M.add(i, j, S);
        }

    return M;
}

SparceVector SparceMatrix::operator *(const SparceVector& V1) {

    SparceVector VV(D);

    for (int i = 0; i < D; ++i) {
        VV.add(i,R[i]*V1);
    }

    return VV;
}

//SparceMatrix& SparceMatrix::operator =(const SparceMatrix& M1) {

//    D = M1.D;
//    V = M1.V;
//    N = M1.N;
//    F = M1.F;
//    C = M1.C;

//    return *this;
//}

//##################//
//                  //
//   SparceVector   //
//                  //
//##################//

SparceVector::SparceVector() {

    D = 0;
    F = -1;
}

SparceVector::SparceVector(int dimension) {

    D = dimension;
    F = -1;
}

SparceVector::SparceVector(const SparceVector &V1) {

    V = V1.V;
    N = V1.N;
    C = V1.C;
    F = V1.F;
    D = V1.D;

}

void SparceVector::fromFile(const char file[]) {

    ifstream in(file);

    double x;

    V.clear();
    C.clear();
    N.clear();
    F = 0;

    in >> D;

    int q = 0;
    for (int i = 0; i < D; ++i) {
        in >> x;
        if (x != 0) {
            C.push_back(i);
            V.push_back(x);
            N.push_back(q+1);
            q++;
        }

    }
    if (N.size() > 0)
        N.at(N.size()-1) = -1;
    in.close();
}

void SparceVector::add(int coordinate, double value) {

    if (abs(value) <= eps)
        return;
    int q = F;
    int last = V.size();
    int i = q;
    bool first = true;

    V.push_back(value);
    C.push_back(coordinate);

    if (q != -1)
        while (C[q] < coordinate) {
            i = q;
            q = N[q];
            first = false;
            if (q == -1)
                break;
        }
    if (first) {
        F = last;
    } else {
        N[i] = last;
    }
    N.push_back(q);
}

void SparceVector::set(int coordinate, double value) { //TODO: make faster

    if (abs(value) <= eps) {
        remove(coordinate);
        return;
    }

    int q = F;

    while (true) {
        if (q == -1)
            break;
        if (C[q] == coordinate) {
            V[q] = value;
            return;
        }
        if (C[q] > coordinate) {
            add(coordinate, value);
            return;
        }
        q = N[q];
    }

    add(coordinate, value);
}

double SparceVector::get(int coordinate) const {

    int q = F;

    while (true) {
        if (q == -1)
            break;
        if (C[q] == coordinate)
            return V[q];
        if (C[q] > coordinate)
            return 0;
        q = N[q];
    }

    return 0;
}

void SparceVector::remove(int coordinate) {

    int q = F;
    int last = -1;

    while (true) {
        if (q == -1)
            break;
        if (C[q] == coordinate) {
            if (last == -1)
                F = N[q];
            else
                N[last] = N[q];
            C.erase(C.begin() + q);
            N.erase(N.begin() + q);
            V.erase(V.begin() + q);
            for (unsigned int i = 0; i < N.size(); ++i)
                if (N[i] >= q)
                    N[i] -= 1;
            if (F >= q)
                F -= 1;
            return;
        }
        last = q;
        q = N[q];
    }
}

void SparceVector::swap(int c1, int c2) { // TODO: make faster

    double t1 = get(c1), t2 = get(c2);

    set(c2, t1);
    set(c1, t2);

}

int SparceVector::dimension() {

    return D;
}

int SparceVector::cellsNum() const {

    return V.size();
}

double SparceVector::normInf() {

    double max = 0;

    for (unsigned int i = 0; i < V.size(); ++i)
        if (abs(V[i]) > max)
            max = abs(V[i]);

    return max;
}

ostream& operator <<(ostream& os, SparceVector& V1) { // TODO: make faster

    int D = V1.dimension();

    for (int i = 0; i < D; ++i) {
        os << V1.get(i) << "\t";
    }
    return os;
}

SparceVector SparceVector::operator +(const SparceVector& V1) { // TODO: make faster

    SparceVector VV = *this;

    for (int i = 0; i < V1.cellsNum(); ++i)
        VV.set(V1.C[i], VV.get(V1.C[i]) + V1.V[i]);

    return VV;
}

SparceVector SparceVector::operator -(const SparceVector& V1) { // TODO: make faster

    SparceVector VV = *this;

    for (int i = 0; i < V1.cellsNum(); ++i)
        VV.set(V1.C[i], VV.get(V1.C[i]) - V1.V[i]);

    return VV;
}

SparceVector SparceVector::operator *(const double x) {

    SparceVector VV = *this;

    for (int i = 0; i < cellsNum(); ++i)
        VV.V[i] *= x;

    return VV;
}

double SparceVector::operator *(const SparceVector& V1) { //TODO: make faster

    double sum = 0;
    for (int i = 0; i < D; ++i)
        sum += get(i)*V1.get(i);
    return sum;
}

