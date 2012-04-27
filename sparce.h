#ifndef SPARCE_H
#define SPARCE_H

#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

struct Coordinates {

    int col;
    int row;
};

class SparceVector {

public:
    SparceVector();
    SparceVector(int dimension);
    SparceVector(const SparceVector& V1);

    void fromFile(const char file[]);

    void remove(int coordinate);
    void set(int coordinate, double value);
    double get(int coordinate) const;
    void add(int coordinate, double value);
    void swap(int c1, int c2);

    int dimension();
    int cellsNum() const;
    double normInf();

    friend ostream& operator <<(ostream& os, SparceVector& V1);
    SparceVector operator -(const SparceVector& V1);
    SparceVector operator +(const SparceVector& V1);
    SparceVector operator *(const double x);
    SparceVector& operator =(const SparceVector& M1);
    double operator *(const SparceVector& V1);

    vector<double> V; //values
    vector<int> N; //next in row
    vector<int> C; //column
    int F; //first
    int D; //dimension
};


class PermutationMatrix {

public:
    PermutationMatrix();
    PermutationMatrix(int dimension);

    void swapCol(int c1, int c2);
    void swapRow(int r1, int r2);
    PermutationMatrix transpose();

    PermutationMatrix& operator =(const PermutationMatrix& M1);
    friend ostream& operator <<(ostream& os, PermutationMatrix& M);

    int D;
    vector<int> M;
};

class SparceMatrix {

public:
    SparceMatrix();
    SparceMatrix(int dimension);
    SparceMatrix(const SparceMatrix& M1);

    void fromFile(const char file[]);

    void remove(int row, int col);
    void set(int row, int col, double value);
    double get(int row, int col) const;
    void add(int row,int col, double value);

    int dimension();
    int cellsNum() const;
    double normM();

    void swapCol(int c1, int c2);
    void swapRow(int r1, int r2);
    void swapColRow(int r1, int r2, int c1, int c2);
    SparceMatrix transpose() const;

    SparceMatrix operator +(const SparceMatrix& M1);
    SparceMatrix operator -(const SparceMatrix& M1);
    SparceMatrix operator *(const double x);
    SparceMatrix operator *(const SparceMatrix& M1);
    SparceMatrix operator *(const PermutationMatrix& M1);
    friend SparceMatrix operator *(const PermutationMatrix& M1, const SparceMatrix& M2);
    SparceVector operator *(const SparceVector& V1);
    SparceMatrix& operator =(const SparceMatrix& M1);
    friend ostream& operator <<(ostream& os, SparceMatrix& M);


    vector<SparceVector> R; // rows
    int D; //dimension
};

struct LUP {

    SparceMatrix L;
    SparceMatrix U;
    PermutationMatrix Pc;
    PermutationMatrix Pr;
};

#endif // SPARCE_H
