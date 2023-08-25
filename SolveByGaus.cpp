//
//  SolveByGaus.cpp
//  Сплайны
//
//  Created by Arsenii Klevakin on 11.07.2023.
//

#include "SolveByGaus.hpp"

#include <iostream>


using namespace std;

//Здесь M[строка][столбец]
vector<double> solveEquations(vector<vector<double>>& A, vector<double>& B) {
    int n = A.size();
    vector<vector<double>> augmentedMatrix(n, vector<double>(n + 1));

    // Создание расширенной матрицы
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            augmentedMatrix[i][j] = A[i][j];
        }
        augmentedMatrix[i][n] = B[i];
    }

    // Прямой ход метода Гаусса
    for (int i = 0; i < n; i++) {
        if (augmentedMatrix[i][i] == 0.0) {
            // Если ведущий элемент равен нулю, меняем строки
            int j;
            for (j = i + 1; j < n; j++) {
                if (augmentedMatrix[j][i] != 0.0) {
                    break;
                }
            }
            if (j == n) {
                // Система уравнений не имеет решения
                return vector<double>();
            }
            for (int k = 0; k <= n; k++) {
                swap(augmentedMatrix[i][k], augmentedMatrix[j][k]);
            }
        }

        for (int j = i + 1; j < n; j++) {
            double ratio = augmentedMatrix[j][i] / augmentedMatrix[i][i];
            for (int k = i; k <= n; k++) {
                augmentedMatrix[j][k] -= ratio * augmentedMatrix[i][k];
            }
        }
    }

    // Обратный ход метода Гаусса
    vector<double> solution(n);
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0.0;
        for (int j = i + 1; j < n; j++) {
            sum += augmentedMatrix[i][j] * solution[j];
        }
        solution[i] = (augmentedMatrix[i][n] - sum) / augmentedMatrix[i][i];
    }

    return solution;
}
