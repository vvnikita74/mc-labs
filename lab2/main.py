from __future__ import annotations

from typing import Optional

def read_matrix(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    matrix = []
    for line in lines:
        row = list(map(float, line.split()))
        if row:
            matrix.append(row)
    return matrix

def read_vector(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    vec = []
    for line in lines:
        val = line.strip()
        if val:
            vec.append(float(val))
    return vec

def relaxation_method(A, b, eps=1e-6, max_iter=10000):
    n = len(b)

    # Подготовка P и c
    P = [[0.0 for _ in range(n)] for _ in range(n)]
    c = [0.0 for _ in range(n)]
    for i in range(n):
        c[i] = b[i] / A[i][i]
        for j in range(n):
            P[i][j] = -1.0 if i == j else -A[i][j] / A[i][i]

    # Начальное приближение
    x = [0.0 for _ in range(n)]
    R = c[:]  # копия c

    iter_count = 0
    while iter_count < max_iter:
        # Находим индекс максимальной по модулю невязки
        max_r = abs(R[0])
        s = 0
        for i in range(1, n):
            if abs(R[i]) > max_r:
                max_r = abs(R[i])
                s = i

        if max_r < eps:
            break

        delta = R[s]          # величина поправки
        x[s] += delta         # обновление переменной

        # Обновляем остальные невязки
        for i in range(n):
            if i != s:
                R[i] += P[i][s] * delta
        R[s] = 0.0

        iter_count += 1

    if iter_count >= max_iter:
        print("Достигнуто максимальное количество итераций.")

    return x, iter_count

def is_square_matrix(matrix):
    n = len(matrix)
    if n == 0:
        return False
    for row in matrix:
        if len(row) != n:
            return False
    return True


def is_diagonally_dominant(
    A: list[list[float]],
    tol: float = 0.0,
) -> bool:
    """
    Проверка диагонального преобладания по строкам:
      |a_ii| > sum_{j != i} |a_ij|
    tol — допуск для вещественных чисел.
    """
    n = len(A)
    for i in range(n):
        diag = abs(A[i][i])
        others = 0.0
        for j in range(n):
            if j != i:
                others += abs(A[i][j])
        if not (diag > others + tol):
            return False
    return True


def determinant(matrix: list[list[float]]) -> float:
    if not matrix or not matrix[0]:
        raise ValueError("Нельзя вычислить определитель пустой матрицы.")
    ln = len(matrix)
    if ln == 1:
        return matrix[0][0]
    if ln == 2:
        return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]

    det = 0.0
    for k in range(ln):
        minor = [matrix[i][:k] + matrix[i][k + 1 :] for i in range(1, ln)]
        det += matrix[0][k] * ((-1.0) ** k) * determinant(minor)
    return det


def validate_inputs(
    A: list[list[float]],
    b: list[float],
    *,
    eps: float,
    max_iter: int,
    dominance_tol: float = 0.0,
) -> None:
    n = len(A)
    if len(b) != n:
        raise ValueError(f"Размерность b ({len(b)}) не совпадает с размерностью A ({n}x{n}).")

    # Проверка ненулевой диагонали
    for i in range(n):
        if A[i][i] == 0:
            raise ValueError(f"Диагональный элемент A[{i}][{i}] равен 0.")

    if not is_diagonally_dominant(A, tol=dominance_tol):
        raise ValueError("Матрица A не удовлетворяет условию диагонального преобладания")

    if not is_square_matrix(A):
        raise ValueError("Матрица A не является квадратной")

    # Проверка невырожденности
    det = determinant(A)
    if abs(det) < 1e-12:
        raise ValueError("Матрица A является вырожденной")


def solve_relaxation(
    *,
    A: Optional[list[list[float]]] = None,
    b: Optional[list[float]] = None,
    A_file: str = "A.txt",
    B_file: str = "B.txt",
    eps: float = 1e-9,
    max_iter: int = 10000,
    dominance_tol: float = 0.0,
) -> tuple[list[float], int]:
    if A is None:
        A = read_matrix(A_file)
    if b is None:
        b = read_vector(B_file)

    validate_inputs(
        A,
        b,
        eps=eps,
        max_iter=max_iter,
    )
    return relaxation_method(A, b, eps=eps, max_iter=max_iter)

if __name__ == "__main__":
    try:
        solution, iterations = solve_relaxation(
            A_file="A.txt",
            B_file="B.txt",
            eps=1e-9,
        )
        for i, val in enumerate(solution, start=1):
            print(f"x{i} = {val:.4g}")
        print(f"Итераций: {iterations}")
    except ValueError as e:
        print(str(e))
        raise SystemExit(1)
