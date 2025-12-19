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
        print("Достигнуто максимальное количество итерация.")

    return x, iter_count

def is_square_matrix(matrix):
    n = len(matrix)
    if n == 0:
        return False
    for row in matrix:
        if len(row) != n:
            return False
    return True

def determinant(matrix: list[list[int]]) -> int:
    det = 0
    match ln := len(matrix):
        case 1:
            return matrix[0][0]
        case 2:
            return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]
        case _:
            for k in range(ln):
                det += matrix[0][k] * (-1) ** k * determinant([matrix[i][:k] + matrix[i][k+1:] for i in range(1,ln)])
    return det

if __name__ == "__main__":
    A = read_matrix('A.txt')
    b = read_vector('B.txt')

    if not (is_square_matrix(A)):
      print("Матрица А не является квадратной")
      exit(1);

    if (determinant(A) == 0):
      print("Матрица А является вырожденной")
      exit(1);

    # Проверка ненулевой диагонали
    for i in range(len(b)):
        if A[i][i] == 0:
            raise ValueError(f"Диагональный элемент A[{i}][{i}] равен 0.")


    solution, iterations = relaxation_method(A, b, eps=1e-9)
    for i, val in enumerate(solution, start=1):
        print(f"x{i} = {val:.4g}")
    print(f"Итераций: {iterations}")
