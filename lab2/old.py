def read_matrix(filename):
    matrix = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line:
                row = list(map(float, line.split()))
                matrix.append(row)
    return matrix

def read_vector(filename):
    vector = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line:
                vector.append(float(line))
    return vector

def vector_norm_diff(x_new, x_old):
    """ Максимальная (бесконечная) норма разности векторов """
    max_diff = 0.0
    for i in range(len(x_new)):
        diff = abs(x_new[i] - x_old[i])
        if diff > max_diff:
            max_diff = diff
    return max_diff

def sor_method(A, b, omega=1.25, tol=1e-6, max_iter=1000):
    n = len(b)
    # Начальное приближение: все нули
    x = [0.0] * n

    for iteration in range(max_iter):
        x_old = x[:]  # Сохраняем копию решения с предыдущей итерации

        # Обновляем каждую компоненту x[i]
        for i in range(n):
            # Сумма по j = 0 до i-1: x[j] (x[j] = x_j^{(k+1)})
            sum_new = 0.0
            for j in range(0, i):
                sum_new += A[i][j] * x[j]

            # Сумма по j = i+1 до n-1: x[j] (x[j] = x_j^{(k)} = x_old[j])
            sum_old = 0.0
            for j in range(i + 1, n):
                sum_old += A[i][j] * x_old[j]

            x_gs = (b[i] - sum_new - sum_old) / A[i][i]

            # Применяем релаксацию
            x[i] = (1 - omega) * x_old[i] + omega * x_gs

        # Проверяем, изменилось ли решение мало?
        if vector_norm_diff(x, x_old) < tol:
            print(f"Сходимость достигнута на итерации {iteration + 1}")
            return x

    print("Достигнуто максимальное число итераций. Решение может быть неточным.")
    return x

def main():
    try:
        A = read_matrix('A.txt')
        b = read_vector('B.txt')
    except FileNotFoundError as e:
        print(f"Ошибка: файл не найден — {e}")
        return
    except Exception as e:
        print(f"Ошибка при чтении данных: {e}")
        return

    n = len(b)
    if len(A) != n:
        print("Ошибка: количество строк в A не совпадает с длиной B")
        return
    for i, row in enumerate(A):
        if len(row) != n:
            print(f"Ошибка: строка {i} матрицы A имеет неверную длину")
            return
        if A[i][i] == 0:
            print(f"Ошибка: диагональный элемент A[{i}][{i}] равен нулю")
            return

    # Параметры метода
    omega = 1.2
    tol = 1e-6
    max_iter = 1000

    solution = sor_method(A, b, omega=omega, tol=tol, max_iter=max_iter)

    print("\nРешение системы:")
    for i, val in enumerate(solution):
        print(f"x[{i + 1}] = {val:.4f}")

if __name__ == "__main__":
    main()
