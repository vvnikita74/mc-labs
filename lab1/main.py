from __future__ import annotations

import math
import cmath

# Метод Ньютона (касательных) для нахождения одного корня уравнения 4-ой степени
def newton_method_single(a: float, b: float, c: float, d: float, k: float, x0: float = 1, tolerance: float = 1e-10, max_iterations: int = 1000) -> float:
    # Пробуем разные начальные точки
    # Сначала пробуем 0.0 (важно для уравнений вида x^n = 0), затем x0, потом плотную сетку
    initial_points: list[float] = [0.0, x0] + [x / 2 for x in range(-200, 201)]

    best_x: float = x0
    best_fx: float = abs(a * x0**4 + b * x0**3 + c * x0**2 + d * x0 + k)

    for start_point in initial_points:
        x: float = start_point

        for _ in range(max_iterations):
            fx: float = a * x**4 + b * x**3 + c * x**2 + d * x + k
            fpx: float = 4 * a * x**3 + 3 * b * x**2 + 2 * c * x + d

            # Если нашли достаточно точное решение
            if abs(fx) < tolerance:
                return x

            if abs(fpx) < 1e-15:
                # Производная близка к нулю, пробуем следующую начальную точку
                break

            x_new: float = x - fx / fpx

            # Проверяем сходимость
            if abs(x_new - x) < tolerance:
                # Проверяем, насколько хорошо это решение
                fx_new: float = a * x_new**4 + b * x_new**3 + c * x_new**2 + d * x_new + k
                if abs(fx_new) < abs(best_fx):
                    best_x = x_new
                    best_fx = abs(fx_new)
                if abs(fx_new) < tolerance:
                    return x_new
                break

            x = x_new

        # Сохраняем лучшее найденное решение
        fx_final: float = a * x**4 + b * x**3 + c * x**2 + d * x + k
        if abs(fx_final) < abs(best_fx):
            best_x = x
            best_fx = abs(fx_final)

    return best_x

# Деление полинома 4-ой степени на (x - root) методом Горнера (рекурсивная реализация)
# a - найденный корень
# b(n) - коэффициент многочлена (нового, 3-ей степени)
# a(n) - изначальный коэффициент
# b(n) = a(n)
# b(n-1) = a(n-1) + b(n) * a
def divide_polynomial(a: float, b: float, c: float, d: float, k: float, root: float) -> tuple[float, float, float, float]:
    def horner_step(coeffs: list[float], root: float, accumulated: list[float] = None) -> list[float]:
        if accumulated is None:
            accumulated = []

        if not coeffs:
            return accumulated

        if not accumulated:
            new_coeff = coeffs[0]
        else:
            new_coeff = coeffs[0] + accumulated[-1] * root

        return horner_step(coeffs[1:], root, accumulated + [new_coeff])

    result = horner_step([a, b, c, d], root)
    return tuple(result)

# Метод Ньютона с понижением степени для нахождения всех корней уравнения 4-ой степени
def newton_method(a: float, b: float, c: float, d: float, k: float, x0: float = 1, tolerance: float = 1e-10, max_iterations: int = 1000) -> list[float | complex]:
    # Находим первый корень методом Ньютона
    root1: float = newton_method_single(a, b, c, d, k, x0, tolerance, max_iterations)
    print(f"Первый корень (метод Ньютона): x1 = {root1:.4f}")

    # Делим полином на (x - root1) для получения кубического полинома
    b3, c3, d3, k3 = divide_polynomial(a, b, c, d, k, root1)
    print(f"Понижение степени: получен кубический полином")
    print(f"  {b3:.4f}x³ + {c3:.4f}x² + {d3:.4f}x + {k3:.4f} = 0")

    # Решаем кубическое уравнение методом Кардано
    cubic_roots: list[float | complex] = cardano_method(b3, c3, d3, k3)

    # Собираем все корни вместе
    all_roots: list[float | complex] = [root1] + cubic_roots

    return all_roots

# Метод Кардано для решения кубического уравнения
def cardano_method(b: float, c: float, d: float, k: float) -> list[float | complex]:
    # Приводим к виду x^3 + px + q = 0
    # Делим на b
    a2: float = c / b
    a1: float = d / b
    a0: float = k / b

    # Подстановка x = t - a2/3 для устранения квадратного члена
    p: float = a1 - a2**2 / 3
    q: float = 2 * a2**3 / 27 - a1 * a2 / 3 + a0

    # Вычисляем дискриминант для формулы Кардано
    Q: float = p / 3
    R: float = q / 2
    D: float = Q**3 + R**2

    roots: list[float | complex]

    if abs(D) < 1e-10:
        # Кратные корни (D = 0)
        if abs(R) < 1e-10:
            roots = [0.0 - a2 / 3] * 3
        else:
            # Формула Кардано для D = 0: x1 = 2*cuberoot(-R), x2 = x3 = -cuberoot(-R)
            neg_R: float = -R
            cube_neg_R_val: float = neg_R ** (1.0/3.0) if neg_R >= 0 else -(abs(neg_R) ** (1.0/3.0))
            cube_neg_R: float = cube_neg_R_val
            r1: float = 2 * cube_neg_R
            r2: float = -cube_neg_R
            roots = [r1 - a2 / 3, r2 - a2 / 3, r2 - a2 / 3]
    elif D > 0:
        # Один действительный корень и два комплексных
        sqrtD: float = math.sqrt(D)
        val_S: float = -R + sqrtD
        val_T: float = -R - sqrtD
        S_val: float = val_S ** (1.0/3.0) if val_S >= 0 else -(abs(val_S) ** (1.0/3.0))
        T_val: float = val_T ** (1.0/3.0) if val_T >= 0 else -(abs(val_T) ** (1.0/3.0))
        S: float = S_val
        T: float = T_val

        x1_sqrtD: float = S + T - a2 / 3
        x2_sqrtD: complex = -(S + T) / 2 - a2 / 3 + 1j * math.sqrt(3) * (S - T) / 2
        x3_sqrtD: complex = -(S + T) / 2 - a2 / 3 - 1j * math.sqrt(3) * (S - T) / 2

        roots = [x1_sqrtD, x2_sqrtD, x3_sqrtD]
    else:
        # Три различных действительных корня
        theta: float = math.acos(-R / math.sqrt(-Q**3))
        x1_theta: float = 2 * math.sqrt(-Q) * math.cos(theta / 3) - a2 / 3
        x2_theta: float = 2 * math.sqrt(-Q) * math.cos((theta + 2 * math.pi) / 3) - a2 / 3
        x3_theta: float = 2 * math.sqrt(-Q) * math.cos((theta + 4 * math.pi) / 3) - a2 / 3

        roots = [x1_theta, x2_theta, x3_theta]

    return roots

# Решение квадратного уравнения cx^2 + dx + k = 0 по формуле Виета
def vieta_quadratic(c: float, d: float, k: float) -> list[float | complex]:
    discriminant: float = d**2 - 4 * c * k

    x1: float | complex
    x2: float | complex

    if discriminant < 0:
        # Комплексные корни
        x1 = (-d + cmath.sqrt(discriminant)) / (2 * c)
        x2 = (-d - cmath.sqrt(discriminant)) / (2 * c)
    else:
        x1 = (-d + math.sqrt(discriminant)) / (2 * c)
        x2 = (-d - math.sqrt(discriminant)) / (2 * c)

    return [x1, x2]

def remove_duplicate_roots(roots: list[float | complex]) -> list[float | complex]:
    unique_roots = []
    seen = set()
    for root in roots:
        if isinstance(root, complex):
            rounded = (round(root.real, 4), round(root.imag, 4))
        else:
            rounded = round(root, 4)

        if rounded not in seen:
            seen.add(rounded)
            unique_roots.append(root)
    return unique_roots

# Решение уравнения 4-ой степени ax^4 + bx^3 + cx^2 + dx + k = 0
def solve_quartic(a: float, b: float, c: float, d: float, k: float) -> list[float | complex] | str:
    # Случай 1: a != 0 - уравнение 4-ой степени
    if a != 0:
        print("Решение уравнения 4-ой степени методом Ньютона с понижением степени")
        roots: list[float | complex] = newton_method(a, b, c, d, k)
        return roots

    # Случай 2: a == 0, b != 0 - кубическое уравнение
    elif b != 0:
        print("Решение кубического уравнения методом Кардано")
        roots: list[float | complex] = cardano_method(b, c, d, k)
        return roots

    # Случай 3: a == 0, b == 0, c != 0 - квадратное уравнение
    elif c != 0:
        print("Решение квадратного уравнения по формуле Виета")
        roots = vieta_quadratic(c, d, k)
        return roots

    # Случай 4: a == 0, b == 0, c == 0, d != 0 - линейное уравнение
    elif d != 0:
        print("Решение линейного уравнения")
        x: float = -k / d
        return [x]

    # Случай 5: все коэффициенты == 0
    else:
        if k == 0:
            print("Любое число является решением (0 = 0)")
            return "Любое число"
        else:
            print("Решений нет (противоречие)")
            return []

def main() -> None:
    print("Решение уравнения 4-ой степени: ax^4 + bx^3 + cx^2 + dx + k = 0\n")

    try:
        a: float = float(input("Введите коэффициент a: "))
        b: float = float(input("Введите коэффициент b: "))
        c: float = float(input("Введите коэффициент c: "))
        d: float = float(input("Введите коэффициент d: "))
        k: float = float(input("Введите коэффициент k: "))

        print(f"\nУравнение: {a}x^4 + {b}x^3 + {c}x^2 + {d}x + {k} = 0\n")

        roots: list[float | complex] | str = solve_quartic(a, b, c, d, k)

        if roots == "Любое число":
            print("Решение: любое число")
        elif len(roots) == 0:
            print("Решений нет")
        else:
            roots = remove_duplicate_roots(roots)

            print("\nКорни уравнения:")
            for i, root in enumerate(roots, 1):
                if isinstance(root, complex):
                    if root.imag == 0:
                        print(f"x{i} = {root.real:.4f}")
                    else:
                        # Форматирование комплексного числа
                        real_part = root.real
                        imag_part = root.imag
                        if abs(real_part) < 1e-10:  # Если действительная часть ~ 0
                            if imag_part > 0:
                                print(f"x{i} = {imag_part:.4f}i")
                            else:
                                print(f"x{i} = {imag_part:.4f}i")
                        else:
                            if imag_part >= 0:
                                print(f"x{i} = {real_part:.4f} + {imag_part:.4f}i")
                            else:
                                print(f"x{i} = {real_part:.4f} - {abs(imag_part):.4f}i")
                else:
                    print(f"x{i} = {root:.4f}")

    except ValueError:
        print("Ошибка: введите числовые значения")
    except Exception as e:
        print(f"Ошибка: {e}")

if __name__ == "__main__":
    main()
