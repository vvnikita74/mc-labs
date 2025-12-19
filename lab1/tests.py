from __future__ import annotations
from .main import newton_method, newton_method_single, cardano_method, vieta_quadratic, solve_quartic

import pytest

def format_root(root: float | complex) -> str:
    """Форматирование корня для красивого вывода"""
    if isinstance(root, complex):
        if abs(root.imag) < 1e-10:
            return f"{root.real:.2f}"
        elif abs(root.real) < 1e-10:
            return f"{root.imag:.2f}i"
        elif root.imag >= 0:
            return f"{root.real:.2f} + {root.imag:.2f}i"
        else:
            return f"{root.real:.2f} - {abs(root.imag):.2f}i"
    else:
        return f"{root:.2f}"

class TestNewtonMethod:
    """Тесты для метода Ньютона с понижением степени (уравнения 4-ой степени, a != 0)"""

    def test_simple_quartic(self) -> None:
        """Тест простого уравнения x^4 - 1 = 0, корни: ±1, ±i"""
        roots = newton_method(1, 0, 0, 0, -1, x0=0.5)
        print(f"\n  Найденные корни:")
        for i, root in enumerate(roots, 1):
            print(f"    x{i} = {format_root(root)}")
        assert len(roots) == 4
        # Проверяем, что каждый корень удовлетворяет уравнению
        for root in roots:
            value = 1 * root**4 - 1
            assert abs(value) < 1e-8

    def test_quartic_with_all_coefficients(self) -> None:
        """Тест уравнения x^4 - 5x^2 + 4 = 0, корни: -2, -1, 1, 2"""
        roots = newton_method(1, 0, -5, 0, 4, x0=1.5)
        print(f"\n  Найденные корни:")
        for i, root in enumerate(roots, 1):
            print(f"    x{i} = {format_root(root)}")
        assert len(roots) == 4
        # Проверяем, что каждый корень удовлетворяет уравнению
        for root in roots:
            value = 1 * root**4 - 5 * root**2 + 4
            assert abs(value) < 1e-7
        # Проверяем наличие реальных корней ±1 и ±2
        real_roots = sorted([r.real if isinstance(r, complex) else r for r in roots if isinstance(r, float) or abs(r.imag) < 1e-9])
        expected = [-2.0, -1.0, 1.0, 2.0]
        assert len(real_roots) == 4
        for i in range(4):
            assert abs(real_roots[i] - expected[i]) < 1e-6

    def test_quartic_x4_minus_16(self) -> None:
        """Тест уравнения x^4 - 16 = 0, корни: ±2, ±2i"""
        roots = newton_method(1, 0, 0, 0, -16, x0=1.5)
        print(f"\n  Найденные корни:")
        for i, root in enumerate(roots, 1):
            print(f"    x{i} = {format_root(root)}")
        assert len(roots) == 4
        # Проверяем, что каждый корень удовлетворяет уравнению
        for root in roots:
            value = 1 * root**4 - 16
            assert abs(value) < 1e-8

    def test_quartic_negative_start(self) -> None:
        """Тест с отрицательной начальной точкой x^4 - 1 = 0"""
        roots = newton_method(1, 0, 0, 0, -1, x0=-0.5)
        print(f"\n  Найденные корни:")
        for i, root in enumerate(roots, 1):
            print(f"    x{i} = {format_root(root)}")
        assert len(roots) == 4
        # Проверяем, что каждый корень удовлетворяет уравнению
        for root in roots:
            value = 1 * root**4 - 1
            assert abs(value) < 1e-8

    def test_quartic_large_coefficient(self) -> None:
        """Тест с большим коэффициентом: 200x^4 - 3200 = 0, корни: ±2, ±2i"""
        roots = newton_method(200, 0, 0, 0, -3200, x0=1.5)
        print(f"\n  Найденные корни:")
        for i, root in enumerate(roots, 1):
            print(f"    x{i} = {format_root(root)}")
        assert len(roots) == 4
        # 200x^4 = 3200 => x^4 = 16 => корни ±2, ±2i
        for root in roots:
            value = 200 * root**4 - 3200
            assert abs(value) < 1e-5

    def test_quartic_small_coefficients(self) -> None:
        """Тест с малыми коэффициентами: 0.001x^4 - 0.016 = 0, корни: ±2, ±2i"""
        roots = newton_method(0.001, 0, 0, 0, -0.016, x0=1.5)
        print(f"\n  Найденные корни:")
        for i, root in enumerate(roots, 1):
            print(f"    x{i} = {format_root(root)}")
        assert len(roots) == 4
        # 0.001x^4 = 0.016 => x^4 = 16 => корни ±2, ±2i
        for root in roots:
            value = 0.001 * root**4 - 0.016
            assert abs(value) < 1e-6

    def test_quartic_mixed_large_small(self) -> None:
        """Тест со смешанными коэффициентами: 4x^4 + 0.001x^2 - 64 = 0"""
        roots = newton_method(4, 0, 0.001, 0, -64, x0=2.0)
        print(f"\n  Найденные корни:")
        for i, root in enumerate(roots, 1):
            print(f"    x{i} = {format_root(root)}")
        assert len(roots) == 4
        # Проверяем, что каждый корень удовлетворяет уравнению
        for root in roots:
            value = 4 * root**4 + 0.001 * root**2 - 64
            assert abs(value) < 1e-5

    def test_quartic_all_coefficients_large(self) -> None:
        """Тест с большими коэффициентами: 100x^4 - 500x^3 + 600x^2 - 200x = 0"""
        roots = newton_method(100, -500, 600, -200, 0, x0=0.5)
        print(f"\n  Найденные корни:")
        for i, root in enumerate(roots, 1):
            print(f"    x{i} = {format_root(root)}")
        assert len(roots) == 4
        # Проверяем, что каждый корень удовлетворяет уравнению
        for root in roots:
            value = 100 * root**4 - 500 * root**3 + 600 * root**2 - 200 * root
            assert abs(value) < 1e-4


class TestCardanoMethod:
    """Тесты для метода Кардано (кубические уравнения, a = 0, b != 0)"""

    def test_simple_cubic(self) -> None:
        """Тест x^3 - 1 = 0, корень x = 1"""
        roots = cardano_method(1, 0, 0, -1)
        print(f"\n  Найденные корни:")
        for i, root in enumerate(roots, 1):
            print(f"    x{i} = {format_root(root)}")
        # Проверяем, что один из корней равен 1
        real_roots = [r for r in roots if isinstance(r, float) or abs(r.imag) < 1e-9]
        assert len(real_roots) > 0
        assert any(abs(r - 1.0) < 1e-9 for r in real_roots)

    def test_cubic_three_real_roots(self) -> None:
        """Тест x^3 - 6x^2 + 11x - 6 = 0, корни: 1, 2, 3"""
        roots = cardano_method(1, -6, 11, -6)
        print(f"\n  Найденные корни:")
        for i, root in enumerate(roots, 1):
            print(f"    x{i} = {format_root(root)}")
        real_roots = sorted([r.real if isinstance(r, complex) else r for r in roots])
        expected = [1.0, 2.0, 3.0]
        for i in range(3):
            assert abs(real_roots[i] - expected[i]) < 1e-6

    def test_cubic_x3_minus_8(self) -> None:
        """Тест x^3 - 8 = 0, один действительный корень x = 2"""
        roots = cardano_method(1, 0, 0, -8)
        print(f"\n  Найденные корни:")
        for i, root in enumerate(roots, 1):
            print(f"    x{i} = {format_root(root)}")
        real_roots = [r for r in roots if isinstance(r, float) or abs(r.imag) < 1e-9]
        complex_roots = [r for r in roots if isinstance(r, complex) and abs(r.imag) >= 1e-9]
        assert len(real_roots) == 1
        assert len(complex_roots) == 2
        assert any(abs(r - 2.0) < 1e-9 for r in real_roots)

    def test_cubic_with_complex_roots(self) -> None:
        """Тест x^3 + 1 = 0, корень x = -1 и два комплексных"""
        roots = cardano_method(1, 0, 0, 1)
        print(f"\n  Найденные корни:")
        for i, root in enumerate(roots, 1):
            print(f"    x{i} = {format_root(root)}")
        real_roots = [r for r in roots if isinstance(r, float) or abs(r.imag) < 1e-9]
        complex_roots = [r for r in roots if isinstance(r, complex) and abs(r.imag) >= 1e-9]
        assert len(real_roots) == 1
        assert len(complex_roots) == 2
        assert any(abs(r - (-1.0)) < 1e-9 for r in real_roots)

    def test_cubic_large_coefficients(self) -> None:
        """Тест с большими коэффициентами: 125x^3 - 1000 = 0, корень x = 2"""
        roots = cardano_method(125, 0, 0, -1000)
        print(f"\n  Найденные корни:")
        for i, root in enumerate(roots, 1):
            print(f"    x{i} = {format_root(root)}")
        real_roots = [r for r in roots if isinstance(r, float) or abs(r.imag) < 1e-9]
        complex_roots = [r for r in roots if isinstance(r, complex) and abs(r.imag) >= 1e-9]
        assert len(real_roots) == 1
        assert len(complex_roots) == 2
        assert any(abs(r - 2.0) < 1e-8 for r in real_roots)

    def test_cubic_small_coefficients(self) -> None:
        """Тест с малыми коэффициентами: 0.001x^3 - 0.008 = 0, корень x = 2"""
        roots = cardano_method(0.001, 0, 0, -0.008)
        print(f"\n  Найденные корни:")
        for i, root in enumerate(roots, 1):
            print(f"    x{i} = {format_root(root)}")
        real_roots = [r for r in roots if isinstance(r, float) or abs(r.imag) < 1e-9]
        complex_roots = [r for r in roots if isinstance(r, complex) and abs(r.imag) >= 1e-9]
        assert len(real_roots) == 1
        assert len(complex_roots) == 2
        assert any(abs(r - 2.0) < 1e-7 for r in real_roots)

    def test_cubic_mixed_coefficients(self) -> None:
        """Тест со смешанными коэффициентами: 50x^3 - 300x^2 + 550x - 300 = 0"""
        roots = cardano_method(50, -300, 550, -300)
        print(f"\n  Найденные корни:")
        for i, root in enumerate(roots, 1):
            print(f"    x{i} = {format_root(root)}")
        # Корни: 1, 2, 3
        real_roots = sorted([r.real if isinstance(r, complex) else r for r in roots])
        expected = [1.0, 2.0, 3.0]
        for i in range(3):
            assert abs(real_roots[i] - expected[i]) < 1e-5

    def test_cubic_negative_large_coefficient(self) -> None:
        """Тест с отрицательным большим коэффициентом: -10x^3 + 270 = 0, корень x = 3"""
        roots = cardano_method(-10, 0, 0, 270)
        print(f"\n  Найденные корни:")
        for i, root in enumerate(roots, 1):
            print(f"    x{i} = {format_root(root)}")
        real_roots = [r for r in roots if isinstance(r, float) or abs(r.imag) < 1e-9]
        complex_roots = [r for r in roots if isinstance(r, complex) and abs(r.imag) >= 1e-9]
        assert len(real_roots) == 1
        assert len(complex_roots) == 2
        assert any(abs(r - 3.0) < 1e-8 for r in real_roots)

    def test_cubic_discriminant_zero(self) -> None:
        """Тест с D = 0 (кратные корни): x^3 - 3x + 2 = 0, корни: -2, 1, 1"""
        roots = cardano_method(1, 0, -3, 2)
        print(f"\n  Найденные корни:")
        for i, root in enumerate(roots, 1):
            print(f"    x{i} = {format_root(root)}")
        # Корни: -2 (простой) и 1 (двукратный)
        real_roots = sorted([r.real if isinstance(r, complex) else r for r in roots])
        # Проверяем наличие корней
        assert abs(real_roots[0] - (-2.0)) < 1e-6
        assert abs(real_roots[1] - 1.0) < 1e-6
        assert abs(real_roots[2] - 1.0) < 1e-6


class TestVietaQuadratic:
    """Тесты для квадратных уравнений (a = 0, b = 0, c != 0)"""

    def test_simple_quadratic(self) -> None:
        """Тест x^2 - 4 = 0, корни: -2, 2"""
        roots = vieta_quadratic(1, 0, -4)
        print(f"\n  Найденные корни:")
        for i, root in enumerate(roots, 1):
            print(f"    x{i} = {format_root(root)}")
        roots_sorted = sorted([r.real if isinstance(r, complex) else r for r in roots])
        assert abs(roots_sorted[0] - (-2.0)) < 1e-9
        assert abs(roots_sorted[1] - 2.0) < 1e-9

    def test_quadratic_x2_minus_1(self) -> None:
        """Тест x^2 - 1 = 0, корни: -1, 1"""
        roots = vieta_quadratic(1, 0, -1)
        print(f"\n  Найденные корни:")
        for i, root in enumerate(roots, 1):
            print(f"    x{i} = {format_root(root)}")
        roots_sorted = sorted([r.real if isinstance(r, complex) else r for r in roots])
        assert abs(roots_sorted[0] - (-1.0)) < 1e-9
        assert abs(roots_sorted[1] - 1.0) < 1e-9

    def test_quadratic_complex_roots(self) -> None:
        """Тест x^2 + 1 = 0, корни: -i, i (комплексные)"""
        roots = vieta_quadratic(1, 0, 1)
        print(f"\n  Найденные корни:")
        for i, root in enumerate(roots, 1):
            print(f"    x{i} = {format_root(root)}")
        assert all(isinstance(r, complex) for r in roots)
        # Проверяем, что корни: i и -i
        assert any(abs(r - 1j) < 1e-9 for r in roots)
        assert any(abs(r - (-1j)) < 1e-9 for r in roots)

    def test_quadratic_with_all_coefficients(self) -> None:
        """Тест x^2 + 3x + 2 = 0, корни: -2, -1"""
        roots = vieta_quadratic(1, 3, 2)
        print(f"\n  Найденные корни:")
        for i, root in enumerate(roots, 1):
            print(f"    x{i} = {format_root(root)}")
        roots_sorted = sorted([r.real if isinstance(r, complex) else r for r in roots])
        assert abs(roots_sorted[0] - (-2.0)) < 1e-9
        assert abs(roots_sorted[1] - (-1.0)) < 1e-9

    def test_quadratic_double_root(self) -> None:
        """Тест x^2 - 2x + 1 = 0, корень: 1 (кратности 2)"""
        roots = vieta_quadratic(1, -2, 1)
        assert abs(roots[0] - 1.0) < 1e-9
        assert abs(roots[1] - 1.0) < 1e-9

    def test_quadratic_large_coefficients(self) -> None:
        """Тест с большими коэффициентами: 100x^2 - 400 = 0, корни: -2, 2"""
        roots = vieta_quadratic(100, 0, -400)
        print(f"\n  Найденные корни:")
        for i, root in enumerate(roots, 1):
            print(f"    x{i} = {format_root(root)}")
        roots_sorted = sorted([r.real if isinstance(r, complex) else r for r in roots])
        assert abs(roots_sorted[0] - (-2.0)) < 1e-8
        assert abs(roots_sorted[1] - 2.0) < 1e-8

    def test_quadratic_small_coefficients(self) -> None:
        """Тест с малыми коэффициентами: 0.0001x^2 - 0.0004 = 0, корни: -2, 2"""
        roots = vieta_quadratic(0.0001, 0, -0.0004)
        print(f"\n  Найденные корни:")
        for i, root in enumerate(roots, 1):
            print(f"    x{i} = {format_root(root)}")
        roots_sorted = sorted([r.real if isinstance(r, complex) else r for r in roots])
        assert abs(roots_sorted[0] - (-2.0)) < 1e-6
        assert abs(roots_sorted[1] - 2.0) < 1e-6

    def test_quadratic_mixed_large_small(self) -> None:
        """Тест со смешанными коэффициентами: 0.5x^2 + 100x + 200 = 0"""
        roots = vieta_quadratic(0.5, 100, 200)
        print(f"\n  Найденные корни:")
        for i, root in enumerate(roots, 1):
            print(f"    x{i} = {format_root(root)}")
        # Проверяем, что корни удовлетворяют уравнению
        for root in roots:
            r = root.real if isinstance(root, complex) else root
            value = 0.5 * r**2 + 100 * r + 200
            assert abs(value) < 1e-7

    def test_quadratic_large_complex_roots(self) -> None:
        """Тест с большими коэффициентами и комплексными корнями: 50x^2 + 200 = 0"""
        roots = vieta_quadratic(50, 0, 200)
        print(f"\n  Найденные корни:")
        for i, root in enumerate(roots, 1):
            print(f"    x{i} = {format_root(root)}")
        assert all(isinstance(r, complex) for r in roots)
        # Корни: ±2i
        imag_parts = sorted([abs(r.imag) for r in roots])
        assert abs(imag_parts[0] - 2.0) < 1e-8
        assert abs(imag_parts[1] - 2.0) < 1e-8


class TestLinearEquations:
    """Тесты для линейных уравнений (a = 0, b = 0, c = 0, d != 0)"""

    def test_linear_large_coefficient(self) -> None:
        """Тест с большим коэффициентом: 200x - 400 = 0, корень x = 2"""
        roots = solve_quartic(0, 0, 0, 200, -400)
        print(f"\n  Найденный корень: x = {format_root(roots[0])}")
        assert isinstance(roots, list)
        assert len(roots) == 1
        assert abs(roots[0] - 2.0) < 1e-9

    def test_linear_small_coefficient(self) -> None:
        """Тест с малым коэффициентом: 0.001x - 0.002 = 0, корень x = 2"""
        roots = solve_quartic(0, 0, 0, 0.001, -0.002)
        print(f"\n  Найденный корень: x = {format_root(roots[0])}")
        assert isinstance(roots, list)
        assert len(roots) == 1
        assert abs(roots[0] - 2.0) < 1e-9

    def test_linear_negative_large(self) -> None:
        """Тест с отрицательным большим коэффициентом: -50x + 150 = 0, корень x = 3"""
        roots = solve_quartic(0, 0, 0, -50, 150)
        print(f"\n  Найденный корень: x = {format_root(roots[0])}")
        assert isinstance(roots, list)
        assert len(roots) == 1
        assert abs(roots[0] - 3.0) < 1e-9

    def test_linear_fractional_result(self) -> None:
        """Тест с дробным результатом: 7x - 3 = 0, корень x = 3/7"""
        roots = solve_quartic(0, 0, 0, 7, -3)
        print(f"\n  Найденный корень: x = {format_root(roots[0])}")
        assert isinstance(roots, list)
        assert len(roots) == 1
        assert abs(roots[0] - 3.0/7.0) < 1e-9


class TestSolveQuartic:
    """Интеграционные тесты для solve_quartic"""

    def test_case_quartic_a_not_zero(self) -> None:
        """Случай 1: a != 0 - уравнение 4-ой степени"""
        roots = solve_quartic(1, 0, 0, 0, -16)
        assert isinstance(roots, list)
        assert len(roots) == 4
        # Проверяем, что все корни удовлетворяют уравнению
        for root in roots:
            value = 1 * root**4 - 16
            assert abs(value) < 1e-7

    def test_case_cubic_b_not_zero(self) -> None:
        """Случай 2: a = 0, b != 0 - кубическое уравнение"""
        roots = solve_quartic(0, 1, 0, 0, -8)
        assert isinstance(roots, list)
        assert len(roots) == 3
        # Проверяем наличие корня x = 2
        real_roots = [r for r in roots if isinstance(r, float) or abs(r.imag) < 1e-9]
        complex_roots = [r for r in roots if isinstance(r, complex) and abs(r.imag) >= 1e-9]
        assert len(real_roots) == 1
        assert len(complex_roots) == 2
        assert any(abs(r - 2.0) < 1e-9 for r in real_roots)

    def test_case_quadratic_c_not_zero(self) -> None:
        """Случай 3: a = 0, b = 0, c != 0 - квадратное уравнение"""
        roots = solve_quartic(0, 0, 1, 0, -4)
        print(f"\n  Найденные корни:")
        for i, root in enumerate(roots, 1):
            print(f"    x{i} = {format_root(root)}")
        assert isinstance(roots, list)
        assert len(roots) == 2
        roots_sorted = sorted([r.real if isinstance(r, complex) else r for r in roots])
        assert abs(roots_sorted[0] - (-2.0)) < 1e-9
        assert abs(roots_sorted[1] - 2.0) < 1e-9

    def test_case_linear_d_not_zero(self) -> None:
        """Случай 4: a = 0, b = 0, c = 0, d != 0 - линейное уравнение"""
        roots = solve_quartic(0, 0, 0, 2, -10)
        assert isinstance(roots, list)
        assert len(roots) == 1
        assert abs(roots[0] - 5.0) < 1e-9

    def test_case_all_zero_k_zero(self) -> None:
        """Случай 5: все коэффициенты = 0, k = 0 - любое число"""
        result = solve_quartic(0, 0, 0, 0, 0)
        assert result == "Любое число"

    def test_case_all_zero_k_not_zero(self) -> None:
        """Случай 5: все коэффициенты = 0, k != 0 - нет решений"""
        result = solve_quartic(0, 0, 0, 0, 5)
        assert isinstance(result, list)
        assert len(result) == 0

    def test_quadratic_complex_roots_full(self) -> None:
        """Квадратное уравнение с комплексными корнями через solve_quartic"""
        roots = solve_quartic(0, 0, 1, 0, 1)
        print(f"\n  Найденные корни:")
        for i, root in enumerate(roots, 1):
            print(f"    x{i} = {format_root(root)}")
        assert isinstance(roots, list)
        assert len(roots) == 2
        assert all(isinstance(r, complex) for r in roots)
        # Корни: i и -i
        assert any(abs(r - 1j) < 1e-9 for r in roots)
        assert any(abs(r - (-1j)) < 1e-9 for r in roots)

    def test_cubic_complex_coefficients(self) -> None:
        """Кубическое уравнение с комплексными корнями x^3 + 27 = 0"""
        roots = solve_quartic(0, 1, 0, 0, 27)
        assert isinstance(roots, list)
        assert len(roots) == 3
        # Один из корней должен быть -3, остальные два комплексные
        real_roots = [r for r in roots if isinstance(r, float) or abs(r.imag) < 1e-9]
        complex_roots = [r for r in roots if isinstance(r, complex) and abs(r.imag) >= 1e-9]
        assert len(real_roots) == 1
        assert len(complex_roots) == 2
        assert any(abs(r - (-3.0)) < 1e-9 for r in real_roots)


class TestEdgeCases:
    """Тесты граничных случаев"""

    def test_very_small_coefficients(self) -> None:
        """Тест с очень малыми коэффициентами"""
        roots = solve_quartic(0, 0, 1e-10, 0, -1e-10)
        assert isinstance(roots, list)
        assert len(roots) == 2

    def test_very_large_coefficients(self) -> None:
        """Тест с очень большими коэффициентами"""
        roots = solve_quartic(0, 0, 1e10, 0, -1e10)
        assert isinstance(roots, list)
        assert len(roots) == 2

    def test_negative_coefficients(self) -> None:
        """Тест с отрицательными коэффициентами"""
        roots = solve_quartic(0, 0, -1, 0, -4)
        assert isinstance(roots, list)
        assert len(roots) == 2
        # Должны быть комплексные корни
        assert all(isinstance(r, complex) for r in roots)

    def test_linear_negative_d(self) -> None:
        """Линейное уравнение с отрицательным d"""
        roots = solve_quartic(0, 0, 0, -3, 9)
        assert isinstance(roots, list)
        assert len(roots) == 1
        assert abs(roots[0] - 3.0) < 1e-9

    def test_quartic_zero_constant(self) -> None:
        """Уравнение 4-ой степени с k = 0"""
        roots = solve_quartic(1, 0, 0, 0, 0)
        assert isinstance(roots, list)
        assert len(roots) == 4
        # x = 0 является одним из корней
        assert any(abs(root) < 1e-7 for root in roots if isinstance(root, float) or abs(root.imag) < 1e-9)

    def test_quartic_negative_leading_coefficient(self) -> None:
        """Уравнение 4-ой степени с отрицательным старшим коэффициентом"""
        roots = solve_quartic(-1, 0, 0, 0, 16)
        assert isinstance(roots, list)
        assert len(roots) == 4
        # -x^4 + 16 = 0 => x^4 = 16 => x = ±2, ±2i
        # Проверяем, что все корни удовлетворяют уравнению
        for root in roots:
            value = -1 * root**4 + 16
            assert abs(value) < 1e-7

    def test_linear_with_zero_k(self) -> None:
        """Линейное уравнение dx = 0, корень x = 0"""
        roots = solve_quartic(0, 0, 0, 5, 0)
        assert isinstance(roots, list)
        assert len(roots) == 1
        assert abs(roots[0]) < 1e-9

    def test_quadratic_perfect_square(self) -> None:
        """Квадратное уравнение (x-3)^2 = x^2 - 6x + 9 = 0"""
        roots = solve_quartic(0, 0, 1, -6, 9)
        assert isinstance(roots, list)
        assert len(roots) == 2
        # Оба корня равны 3
        assert abs(roots[0] - 3.0) < 1e-9
        assert abs(roots[1] - 3.0) < 1e-9

    def test_cubic_with_negative_b(self) -> None:
        """Кубическое уравнение с отрицательным b"""
        roots = solve_quartic(0, -1, 0, 0, 8)
        assert isinstance(roots, list)
        assert len(roots) == 3
        # -x^3 + 8 = 0 => x^3 = 8 => x = 2
        real_roots = [r for r in roots if isinstance(r, float) or abs(r.imag) < 1e-9]
        complex_roots = [r for r in roots if isinstance(r, complex) and abs(r.imag) >= 1e-9]
        assert len(real_roots) == 1
        assert len(complex_roots) == 2
        assert any(abs(r - 2.0) < 1e-9 for r in real_roots)

    def test_quartic_with_mixed_coefficients(self) -> None:
        """Уравнение 4-ой степени со смешанными коэффициентами"""
        roots = solve_quartic(2, -3, 0, 5, -7)
        assert isinstance(roots, list)
        assert len(roots) == 4
        # Проверяем, что все найденные корни удовлетворяют уравнению
        for root in roots:
            value = 2 * root**4 - 3 * root**3 + 5 * root - 7
            assert abs(value) < 1e-5


class TestValidation:
    """Тесты для проверки корректности вычислений"""

    def test_verify_quartic_root(self) -> None:
        """Проверка, что все найденные корни действительно удовлетворяют уравнению"""
        a, b, c, d, k = 1, -2, -5, 6, 0
        roots = solve_quartic(a, b, c, d, k)
        assert len(roots) == 4
        # Проверяем, что все корни удовлетворяют уравнению
        for root in roots:
            value = a * root**4 + b * root**3 + c * root**2 + d * root + k
            assert abs(value) < 1e-5

    def test_verify_cubic_roots(self) -> None:
        """Проверка всех корней кубического уравнения"""
        roots = solve_quartic(0, 1, -3, 3, -1)  # (x-1)^3 = x^3 - 3x^2 + 3x - 1
        # Все три корня должны быть близки к 1
        for root in roots:
            if isinstance(root, complex):
                r = root.real if abs(root.imag) < 1e-9 else root
            else:
                r = root
            if isinstance(r, (int, float)):
                # Проверяем уравнение
                value = 1 * r**3 - 3 * r**2 + 3 * r - 1
                assert abs(value) < 1e-6

    def test_verify_quadratic_discriminant_positive(self) -> None:
        """Проверка квадратного уравнения с положительным дискриминантом"""
        roots = solve_quartic(0, 0, 1, -5, 6)  # x^2 - 5x + 6 = 0, корни: 2, 3
        assert len(roots) == 2
        # Проверяем оба корня
        for root in roots:
            r = root.real if isinstance(root, complex) else root
            value = r**2 - 5*r + 6
            assert abs(value) < 1e-9

    def test_verify_quadratic_discriminant_negative(self) -> None:
        """Проверка квадратного уравнения с отрицательным дискриминантом"""
        roots = solve_quartic(0, 0, 1, 2, 5)  # x^2 + 2x + 5 = 0
        print(f"\n  Найденные корни:")
        for i, root in enumerate(roots, 1):
            print(f"    x{i} = {format_root(root)}")
        assert len(roots) == 2
        assert all(isinstance(r, complex) for r in roots)
        # Проверяем комплексные корни
        for root in roots:
            value = root**2 + 2*root + 5
            assert abs(value) < 1e-9


if __name__ == "__main__":
    _ = pytest.main([__file__, "-v"])
