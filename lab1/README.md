# Лабораторная работа №1: Решение уравнения 4-ой степени

## Описание задачи

Необходимо разработать программу решения уравнения 4-ой степени:

```
ax⁴ + bx³ + cx² + dx + k = 0
```

### Требования к реализации

1. **При a ≠ 0**: Реализовать **метод касательных (Ньютона)** для получения первого корня уравнения
2. **При a = 0, b ≠ 0**: Реализовать **метод Кардано** для кубического уравнения
3. **При a = 0, b = 0, c ≠ 0**: Решать квадратное уравнение по **формуле Виета**
4. **При a = 0, b = 0, c = 0, d ≠ 0**: Решать линейное уравнение: `x = -k / d`
5. **Если все коэффициенты = 0**:
   - При k = 0: любое число является решением
   - При k ≠ 0: решений нет

---

## Структура проекта

```
lab1/
├── main.py                    # Основная программа
├── tests.py                   # Тестовый набор (35 тестов)
├── algorithm_flowchart.drawio # Блок-схема алгоритма
└── README.md                  # Документация
```

---

## Описание функций

### 1. `newton_method(a, b, c, d, k, x0=1, tolerance=1e-10, max_iterations=1000)`

**Метод Ньютона (касательных)** для нахождения корня уравнения 4-ой степени.

**Формулы:**
- f(x) = ax⁴ + bx³ + cx² + dx + k
- f'(x) = 4ax³ + 3bx² + 2cx + d
- x_{n+1} = x_n - f(x_n) / f'(x_n)

**Особенности реализации:**
- Пробует несколько начальных точек для обеспечения сходимости
- Возвращает наилучшее найденное приближение
- Обрабатывает случаи, когда производная близка к нулю

### 2. `cardano_method(b, c, d, k)`

**Метод Кардано** для решения кубического уравнения bx³ + cx² + dx + k = 0.

**Алгоритм:**
1. Приведение к виду x³ + px + q = 0
2. Вычисление дискриминанта: D = (q/2)² + (p/3)³
3. В зависимости от знака D:
   - D > 0: один действительный и два комплексных корня
   - D = 0: кратные корни
   - D < 0: три различных действительных корня

**Возвращает:** список из 3 корней (действительных или комплексных)

### 3. `vieta_quadratic(c, d, k)`

Решение квадратного уравнения cx² + dx + k = 0 **по формуле Виета**.

**Формулы:**
- D = d² - 4ck
- x₁ = (-d + √D) / (2c)
- x₂ = (-d - √D) / (2c)

**Возвращает:** список из 2 корней (действительных или комплексных)

### 4. `solve_quartic(a, b, c, d, k)`

**Главная функция**, которая определяет тип уравнения и вызывает соответствующий метод решения.

**Возвращает:**
- Список корней (float или complex)
- Строку "Любое число" (если все коэф. = 0 и k = 0)
- Пустой список (если решений нет)

---

## Тестовый набор

### Всего тестов: **35**

#### 1. TestNewtonMethod (4 теста)
- `test_simple_quartic` - простое уравнение x⁴ - 1 = 0
- `test_quartic_with_all_coefficients` - уравнение с различными коэффициентами
- `test_quartic_x4_minus_16` - уравнение x⁴ - 16 = 0
- `test_quartic_negative_start` - тест с отрицательной начальной точкой

#### 2. TestCardanoMethod (4 теста)
- `test_simple_cubic` - x³ - 1 = 0
- `test_cubic_three_real_roots` - три действительных корня
- `test_cubic_x3_minus_8` - x³ - 8 = 0
- `test_cubic_with_complex_roots` - комплексные корни

#### 3. TestVietaQuadratic (5 тестов)
- `test_simple_quadratic` - x² - 4 = 0
- `test_quadratic_x2_minus_1` - x² - 1 = 0
- `test_quadratic_complex_roots` - комплексные корни
- `test_quadratic_with_all_coefficients` - x² + 3x + 2 = 0
- `test_quadratic_double_root` - кратный корень

#### 4. TestSolveQuartic (8 тестов)
Интеграционные тесты для всех случаев:
- Уравнение 4-ой степени (a ≠ 0)
- Кубическое уравнение (b ≠ 0)
- Квадратное уравнение (c ≠ 0)
- Линейное уравнение (d ≠ 0)
- Все коэффициенты = 0 (два случая)
- Комплексные корни

#### 5. TestEdgeCases (10 тестов)
Граничные случаи:
- Очень малые коэффициенты (1e-10)
- Очень большие коэффициенты (1e10)
- Отрицательные коэффициенты
- Отрицательный старший коэффициент
- Нулевой свободный член
- Полный квадрат
- Смешанные коэффициенты

#### 6. TestValidation (4 теста)
Проверка корректности найденных корней:
- Подстановка корня в исходное уравнение
- Проверка всех корней кубического уравнения
- Квадратное уравнение с положительным дискриминантом
- Квадратное уравнение с отрицательным дискриминантом

---

## Запуск программы

### Интерактивный режим:

```bash
cd lab1
python main.py
```

Программа запросит ввод коэффициентов a, b, c, d, k.

### Запуск тестов:

```bash
# Все тесты
python -m pytest tests.py -v

# Конкретный класс тестов
python -m pytest tests.py::TestNewtonMethod -v

# Конкретный тест
python -m pytest tests.py::TestNewtonMethod::test_simple_quartic -v
```

---

## Примеры использования

### Пример 1: Уравнение 4-ой степени

```python
# x⁴ - 16 = 0, корни: ±2, ±2i
roots = solve_quartic(1, 0, 0, 0, -16)
# Результат: [2.0] (первый действительный корень)
```

### Пример 2: Кубическое уравнение

```python
# x³ - 6x² + 11x - 6 = 0, корни: 1, 2, 3
roots = solve_quartic(0, 1, -6, 11, -6)
# Результат: [1.0, 2.0, 3.0]
```

### Пример 3: Квадратное уравнение

```python
# x² - 4 = 0, корни: -2, 2
roots = solve_quartic(0, 0, 1, 0, -4)
# Результат: [-2.0, 2.0]
```

### Пример 4: Линейное уравнение

```python
# 2x - 10 = 0, корень: 5
roots = solve_quartic(0, 0, 0, 2, -10)
# Результат: [5.0]
```

### Пример 5: Вырожденные случаи

```python
# Любое число: 0 = 0
result = solve_quartic(0, 0, 0, 0, 0)
# Результат: "Любое число"

# Нет решений: 5 = 0
result = solve_quartic(0, 0, 0, 0, 5)
# Результат: []
```

---

## Блок-схема алгоритма

Визуальное представление алгоритма доступно в файле **`algorithm_flowchart.drawio`**.

Его можно открыть с помощью:
- [draw.io](https://app.diagrams.net/) (онлайн)
- [VS Code с расширением Draw.io Integration](https://marketplace.visualstudio.com/items?itemName=hediet.vscode-drawio)

### Структура блок-схемы:

```
Начало
   ↓
Ввод параметров (a, b, c, d, k)
   ↓
a ≠ 0? ──[Да]──→ Метод Ньютона ──→ Результат: x₁
   ↓[Нет]
b ≠ 0? ──[Да]──→ Метод Кардано ──→ Результат: x₁, x₂, x₃
   ↓[Нет]
c ≠ 0? ──[Да]──→ Формула Виета ──→ Результат: x₁, x₂
   ↓[Нет]
d ≠ 0? ──[Да]──→ Линейное уравнение ──→ Результат: x
   ↓[Нет]
k = 0? ──[Да]──→ Любое число
   ↓[Нет]
Решений нет
   ↓
Вывод результата
   ↓
Конец
```

---

## Результаты тестирования

```
============================= test session starts ==============================
collected 35 items

lab1/tests.py::TestNewtonMethod::test_simple_quartic PASSED              [  2%]
lab1/tests.py::TestNewtonMethod::test_quartic_with_all_coefficients PASSED [  5%]
lab1/tests.py::TestNewtonMethod::test_quartic_x4_minus_16 PASSED         [  8%]
lab1/tests.py::TestNewtonMethod::test_quartic_negative_start PASSED      [ 11%]
lab1/tests.py::TestCardanoMethod::test_simple_cubic PASSED               [ 14%]
lab1/tests.py::TestCardanoMethod::test_cubic_three_real_roots PASSED     [ 17%]
lab1/tests.py::TestCardanoMethod::test_cubic_x3_minus_8 PASSED           [ 20%]
lab1/tests.py::TestCardanoMethod::test_cubic_with_complex_roots PASSED   [ 22%]
lab1/tests.py::TestVietaQuadratic::test_simple_quadratic PASSED          [ 25%]
lab1/tests.py::TestVietaQuadratic::test_quadratic_x2_minus_1 PASSED      [ 28%]
lab1/tests.py::TestVietaQuadratic::test_quadratic_complex_roots PASSED   [ 31%]
lab1/tests.py::TestVietaQuadratic::test_quadratic_with_all_coefficients PASSED [ 34%]
lab1/tests.py::TestVietaQuadratic::test_quadratic_double_root PASSED     [ 37%]
lab1/tests.py::TestSolveQuartic::test_case_quartic_a_not_zero PASSED     [ 40%]
lab1/tests.py::TestSolveQuartic::test_case_cubic_b_not_zero PASSED       [ 42%]
lab1/tests.py::TestSolveQuartic::test_case_quadratic_c_not_zero PASSED   [ 45%]
lab1/tests.py::TestSolveQuartic::test_case_linear_d_not_zero PASSED      [ 48%]
lab1/tests.py::TestSolveQuartic::test_case_all_zero_k_zero PASSED        [ 51%]
lab1/tests.py::TestSolveQuartic::test_case_all_zero_k_not_zero PASSED    [ 54%]
lab1/tests.py::TestSolveQuartic::test_quadratic_complex_roots_full PASSED [ 57%]
lab1/tests.py::TestSolveQuartic::test_cubic_complex_coefficients PASSED  [ 60%]
lab1/tests.py::TestEdgeCases::test_very_small_coefficients PASSED        [ 62%]
lab1/tests.py::TestEdgeCases::test_very_large_coefficients PASSED        [ 65%]
lab1/tests.py::TestEdgeCases::test_negative_coefficients PASSED          [ 68%]
lab1/tests.py::TestEdgeCases::test_linear_negative_d PASSED              [ 71%]
lab1/tests.py::TestEdgeCases::test_quartic_zero_constant PASSED          [ 74%]
lab1/tests.py::TestEdgeCases::test_quartic_negative_leading_coefficient PASSED [ 77%]
lab1/tests.py::TestEdgeCases::test_linear_with_zero_k PASSED             [ 80%]
lab1/tests.py::TestEdgeCases::test_quadratic_perfect_square PASSED       [ 82%]
lab1/tests.py::TestEdgeCases::test_cubic_with_negative_b PASSED          [ 85%]
lab1/tests.py::TestEdgeCases::test_quartic_with_mixed_coefficients PASSED [ 88%]
lab1/tests.py::TestValidation::test_verify_quartic_root PASSED           [ 91%]
lab1/tests.py::TestValidation::test_verify_cubic_roots PASSED            [ 94%]
lab1/tests.py::TestValidation::test_verify_quadratic_discriminant_positive PASSED [ 97%]
lab1/tests.py::TestValidation::test_verify_quadratic_discriminant_negative PASSED [100%]

============================== 35 passed in 0.05s ========================
```

**✅ Все 35 тестов успешно пройдены!**

---

## Особенности реализации

### 1. Метод Ньютона
- Использует множественные начальные точки для повышения надежности
- Обрабатывает случаи, когда производная близка к нулю
- Возвращает наилучшее приближение из всех попыток

### 2. Метод Кардано
- Полная реализация для всех случаев дискриминанта
- Корректная обработка кратных корней
- Правильное вычисление комплексных корней

### 3. Формула Виета
- Использует комплексные числа для отрицательного дискриминанта
- Точное вычисление корней

### 4. Обработка особых случаев
- Вырожденные уравнения (все коэффициенты = 0)
- Очень малые и очень большие коэффициенты
- Отрицательные коэффициенты

---

## Требования

- Python 3.9+
- pytest 8.4.2+
- Стандартные библиотеки: math, cmath

## Автор

Лабораторная работа выполнена в рамках курса численных методов.

