from __future__ import annotations

import os
import re
import subprocess
import sys
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]


def matvec(A: list[list[float]], x: list[float]) -> list[float]:
    return [sum(A[i][j] * x[j] for j in range(len(x))) for i in range(len(A))]


def run_main(tmp_path: Path, *, a_txt: str, b_txt: str) -> subprocess.CompletedProcess[str]:
    (tmp_path / "A.txt").write_text(a_txt, encoding="utf-8")
    (tmp_path / "B.txt").write_text(b_txt, encoding="utf-8")

    env = os.environ.copy()
    # Ensure `python -m lab2.main` can import the local package even when cwd=tmp_path
    env["PYTHONPATH"] = str(REPO_ROOT) + os.pathsep + env.get("PYTHONPATH", "")

    return subprocess.run(
        [sys.executable, "-m", "lab2.main"],
        cwd=str(tmp_path),
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )


def parse_stdout(stdout: str) -> tuple[list[float], int]:
    # Expected lines:
    #   x1 = 0.1429
    #   ...
    #   Итераций: 22
    xs: dict[int, float] = {}
    for m in re.finditer(
        r"^x(?P<i>\d+)\s*=\s*(?P<v>[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)\s*$",
        stdout,
        flags=re.MULTILINE,
    ):
        xs[int(m.group("i"))] = float(m.group("v"))

    it_m = re.search(r"Итераций:\s*(\d+)", stdout)
    if it_m is None:
        raise AssertionError(f"Не удалось найти строку 'Итераций: ...' в выводе:\n{stdout}")

    if not xs:
        raise AssertionError(f"Не удалось найти строки вида 'xN = ...' в выводе:\n{stdout}")

    n = max(xs.keys())
    x_list = [xs[i] for i in range(1, n + 1)]
    return x_list, int(it_m.group(1))


class TestMainProgram:
    def test_success_identity_2x2(self, tmp_path: Path) -> None:
        res = run_main(tmp_path, a_txt="1 0\n0 1\n", b_txt="7\n-3\n")
        assert res.returncode == 0
        x, it = parse_stdout(res.stdout)
        assert x == [7.0, -3.0]
        assert it >= 0

    def test_success_known_2x2(self, tmp_path: Path) -> None:
        # 4x + y = 9
        # 2x + 3y = 13
        # x=1.4, y=3.4
        res = run_main(tmp_path, a_txt="4 1\n2 3\n", b_txt="9\n13\n")
        assert res.returncode == 0
        x, it = parse_stdout(res.stdout)
        assert abs(x[0] - 1.4) < 1e-6
        assert abs(x[1] - 3.4) < 1e-6
        assert it > 0

    def test_error_non_square_matrix(self, tmp_path: Path) -> None:
        res = run_main(tmp_path, a_txt="1 2 3\n4 5 6\n", b_txt="1\n2\n")
        assert res.returncode == 1
        assert "Матрица A не является квадратной" in res.stdout

    def test_error_size_mismatch(self, tmp_path: Path) -> None:
        res = run_main(tmp_path, a_txt="1 0\n0 1\n", b_txt="1\n2\n3\n")
        assert res.returncode == 1
        assert "Размерность b" in res.stdout

    def test_error_zero_diagonal(self, tmp_path: Path) -> None:
        res = run_main(tmp_path, a_txt="0 1\n1 1\n", b_txt="1\n2\n")
        assert res.returncode == 1
        assert "Диагональный элемент A[0][0] равен 0." in res.stdout

    def test_error_not_diagonally_dominant(self, tmp_path: Path) -> None:
        res = run_main(tmp_path, a_txt="1 2\n2 1\n", b_txt="1\n1\n")
        assert res.returncode == 1
        assert "Матрица A не удовлетворяет условию диагонального преобладания" in res.stdout


@pytest.mark.parametrize(
    "A,x_true",
    [
        # 3x3 (5 шт.)
        (
            [[10.0, 1.0, 1.0], [2.0, 10.0, 1.0], [2.0, 2.0, 10.0]],
            [1.0, 2.0, 3.0],
        ),
        (
            [[8.0, -1.0, 0.5], [1.0, 9.0, -1.0], [0.25, 1.0, 7.0]],
            [-2.0, 1.0, 0.5],
        ),
        (
            [[12.0, 0.5, -0.5], [-1.0, 11.0, 1.0], [0.5, -1.0, 9.0]],
            [3.0, -1.0, 2.0],
        ),
        (
            [[6.0, 1.0, 0.0], [1.0, 7.0, 1.0], [0.0, 1.0, 8.0]],
            [1.5, -2.0, 0.25],
        ),
        (
            [[15.0, 1.0, 2.0], [1.0, 14.0, -1.0], [2.0, -1.0, 13.0]],
            [0.0, 1.0, -1.0],
        ),
        # 4x4 (5 шт.)
        (
            [
                [20.0, 1.0, 0.0, -1.0],
                [2.0, 18.0, 1.0, 0.0],
                [0.0, -1.0, 16.0, 2.0],
                [1.0, 0.0, 1.0, 17.0],
            ],
            [1.0, -2.0, 0.5, 3.0],
        ),
        (
            [
                [12.0, -1.0, 1.0, 0.5],
                [1.0, 13.0, -1.0, 1.0],
                [0.5, 1.0, 14.0, -1.0],
                [-1.0, 0.5, 1.0, 15.0],
            ],
            [-1.0, 2.0, -3.0, 0.25],
        ),
        (
            [
                [9.0, 1.0, 1.0, 0.0],
                [1.0, 10.0, 0.0, 1.0],
                [1.0, 0.0, 11.0, 1.0],
                [0.0, 1.0, 1.0, 12.0],
            ],
            [2.0, 1.0, -1.0, 0.5],
        ),
        (
            [
                [25.0, 2.0, -1.0, 0.0],
                [1.0, 22.0, 2.0, -1.0],
                [-1.0, 1.0, 24.0, 2.0],
                [0.0, -1.0, 1.0, 23.0],
            ],
            [0.0, 0.0, 1.0, -1.0],
        ),
        (
            [
                [14.0, 1.0, -1.0, 1.0],
                [2.0, 15.0, 1.0, -1.0],
                [-1.0, 1.0, 16.0, 1.0],
                [1.0, -1.0, 2.0, 17.0],
            ],
            [3.0, -2.0, 1.0, -0.5],
        ),
        # 5x5 (5 шт.)
        (
            [
                [30.0, 1.0, 0.0, -1.0, 0.5],
                [2.0, 28.0, 1.0, 0.0, -1.0],
                [0.0, -1.0, 27.0, 2.0, 1.0],
                [1.0, 0.0, 1.0, 29.0, -1.0],
                [-0.5, 1.0, 0.0, 2.0, 26.0],
            ],
            [1.0, -2.0, 0.5, 3.0, -1.0],
        ),
        (
            [
                [18.0, 1.0, -1.0, 0.0, 0.5],
                [1.0, 19.0, 1.0, -1.0, 0.0],
                [-1.0, 1.0, 20.0, 1.0, -1.0],
                [0.0, -1.0, 1.0, 21.0, 1.0],
                [0.5, 0.0, -1.0, 1.0, 22.0],
            ],
            [-1.0, 2.0, -3.0, 0.25, 1.5],
        ),
        (
            [
                [16.0, 1.0, 1.0, 0.0, 0.0],
                [1.0, 17.0, 0.0, 1.0, 0.0],
                [1.0, 0.0, 18.0, 1.0, 0.0],
                [0.0, 1.0, 1.0, 19.0, 1.0],
                [0.0, 0.0, 0.0, 1.0, 20.0],
            ],
            [2.0, 1.0, -1.0, 0.5, -2.5],
        ),
        (
            [
                [40.0, -2.0, 1.0, 0.0, -1.0],
                [1.0, 38.0, -2.0, 1.0, 0.0],
                [-1.0, 1.0, 39.0, -2.0, 1.0],
                [0.0, -1.0, 1.0, 37.0, -2.0],
                [2.0, 0.0, -1.0, 1.0, 36.0],
            ],
            [0.0, 0.0, 1.0, -1.0, 2.0],
        ),
        (
            [
                [22.0, 1.0, 0.0, -1.0, 1.0],
                [2.0, 23.0, 1.0, 0.0, -1.0],
                [0.0, 2.0, 24.0, 1.0, 0.0],
                [-1.0, 0.0, 2.0, 25.0, 1.0],
                [1.0, -1.0, 0.0, 2.0, 26.0],
            ],
            [3.0, -2.0, 1.0, -0.5, 0.75],
        ),
    ],
)
def test_main_program_compares_roots(tmp_path: Path, A: list[list[float]], x_true: list[float]) -> None:
    b = matvec(A, x_true)
    a_txt = "\n".join(" ".join(str(v) for v in row) for row in A) + "\n"
    b_txt = "\n".join(str(v) for v in b) + "\n"

    res = run_main(tmp_path, a_txt=a_txt, b_txt=b_txt)
    assert res.returncode == 0, f"stdout:\n{res.stdout}\nstderr:\n{res.stderr}"

    x, _ = parse_stdout(res.stdout)
    assert len(x) == len(x_true)
    # main.py печатает значения с форматом {:.4g}, поэтому точность ограничена выводом
    for got, exp in zip(x, x_true):
        assert abs(got - exp) <= 2e-3


