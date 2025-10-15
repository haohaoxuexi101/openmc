"""Lesson 12: Log-Log 插值算法

OpenMC 处理中子截面时大量使用对数-对数插值（log-log interpolation），
以保持在数量级跨越巨大时的精度。本课程构建 ``LogLogInterpolator``，
演示其在截面缩放、温度表格外推中的应用，重点关注以下算法细节：

* 使用 `numpy.log` 将原始数据转换到对数空间，转换后再执行线性插值；
* 处理表格外 extrapolation 时，保持与 OpenMC 一致的幂律外推；
* 支持批量向量输入，以匹配 OpenMC 对矢量化运算的追求。

示例中我们以一段理想化的 `1/v` 吸收截面为例，对比线性插值和 log-log
插值的误差差异。
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Union

import numpy as np


@dataclass
class LogLogInterpolator:
    energies: np.ndarray
    values: np.ndarray

    def __post_init__(self) -> None:
        if np.any(self.energies <= 0) or np.any(self.values <= 0):
            raise ValueError("log-log interpolation requires positive grids and values")
        order = np.argsort(self.energies)
        self.energies = self.energies[order]
        self.values = self.values[order]
        self._log_e = np.log(self.energies)
        self._log_v = np.log(self.values)

    def __call__(self, energy: Union[Iterable[float], float]) -> np.ndarray:
        x = np.atleast_1d(energy).astype(float)
        log_x = np.log(x)
        log_y = np.interp(log_x, self._log_e, self._log_v, left=None, right=None)

        below = log_x < self._log_e[0]
        if np.any(below):
            slope = (self._log_v[1] - self._log_v[0]) / (self._log_e[1] - self._log_e[0])
            log_y[below] = self._log_v[0] + slope * (log_x[below] - self._log_e[0])

        above = log_x > self._log_e[-1]
        if np.any(above):
            slope = (
                self._log_v[-1] - self._log_v[-2]
            ) / (self._log_e[-1] - self._log_e[-2])
            log_y[above] = self._log_v[-1] + slope * (log_x[above] - self._log_e[-1])

        return np.exp(log_y)


def linear_interp(energies: Iterable[float], xs: Iterable[float], energy: Iterable[float]) -> np.ndarray:
    return np.interp(energy, energies, xs)


if __name__ == "__main__":
    grid = np.array([1e-5, 1e-3, 1e-1, 1e1])
    xs = 1.0 / np.sqrt(grid)  # 理想化的 1/v 截面

    log_interp = LogLogInterpolator(grid, xs)

    sample = np.logspace(-6, 2, num=9)
    linear_values = linear_interp(grid, xs, sample)
    loglog_values = log_interp(sample)

    true_values = 1.0 / np.sqrt(sample)

    linear_error = np.abs(linear_values - true_values)
    loglog_error = np.abs(loglog_values - true_values)

    print("energy (eV) | true    | linear  | log-log | linear err | log-log err")
    for e, t, l, g, el, eg in zip(sample, true_values, linear_values, loglog_values, linear_error, loglog_error):
        print(f"{e:11.4e} | {t:7.3f} | {l:7.3f} | {g:7.3f} | {el:10.3e} | {eg:11.3e}")
