"""Lesson 2: 对象相等性与数组比较

OpenMC 的 `EqualityMixin` 统一了不同对象的相等比较逻辑，会自动处理 `numpy.ndarray`
等不可直接比较的成员。本示例提炼该思路，让你了解如何在数据类之间复用相等性判断。
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np


class EqualityMixin:
    """提供对数值和数组友好的 ``__eq__`` 实现。"""

    def __eq__(self, other: Any) -> bool:  # pragma: no cover - 演示代码
        if self.__class__ is not other.__class__:
            return NotImplemented

        for name, value in self.__dict__.items():
            other_value = getattr(other, name)
            if isinstance(value, np.ndarray):
                if not np.array_equal(value, other_value):
                    return False
            else:
                if value != other_value:
                    return False
        return True


@dataclass
class Spectrum(EqualityMixin):
    """计数器能谱，包含数组属性。"""

    groups: np.ndarray
    values: np.ndarray


@dataclass
class MaterialVector(EqualityMixin):
    """材料向量，混合标量和数组。"""

    name: str
    fractions: np.ndarray


if __name__ == "__main__":
    group_bounds = np.array([0.0, 1.0, 10.0])
    s1 = Spectrum(groups=group_bounds, values=np.array([0.2, 0.8]))
    s2 = Spectrum(groups=group_bounds.copy(), values=np.array([0.2, 0.8]))
    s3 = Spectrum(groups=group_bounds, values=np.array([0.1, 0.9]))

    print(f"s1 == s2 ? {s1 == s2}")
    print(f"s1 == s3 ? {s1 == s3}")

    mv1 = MaterialVector(name="fuel", fractions=np.array([0.97, 0.03]))
    mv2 = MaterialVector(name="fuel", fractions=np.array([0.97, 0.03]))
    mv3 = MaterialVector(name="moderator", fractions=np.array([0.97, 0.03]))

    print(f"mv1 == mv2 ? {mv1 == mv2}")
    print(f"mv1 == mv3 ? {mv1 == mv3}")
