"""Lesson 4: 滤器元类、命名约定与真实筛选逻辑

OpenMC 的 `Filter` 体系不仅依赖元类保证命名一致，还定义了**如何根据粒子轨迹
映射到具体的 tally bin**。本课程实现了一个可运行的迷你滤器框架，覆盖：

* `FilterMeta`：确保所有滤器以 ``Filter`` 结尾，并自动生成 `short_name`，与 OpenMC
  中相同的 API 行为保持一致。
* `BaseFilter`：提供和 OpenMC `Filter` 基类类似的方法，如 `num_bins`、
  `bins`、`apply`，能够把事件映射成 ``(bin_index, weight)``。
* `EnergyFilter` 与 `CellFilter`：分别按照能量区间和几何单元 ID 进行匹配，真实地
  将粒子事件映射到 tally 维度。
* `FilterPipeline`：和 OpenMC `Tally._mean` 中的逻辑呼应，可把多个滤器联用，生成
  多维 bin 索引。

通过直接运行脚本可以看到与 OpenMC Python API 等价的行为：在离散化能量段、几何
单元的双重滤器下，事件被映射到组合 bin 索引，并给出权重。
"""

from __future__ import annotations

from bisect import bisect_right
from dataclasses import dataclass
from typing import Any, Iterable, Sequence


class FilterMeta(type):
    """要求子类名称以 ``Filter`` 结尾，并生成 ``short_name``。"""

    def __new__(mcls, name: str, bases: tuple[type, ...], namespace: dict[str, Any]):
        if not name.endswith("Filter") and name != "BaseFilter":
            raise TypeError("滤器类名称必须以 'Filter' 结尾")

        cls = super().__new__(mcls, name, bases, namespace)
        if name != "BaseFilter":
            cls.short_name = name.replace("Filter", "").lower()
        return cls


class BaseFilter(metaclass=FilterMeta):
    """和 OpenMC 中等价的滤器基类。"""

    short_name: str

    def __init__(self) -> None:
        self._num_bins: int = 0

    @property
    def num_bins(self) -> int:
        return self._num_bins

    def describe(self) -> str:
        return f"{self.__class__.__name__}(short_name={self.short_name}, bins={self.num_bins})"

    def apply(self, event: "TallyEvent") -> list[tuple[int, float]]:
        """返回 ``(bin_index, weight)`` 列表，仿照 OpenMC `Filter.apply`."""

        raise NotImplementedError


class EnergyFilter(BaseFilter):
    """能量滤器，按照能量区间映射粒子事件。"""

    def __init__(self, bins: Sequence[float]):
        super().__init__()
        if len(bins) < 2:
            raise ValueError("能量滤器至少需要两个边界")
        if sorted(bins) != list(bins):
            raise ValueError("能量边界必须递增")
        self._bins = list(bins)
        # OpenMC 在构造时记录 bin 数量
        self._num_bins = len(self._bins) - 1

    @property
    def bins(self) -> Sequence[float]:
        return self._bins

    def apply(self, event: "TallyEvent") -> list[tuple[int, float]]:
        idx = bisect_right(self._bins, event.energy) - 1
        if 0 <= idx < self._num_bins:
            return [(idx, 1.0)]
        return []


class CellFilter(BaseFilter):
    """根据几何 Cell ID 匹配事件。"""

    def __init__(self, cell_ids: Iterable[int]):
        super().__init__()
        cells = list(dict.fromkeys(cell_ids))
        if not cells:
            raise ValueError("CellFilter 至少需要一个 cell ID")
        self._cells = cells
        self._num_bins = len(self._cells)

    @property
    def cells(self) -> Sequence[int]:
        return self._cells

    def apply(self, event: "TallyEvent") -> list[tuple[int, float]]:
        try:
            idx = self._cells.index(event.cell_id)
        except ValueError:
            return []
        return [(idx, 1.0)]


@dataclass
class TallyEvent:
    """模拟一次粒子在 tally 中被采样的事件。"""

    energy: float
    cell_id: int
    weight: float = 1.0


class FilterPipeline:
    """组合多个滤器，产生和 OpenMC 中一致的多维 bin 索引。"""

    def __init__(self, filters: Sequence[BaseFilter]):
        self.filters = list(filters)

    def apply(self, event: TallyEvent) -> list[tuple[tuple[int, ...], float]]:
        indices = [f.apply(event) for f in self.filters]
        if any(not matches for matches in indices):
            return []

        combined: list[tuple[tuple[int, ...], float]] = [(tuple(), event.weight)]
        for matches in indices:
            new_combined: list[tuple[tuple[int, ...], float]] = []
            for base_bins, base_weight in combined:
                for idx, weight in matches:
                    new_combined.append((base_bins + (idx,), base_weight * weight))
            combined = new_combined
        return combined


if __name__ == "__main__":
    energy_filter = EnergyFilter(bins=[0.0, 1e3, 20e6])
    cell_filter = CellFilter(cell_ids=[1, 2, 3])

    event = TallyEvent(energy=800.0, cell_id=2, weight=0.5)
    pipeline = FilterPipeline([cell_filter, energy_filter])

    print(energy_filter.describe())
    print(cell_filter.describe())
    print("事件映射结果:", pipeline.apply(event))

    try:
        class BadName(BaseFilter):
            pass
    except TypeError as exc:
        print(f"创建失败: {exc}")
