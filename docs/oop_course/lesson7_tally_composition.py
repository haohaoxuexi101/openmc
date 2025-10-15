"""Lesson 7: 计数器组合操作

OpenMC 的 `Tally` 提供了滤器、得分、核素等列表的组合操作，还支持加减乘除等算术。
本示例使用 Lesson 6 中的 `CheckedList` 构建一个精简版计数器，展示如何在对象上定义
领域特定操作，让用户以自然的方式组合计算。"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, List

from lesson6_checked_list import CheckedList, Score


@dataclass
class Tally:
    name: str
    filters: CheckedList[str]
    scores: CheckedList[Score]
    nuclides: CheckedList[str]
    values: List[float]

    def __post_init__(self) -> None:
        if len(self.scores) != len(self.values):
            raise ValueError("得分数量需与值数量一致")

    def combine(self, other: "Tally", weight: float = 1.0) -> "Tally":
        """模仿 OpenMC 的加权相加。"""

        if self.scores._items != other.scores._items:  # 直接访问演示
            raise ValueError("仅支持相同得分的组合")

        new_values = [a + weight * b for a, b in zip(self.values, other.values)]
        return Tally(
            name=f"{self.name}+{other.name}",
            filters=self.filters,
            scores=self.scores,
            nuclides=self.nuclides,
            values=new_values,
        )

    def scale(self, factor: float) -> "Tally":
        return Tally(
            name=f"{self.name}*{factor}",
            filters=self.filters,
            scores=self.scores,
            nuclides=self.nuclides,
            values=[v * factor for v in self.values],
        )

    def describe(self) -> str:
        parts = [f"Tally({self.name})"]
        parts.append(f"  Filters: {list(self.filters)}")
        parts.append(f"  Scores: {[score.name for score in self.scores]}")
        parts.append(f"  Nuclides: {list(self.nuclides)}")
        parts.append(f"  Values: {self.values}")
        return "\n".join(parts)


def make_checked_strings(items: Iterable[str]) -> CheckedList[str]:
    return CheckedList(str, list(items))


if __name__ == "__main__":
    filters = make_checked_strings(["cell 1"])
    scores = CheckedList(Score, [Score("flux"), Score("fission")])
    nuclides = make_checked_strings(["U235", "U238"])

    tally1 = Tally(name="tally-1", filters=filters, scores=scores, nuclides=nuclides, values=[1.0, 2.0])
    tally2 = Tally(name="tally-2", filters=filters, scores=scores, nuclides=nuclides, values=[0.5, 0.7])

    combined = tally1.combine(tally2, weight=2.0)
    scaled = tally1.scale(10.0)

    print(tally1.describe())
    print(combined.describe())
    print(scaled.describe())
