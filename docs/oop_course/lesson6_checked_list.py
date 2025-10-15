"""Lesson 6: 类型安全容器

OpenMC 定义了 `CheckedList` 来确保列表元素始终满足类型要求。本示例展示如何通过
自定义容器类在添加、扩展时执行校验，从而在 API 层构建数据安全网。
"""

from __future__ import annotations

from typing import Generic, Iterable, Iterator, List, Sequence, TypeVar

T = TypeVar("T")


class CheckedList(Generic[T]):
    """限制元素类型的列表。"""

    def __init__(self, item_type: type[T], data: Sequence[T] | None = None):
        self.item_type = item_type
        self._items: List[T] = []
        if data is not None:
            self.extend(data)

    def _check(self, value: T) -> None:
        if not isinstance(value, self.item_type):
            raise TypeError(f"元素必须是 {self.item_type.__name__}，当前为 {type(value).__name__}")

    def append(self, value: T) -> None:
        self._check(value)
        self._items.append(value)

    def extend(self, values: Iterable[T]) -> None:
        for value in values:
            self.append(value)

    def insert(self, index: int, value: T) -> None:
        self._check(value)
        self._items.insert(index, value)

    def __iter__(self) -> Iterator[T]:
        return iter(self._items)

    def __len__(self) -> int:
        return len(self._items)

    def __getitem__(self, item):
        return self._items[item]

    def __repr__(self) -> str:
        return f"CheckedList(item_type={self.item_type.__name__}, items={self._items})"


class Score:
    """计数器得分，作为容器的元素类型。"""

    def __init__(self, name: str):
        self.name = name

    def __repr__(self) -> str:  # pragma: no cover - 演示
        return f"Score({self.name})"


if __name__ == "__main__":
    scores = CheckedList(Score)
    scores.append(Score("flux"))
    scores.extend([Score("fission"), Score("absorption")])
    print(scores)

    try:
        scores.append("not-a-score")  # type: ignore[arg-type]
    except TypeError as exc:
        print(f"类型校验失败: {exc}")
