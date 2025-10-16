"""Lesson 3: 几何抽象基类

OpenMC 使用 `UniverseBase` 等抽象基类统一几何构件的接口。本示例展示如何利用
`abc` 模块定义抽象方法，并在派生类中实现具体逻辑，确保扩展时接口一致。
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Dict, Iterable, List


class UniverseBase(ABC):
    """抽象宇宙基类，定义公共接口。"""

    @abstractmethod
    def add_cell(self, cell: "Cell") -> None:
        raise NotImplementedError

    @abstractmethod
    def cells(self) -> Iterable["Cell"]:
        raise NotImplementedError

    @abstractmethod
    def describe(self, indent: int = 0) -> str:
        """返回层次化描述。"""


@dataclass
class Cell:
    """简化的几何单元。"""

    name: str
    fill: str

    def describe(self) -> str:
        return f"Cell(name={self.name}, fill={self.fill})"


@dataclass
class SimpleUniverse(UniverseBase):
    """具体实现：线性存储所有 Cell。"""

    name: str
    _cells: List[Cell] = field(default_factory=list)

    def add_cell(self, cell: Cell) -> None:
        print(f"[SimpleUniverse] Adding cell {cell.name}")
        self._cells.append(cell)

    def cells(self) -> Iterable[Cell]:
        return iter(self._cells)

    def describe(self, indent: int = 0) -> str:
        spaces = " " * indent
        description = [f"{spaces}Universe({self.name})"]
        for cell in self._cells:
            description.append(f"{spaces}  - {cell.describe()}")
        return "\n".join(description)


@dataclass
class LatticeUniverse(UniverseBase):
    """具体实现：网格化存储 Cell，模仿 OpenMC `Lattice`。"""

    name: str
    width: int
    height: int
    _grid: Dict[tuple[int, int], Cell] = field(default_factory=dict)

    def add_cell(self, cell: Cell, position: tuple[int, int]) -> None:
        if position in self._grid:
            raise ValueError(f"位置 {position} 已有 Cell")
        print(f"[LatticeUniverse] Placing cell {cell.name} at {position}")
        self._grid[position] = cell

    # 为了满足抽象接口，提供不带 position 的方法代理
    def cells(self) -> Iterable[Cell]:
        return self._grid.values()

    def describe(self, indent: int = 0) -> str:
        spaces = " " * indent
        description = [f"{spaces}LatticeUniverse({self.name})"]
        for (i, j), cell in sorted(self._grid.items()):
            description.append(f"{spaces}  - ({i}, {j}) => {cell.describe()}")
        return "\n".join(description)


if __name__ == "__main__":
    cell_fuel = Cell(name="fuel", fill="UO2")
    cell_moderator = Cell(name="moderator", fill="Water")

    root_universe = SimpleUniverse(name="root")
    root_universe.add_cell(cell_fuel)
    root_universe.add_cell(cell_moderator)

    lattice = LatticeUniverse(name="assembly", width=2, height=2)
    lattice.add_cell(cell_fuel, position=(0, 0))
    lattice.add_cell(cell_moderator, position=(0, 1))

    print(root_universe.describe())
    print(lattice.describe())
