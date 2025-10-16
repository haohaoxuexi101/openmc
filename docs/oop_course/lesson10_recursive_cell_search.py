"""Lesson 10: 几何递归搜索算法

本示例实现了一个精简版本的 ``Geometry.find_cell`` 算法，用以演示
OpenMC 在 CSG 几何树中定位粒子所属栅元的做法。核心思想是：

1. ``BoundingBox`` 先利用快速包围盒测试剔除不可能的栅元；
2. 若栅元填充了子 ``Universe``，则递归进入子树继续搜索；
3. 若命中终端栅元，则返回栅元以及访问路径，便于调试。

示例中我们构建了一个简单的二维模型：
``root`` 宇宙包含 ``fuel`` 和 ``moderator`` 两个栅元，``fuel`` 内又嵌套了
``pin`` 子宇宙（含燃料芯块和包壳）。算法将展示粒子位于不同坐标时的
搜索路径。
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, List, Optional, Tuple


@dataclass
class BoundingBox:
    """轴对齐包围盒，用于快速排除点不在栅元内的情况。"""

    xmin: float
    xmax: float
    ymin: float
    ymax: float

    def contains(self, point: Tuple[float, float]) -> bool:
        x, y = point
        return self.xmin <= x <= self.xmax and self.ymin <= y <= self.ymax


class Region:
    """代表一个布尔几何区域，封装 ``contains`` 逻辑。"""

    def __init__(self, predicate: Callable[[Tuple[float, float]], bool]):
        self._predicate = predicate

    def contains(self, point: Tuple[float, float]) -> bool:
        return self._predicate(point)


class Cell:
    """栅元，可填充另一宇宙，或为空（表示材料/终端单元）。"""

    def __init__(
        self,
        name: str,
        region: Region,
        bounding_box: BoundingBox,
        fill: Optional["Universe"] = None,
    ) -> None:
        self.name = name
        self.region = region
        self.bounding_box = bounding_box
        self.fill = fill

    def contains(self, point: Tuple[float, float]) -> bool:
        return self.bounding_box.contains(point) and self.region.contains(point)


class Universe:
    """宇宙，维护一组栅元并提供递归搜索能力。"""

    def __init__(self, name: str):
        self.name = name
        self.cells: List[Cell] = []

    def add_cell(self, cell: Cell) -> None:
        self.cells.append(cell)

    def find_cell(self, point: Tuple[float, float], path: Optional[List[str]] = None) -> Tuple[Cell, List[str]]:
        """递归查找包含 ``point`` 的终端栅元及访问路径。"""

        if path is None:
            path = [self.name]
        else:
            path.append(self.name)

        for cell in self.cells:
            if not cell.contains(point):
                continue

            path.append(cell.name)
            if cell.fill is None:
                return cell, path
            return cell.fill.find_cell(point, path)

        raise LookupError(f"Point {point} not found in universe '{self.name}'")


# --- 构建示例几何 ---------------------------------------------------------

def circle_region(radius: float) -> Region:
    return Region(lambda p: p[0] ** 2 + p[1] ** 2 <= radius**2)


def annulus_region(r_inner: float, r_outer: float) -> Region:
    return Region(
        lambda p: r_inner**2 <= p[0] ** 2 + p[1] ** 2 <= r_outer**2
    )


pin = Universe("pin")
pin.add_cell(
    Cell(
        name="fuel pellet",
        region=circle_region(0.4),
        bounding_box=BoundingBox(-0.4, 0.4, -0.4, 0.4),
    )
)
pin.add_cell(
    Cell(
        name="cladding",
        region=annulus_region(0.4, 0.45),
        bounding_box=BoundingBox(-0.45, 0.45, -0.45, 0.45),
    )
)

root = Universe("root")
root.add_cell(
    Cell(
        name="fuel",  # 填充子宇宙
        region=circle_region(0.45),
        bounding_box=BoundingBox(-0.45, 0.45, -0.45, 0.45),
        fill=pin,
    )
)
root.add_cell(
    Cell(
        name="moderator",
        region=Region(lambda p: p[0] ** 2 + p[1] ** 2 > 0.45**2),
        bounding_box=BoundingBox(-1.0, 1.0, -1.0, 1.0),
    )
)

# --- 运行演示 -------------------------------------------------------------

POINTS = [(0.1, 0.1), (0.43, 0.0), (0.8, 0.0)]

for pt in POINTS:
    try:
        cell, search_path = root.find_cell(pt)
    except LookupError as err:
        print(err)
    else:
        print(f"Point {pt} located in cell '{cell.name}' via path -> {' > '.join(search_path)}")
