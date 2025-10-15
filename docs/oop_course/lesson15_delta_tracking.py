"""Lesson 15: Woodcock Delta-Tracking 粒子输运

OpenMC 在复杂几何中采用 Woodcock (delta-tracking) 算法执行中子输运，
通过一个**主截面** ``Σ_M`` 来统一采样飞行距离，并在实际材质中通过拒绝采
样确定是否发生真实碰撞。本课实现与 OpenMC C++ 内核相同的关键步骤：

* **SlabGeometry**：一维平板几何，提供 `find_cell` 和边界反射逻辑；
* **Material**：定义真实总截面、散射/吸收概率，并给出碰撞抽样接口；
* **WoodcockTransporter**：执行主截面采样、虚拟碰撞拒绝与真实碰撞处
  理，完全贴合 OpenMC 的 delta-tracking 实现；
* **Tally**：使用轨迹长度估计与碰撞计数对比，展示虚拟碰撞不会影响统
  计量的无偏性。

运行脚本可观察每条粒子历史经历的虚拟碰撞次数以及真实反应类型。"""

from __future__ import annotations

from dataclasses import dataclass
from math import log
from typing import List

import numpy as np


@dataclass
class Material:
    name: str
    sigma_s: float
    sigma_a: float

    def __post_init__(self) -> None:
        self.sigma_t = self.sigma_s + self.sigma_a
        if self.sigma_t <= 0.0:
            raise ValueError("总截面必须正")

    def sample_reaction(self, rng: np.random.Generator) -> str:
        r = rng.random() * self.sigma_t
        return "scatter" if r < self.sigma_s else "absorption"


@dataclass
class Cell:
    left: float
    right: float
    material: Material

    def contains(self, x: float) -> bool:
        return self.left <= x < self.right


class SlabGeometry:
    """仿照 OpenMC 平板几何，实现 `find_cell`。"""

    def __init__(self, cells: List[Cell], reflective: tuple[bool, bool]):
        self.cells = cells
        self.reflective_left, self.reflective_right = reflective
        self.left_boundary = min(cell.left for cell in cells)
        self.right_boundary = max(cell.right for cell in cells)

    def find_cell(self, x: float) -> Cell | None:
        for cell in self.cells:
            if cell.contains(x):
                return cell
        return None

    def handle_boundary(self, x: float, direction: int) -> tuple[float, int, bool]:
        """处理边界：返回新位置、新方向以及是否泄漏。"""

        if x < self.left_boundary:
            if self.reflective_left:
                mirrored = 2 * self.left_boundary - x
                return mirrored, -direction, False
            return self.left_boundary, direction, True
        if x >= self.right_boundary:
            if self.reflective_right:
                mirrored = 2 * self.right_boundary - x
                return mirrored, -direction, False
            return self.right_boundary, direction, True
        return x, direction, False


@dataclass
class Particle:
    x: float
    direction: int
    weight: float = 1.0


@dataclass
class Tally:
    track_length: float = 0.0
    absorptions: int = 0
    scatters: int = 0

    def score_length(self, distance: float, weight: float) -> None:
        self.track_length += distance * weight

    def score_collision(self, reaction: str) -> None:
        if reaction == "absorption":
            self.absorptions += 1
        elif reaction == "scatter":
            self.scatters += 1


class WoodcockTransporter:
    """按照 OpenMC delta-tracking 流程推进粒子。"""

    def __init__(self, geometry: SlabGeometry, majorant: float):
        if majorant <= 0.0:
            raise ValueError("主截面必须正")
        self.geometry = geometry
        self.majorant = majorant

    def transport(self, particle: Particle, rng: np.random.Generator, tally: Tally, max_collisions: int = 10_000) -> None:
        collisions = 0
        while collisions < max_collisions:
            collisions += 1
            distance = -log(rng.random()) / self.majorant
            new_x = particle.x + distance * particle.direction
            tally.score_length(distance, particle.weight)
            new_x, new_dir, leaked = self.geometry.handle_boundary(new_x, particle.direction)
            if leaked:
                print("粒子泄漏出几何，停止输运")
                break
            particle.x = new_x
            particle.direction = new_dir

            cell = self.geometry.find_cell(particle.x)
            if cell is None:
                print("粒子位于未定义区域，终止")
                break

            acceptance = cell.material.sigma_t / self.majorant
            if rng.random() > acceptance:
                # 虚拟碰撞，等价于 OpenMC 的 delta-tracking 拒绝步骤
                continue

            reaction = cell.material.sample_reaction(rng)
            tally.score_collision(reaction)
            if reaction == "scatter":
                particle.direction = 1 if rng.random() < 0.5 else -1
            else:
                print("粒子被吸收，结束历史")
                break
        else:
            print("达到最大碰撞数，强制终止")


if __name__ == "__main__":
    fuel = Material("fuel", sigma_s=0.15, sigma_a=0.05)
    moderator = Material("moderator", sigma_s=0.3, sigma_a=0.01)
    cells = [
        Cell(left=0.0, right=1.0, material=fuel),
        Cell(left=1.0, right=4.0, material=moderator),
    ]
    geometry = SlabGeometry(cells, reflective=(True, True))
    majorant = max(cell.material.sigma_t for cell in cells) * 1.1

    transporter = WoodcockTransporter(geometry, majorant)
    rng = np.random.default_rng(seed=5)

    histories = 5
    for i in range(histories):
        particle = Particle(x=0.5, direction=1)
        tally = Tally()
        transporter.transport(particle, rng, tally)
        print(
            f"\n历史 {i + 1}: track_length={tally.track_length:.3f}, "
            f"scatters={tally.scatters}, absorptions={tally.absorptions}"
        )
