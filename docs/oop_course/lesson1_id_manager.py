"""Lesson 1: 唯一 ID 管理与 Mixin 复用

本示例以 OpenMC `IDManagerMixin` 的逻辑为模板，展示如何在多个几何元素类之间共享
统一的 ID 分配策略。通过面向对象的 Mixin，我们可以将 ID 生命周期管理与业务逻辑解耦，
并在派生类中获得自动编号、显式保留和重置等功能。
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import ClassVar, Dict, Optional, Set


class IDManagerMixin:
    """模仿 OpenMC 的 ID 管理混入类。

    子类通过继承该 Mixin 获得统一的 ID 分配/回收逻辑。我们保留了 OpenMC 的核心接口：

    * ``next_id``: 下一个可用 ID 的类属性
    * ``used_ids``: 已分配 ID 的集合
    * ``reset_ids``: 重置所有状态
    * ``reserve_id``: 保留外部指定的 ID，避免自动分配时冲突
    * ``get_next_id``: 自动生成 ID
    """

    next_id: ClassVar[int] = 1
    used_ids: ClassVar[Set[int]] = set()

    @classmethod
    def reset_ids(cls) -> None:
        """重置计数器，便于测试与重新建模。"""

        print(f"[IDManagerMixin] Resetting IDs for {cls.__name__}")
        cls.next_id = 1
        cls.used_ids.clear()

    @classmethod
    def reserve_id(cls, uid: int) -> None:
        """显式保留一个 ID，模仿从外部读取模型时的行为。"""

        if uid in cls.used_ids:
            raise ValueError(f"ID {uid} 已被 {cls.__name__} 使用")
        cls.used_ids.add(uid)
        cls.next_id = max(cls.next_id, uid + 1)
        print(f"[IDManagerMixin] Reserved ID {uid} for {cls.__name__}")

    @classmethod
    def get_next_id(cls) -> int:
        """返回下一个可用 ID，并推进计数器。"""

        while cls.next_id in cls.used_ids:
            cls.next_id += 1
        uid = cls.next_id
        cls.used_ids.add(uid)
        cls.next_id += 1
        print(f"[IDManagerMixin] Assigned ID {uid} for {cls.__name__}")
        return uid


@dataclass
class Surface(IDManagerMixin):
    """几何面，仿照 OpenMC 的 `Surface` 对象。"""

    name: str
    id: Optional[int] = field(default=None)

    def __post_init__(self) -> None:
        if self.id is None:
            self.id = self.get_next_id()
        else:
            self.reserve_id(self.id)


@dataclass
class Material(IDManagerMixin):
    """材料对象，展示 Mixin 的复用性。"""

    name: str
    density: float
    id: Optional[int] = field(default=None)

    def __post_init__(self) -> None:
        if self.id is None:
            self.id = self.get_next_id()
        else:
            self.reserve_id(self.id)


class ModelRegistry:
    """辅助类：跟踪所有构建对象，便于展示 ID 管理效果。"""

    def __init__(self) -> None:
        self.surfaces: Dict[int, Surface] = {}
        self.materials: Dict[int, Material] = {}

    def add_surface(self, surface: Surface) -> None:
        self.surfaces[surface.id] = surface

    def add_material(self, material: Material) -> None:
        self.materials[material.id] = material


if __name__ == "__main__":
    # 1) 自动分配 ID
    ModelRegistry()
    s1 = Surface(name="vacuum")
    s2 = Surface(name="fuel-clad")

    # 2) 手动保留 ID
    uo2 = Material(name="UO2", density=10.4, id=100)
    water = Material(name="Water", density=1.0)

    print(f"分配的 Surface ID: {[s1.id, s2.id]}")
    print(f"材料 ID: {[uo2.id, water.id]}")

    # 3) 重置
    Material.reset_ids()
    Surface.reset_ids()
