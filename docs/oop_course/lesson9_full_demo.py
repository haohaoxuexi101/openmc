"""Lesson 9: 综合示例

本章节将前面所有部件组合成一个迷你版的 OpenMC 建模流程。你将看到如何使用 ID 管理、
抽象基类、属性封装、类型安全容器以及序列化，将一份几何 + 材料 + 计数器的配置导出为
XML/JSON 描述，同时执行简单的算术操作。
"""

from __future__ import annotations

import xml.etree.ElementTree as ET

from lesson1_id_manager import Material, Surface
from lesson2_equality_mixin import Spectrum
from lesson3_geometry_base import Cell, LatticeUniverse, SimpleUniverse
from lesson4_filter_meta import CellFilter, EnergyFilter
from lesson5_material_properties import Material as CheckedMaterial
from lesson6_checked_list import CheckedList, Score
from lesson7_tally_composition import Tally, make_checked_strings
from lesson8_serialization import BoundingBox, Geometry

import numpy as np


def build_geometry() -> SimpleUniverse:
    fuel = Cell(name="fuel", fill="UO2")
    moderator = Cell(name="moderator", fill="Water")

    root = SimpleUniverse(name="core")
    root.add_cell(fuel)
    root.add_cell(moderator)

    lattice = LatticeUniverse(name="assembly", width=1, height=2)
    lattice.add_cell(fuel, position=(0, 0))
    lattice.add_cell(moderator, position=(0, 1))
    print(lattice.describe(indent=2))

    return root


def create_materials() -> list[CheckedMaterial]:
    fuel = CheckedMaterial(name="UO2", density=10.4, temperature=1200.0)
    fuel.volume = 1.5
    fuel.depletable = True

    water = CheckedMaterial(name="Water", density=1.0, temperature=600.0)
    water.volume = 3.0
    return [fuel, water]


def setup_tallies() -> list[Tally]:
    filters = make_checked_strings(["cell 1"])
    energy_filter = EnergyFilter(bins=[0.0, 1e3, 20e6])
    cell_filter = CellFilter(cell_ids=[1])

    scores = CheckedList(Score, [Score("flux"), Score("fission")])
    nuclides = make_checked_strings(["U235"])

    tally_main = Tally(
        name="main",
        filters=filters,
        scores=scores,
        nuclides=nuclides,
        values=[1.2, 3.4],
    )
    tally_aux = tally_main.scale(0.5)
    combined = tally_main.combine(tally_aux)

    print("计数器描述：")
    print(combined.describe())
    print(f"能量滤器 short_name = {energy_filter.short_name}")
    print(f"Cell 滤器 short_name = {cell_filter.short_name}")

    return [tally_main, combined]


def export_geometry() -> None:
    box = BoundingBox(xmin=0.0, xmax=1.0, ymin=-1.0, ymax=1.0)
    geometry = Geometry(name="pin", bounding_box=box)
    xml_string = ET.tostring(geometry.to_xml(), encoding="unicode")
    print("几何 XML:")
    print(xml_string)
    print("包围盒 JSON:")
    print(box.to_json())


def check_spectrum_equality() -> None:
    energy_groups = np.array([0.0, 1.0, 20.0])
    ref = Spectrum(groups=energy_groups, values=np.array([0.3, 0.7]))
    trial = Spectrum(groups=energy_groups.copy(), values=np.array([0.3, 0.7]))
    print(f"能谱对象是否相等? {ref == trial}")


def demonstrate_id_management() -> None:
    Surface.reset_ids()
    Material.reset_ids()
    s1 = Surface(name="vacuum")
    s2 = Surface(name="fuel-clad")
    uo2 = Material(name="UO2", density=10.4)
    print(f"Surface ID: {s1.id}, {s2.id}; Material ID: {uo2.id}")


if __name__ == "__main__":
    build_geometry()
    materials = create_materials()
    for mat in materials:
        print(mat.describe())

    setup_tallies()
    export_geometry()
    check_spectrum_equality()
    demonstrate_id_management()
