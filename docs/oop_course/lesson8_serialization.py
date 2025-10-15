"""Lesson 8: 对象协议与序列化

OpenMC 对象通常提供 ``to_xml_element``、``from_hdf5`` 等方法，将 Python 数据结构与文件
格式连接。本示例使用标准库 `xml.etree.ElementTree` 与 `json` 展示如何实现双向序列化，
同时覆盖 `__repr__`、`__contains__` 等友好协议。
"""

from __future__ import annotations

import json
import xml.etree.ElementTree as ET
from dataclasses import dataclass, asdict
from typing import Dict


@dataclass
class BoundingBox:
    """可序列化的几何包围盒。"""

    xmin: float
    xmax: float
    ymin: float
    ymax: float

    def to_xml_element(self) -> ET.Element:
        element = ET.Element("bounding_box")
        for key, value in asdict(self).items():
            element.set(key, str(value))
        return element

    @classmethod
    def from_xml_element(cls, element: ET.Element) -> "BoundingBox":
        kwargs = {key: float(value) for key, value in element.attrib.items()}
        return cls(**kwargs)

    def to_json(self) -> str:
        return json.dumps(asdict(self))

    @classmethod
    def from_json(cls, data: str) -> "BoundingBox":
        kwargs = json.loads(data)
        return cls(**kwargs)

    def __repr__(self) -> str:
        return (
            "BoundingBox(" + ", ".join(f"{k}={v}" for k, v in asdict(self).items()) + ")"
        )

    def __contains__(self, point: tuple[float, float]) -> bool:
        x, y = point
        return self.xmin <= x <= self.xmax and self.ymin <= y <= self.ymax


@dataclass
class Geometry:
    """包含包围盒的几何对象。"""

    name: str
    bounding_box: BoundingBox

    def to_xml(self) -> ET.Element:
        root = ET.Element("geometry", name=self.name)
        root.append(self.bounding_box.to_xml_element())
        return root

    def to_hdf5_payload(self) -> Dict[str, float]:
        """示例化 HDF5 序列化：OpenMC 实际使用 h5py，我们用 dict 表示。"""

        return asdict(self.bounding_box)

    def __contains__(self, point: tuple[float, float]) -> bool:
        return point in self.bounding_box


if __name__ == "__main__":
    box = BoundingBox(xmin=0.0, xmax=1.0, ymin=-1.0, ymax=1.0)
    geometry = Geometry(name="pin", bounding_box=box)

    xml_element = geometry.to_xml()
    xml_string = ET.tostring(xml_element, encoding="unicode")
    print("XML 表示:")
    print(xml_string)

    json_data = box.to_json()
    print("JSON 表示:")
    print(json_data)

    restored_box = BoundingBox.from_json(json_data)
    print(f"还原后的包围盒: {restored_box}")

    print(f"点 (0.5, 0.0) 是否在包围盒内? { (0.5, 0.0) in box }")
