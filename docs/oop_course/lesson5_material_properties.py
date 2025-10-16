"""Lesson 5: 属性封装与输入校验

OpenMC 的 `Material`、`Cell` 等类通过 `@property` 将属性管理与校验集中化。本示例仿照
这种写法，实现材料温度、体积等字段的验证逻辑，确保对象始终处于合法状态。
"""

from __future__ import annotations

from typing import Optional


def _check_positive(value: float, name: str) -> float:
    if value <= 0:
        raise ValueError(f"{name} 必须为正数")
    return value


def _check_optional_positive(value: Optional[float], name: str) -> Optional[float]:
    if value is None:
        return None
    return _check_positive(value, name)


class Material:
    """带有属性校验的材料对象。"""

    def __init__(self, name: str, density: float, temperature: float = 900.0):
        self.name = name
        self.density = density
        self.temperature = temperature
        self._volume: Optional[float] = None
        self._depletable = False

    @property
    def density(self) -> float:
        return self._density

    @density.setter
    def density(self, value: float) -> None:
        self._density = _check_positive(value, "density")

    @property
    def temperature(self) -> float:
        return self._temperature

    @temperature.setter
    def temperature(self, value: float) -> None:
        self._temperature = _check_positive(value, "temperature")

    @property
    def volume(self) -> Optional[float]:
        return self._volume

    @volume.setter
    def volume(self, value: Optional[float]) -> None:
        self._volume = _check_optional_positive(value, "volume")

    @property
    def depletable(self) -> bool:
        return self._depletable

    @depletable.setter
    def depletable(self, value: bool) -> None:
        if not isinstance(value, bool):
            raise TypeError("depletable 必须是布尔值")
        self._depletable = value

    def describe(self) -> str:
        return (
            f"Material(name={self.name}, density={self.density}, temperature={self.temperature}, "
            f"volume={self.volume}, depletable={self.depletable})"
        )


if __name__ == "__main__":
    fuel = Material(name="UO2", density=10.4, temperature=1200.0)
    fuel.volume = 1.5
    fuel.depletable = True
    print(fuel.describe())

    try:
        fuel.temperature = -50
    except ValueError as exc:
        print(f"温度设置失败: {exc}")

    try:
        fuel.depletable = "yes"
    except TypeError as exc:
        print(f"耗尽设置失败: {exc}")
