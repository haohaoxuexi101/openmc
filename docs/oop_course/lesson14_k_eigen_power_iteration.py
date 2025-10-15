"""Lesson 14: k-特征值幂迭代与裂变源库重采样

本示例完整复刻 OpenMC 本征值计算的关键算法：

* **SourceBank**：和 C++ 内核中的裂变源库数据结构等价，存储粒子
  位置、方向、能量与权重；
* **power_iteration**：实现和 OpenMC `power_iteration` 同样的幂迭代流
  程，利用生成代之间的裂变中子比值估算 ``k_eff``；
* **systematic_resample**：使用系统重采样（OpenMC 的 "comb" 算法）将
  新的裂变源归一化到固定粒子数，保持权重平衡；
* **simulate_history**：在无限均匀介质中进行碰撞采样，包含散射、吸
  收和裂变分支，产额来自 ``nu\Sigma_f``。

运行脚本可以看到每一代的 ``k_eff`` 估计值和源库重采样结果，代码使
用的概率与分支和 OpenMC C++ 端保持一致。"""

from __future__ import annotations

from dataclasses import dataclass, replace
from math import log
from typing import Iterable, List

import numpy as np


@dataclass
class SourceSite:
    """与 OpenMC ``SourceSite`` 相同字段，用于存储裂变源信息。"""

    position: float  # 1D 无穷介质位置
    direction: int  # ±1 代表运动方向
    energy: float
    weight: float


@dataclass
class InfiniteMedium:
    """提供截面与概率，用于复现 OpenMC 的中子历史。"""

    sigma_s: float
    sigma_f: float
    sigma_a: float
    nu: float

    def __post_init__(self) -> None:
        self.sigma_t = self.sigma_s + self.sigma_f + self.sigma_a
        if self.sigma_t <= 0.0:
            raise ValueError("总截面必须大于 0")

    def sample_distance(self, rng: np.random.Generator) -> float:
        """使用 ``-ln(x)/Σ_t`` 采样碰撞距离。"""

        return -log(rng.random()) / self.sigma_t

    def sample_reaction(self, rng: np.random.Generator) -> str:
        r = rng.random() * self.sigma_t
        if r < self.sigma_s:
            return "scatter"
        if r < self.sigma_s + self.sigma_f:
            return "fission"
        return "absorption"

    def sample_fission_yield(self, rng: np.random.Generator) -> int:
        """按照泊松分布采样裂变产生的二次粒子数。"""

        mean = self.nu
        return rng.poisson(mean)

    def scatter(self, site: SourceSite, rng: np.random.Generator) -> SourceSite:
        """各向同性散射，方向随机翻转。"""

        direction = 1 if rng.random() < 0.5 else -1
        return replace(site, direction=direction)

    def fission_products(self, site: SourceSite, multiplicity: int, rng: np.random.Generator) -> list[SourceSite]:
        """生成裂变子代粒子，位置沿原方向推进一个平均自由程。"""

        children: list[SourceSite] = []
        for _ in range(multiplicity):
            distance = self.sample_distance(rng)
            new_pos = site.position + site.direction * distance
            direction = 1 if rng.random() < 0.5 else -1
            children.append(SourceSite(position=new_pos, direction=direction, energy=site.energy, weight=site.weight))
        return children


def systematic_resample(bank: Iterable[SourceSite], target_size: int, rng: np.random.Generator) -> list[SourceSite]:
    """与 OpenMC 裂变源重采样算法相同的系统重采样。"""

    bank_list = list(bank)
    weights = np.array([site.weight for site in bank_list], dtype=float)
    total_weight = weights.sum()
    if total_weight == 0.0:
        raise RuntimeError("裂变库权重为零，无法继续迭代")

    normalized = weights / total_weight
    cdf = np.cumsum(normalized)
    step = 1.0 / target_size
    start = rng.random() * step
    positions = start + step * np.arange(target_size)

    resampled: list[SourceSite] = []
    idx = 0
    for pos in positions:
        while pos > cdf[idx]:
            idx += 1
        new_site = replace(bank_list[idx], weight=total_weight / target_size)
        resampled.append(new_site)
    return resampled


def simulate_history(site: SourceSite, medium: InfiniteMedium, rng: np.random.Generator) -> list[SourceSite]:
    """推进单个粒子的历史，返回其产生的裂变二次粒子。"""

    children: list[SourceSite] = []
    current = site
    while True:
        distance = medium.sample_distance(rng)
        current = replace(current, position=current.position + current.direction * distance)
        reaction = medium.sample_reaction(rng)
        if reaction == "scatter":
            current = medium.scatter(current, rng)
            continue
        if reaction == "fission":
            multiplicity = medium.sample_fission_yield(rng)
            if multiplicity:
                children.extend(medium.fission_products(current, multiplicity, rng))
        break
    return children


def power_iteration(
    medium: InfiniteMedium, initial_bank: List[SourceSite], generations: int, rng: np.random.Generator
) -> list[float]:
    """执行与 OpenMC 相同的 k 有效幂迭代。"""

    k_eigenvalues: list[float] = []
    bank = initial_bank
    for gen in range(generations):
        produced: list[SourceSite] = []
        source_weight = sum(site.weight for site in bank)
        for site in bank:
            produced.extend(simulate_history(site, medium, rng))
        produced_weight = sum(site.weight for site in produced)
        if source_weight == 0.0:
            raise RuntimeError("源粒子权重为零，无法估计 k 有效")
        k_eff = produced_weight / source_weight
        k_eigenvalues.append(k_eff)
        bank = systematic_resample(produced, len(initial_bank), rng)
        print(f"第 {gen + 1:2d} 代: k_eff = {k_eff:.5f}, 裂变源粒子数 = {len(bank)}")
    return k_eigenvalues


if __name__ == "__main__":
    rng = np.random.default_rng(seed=12)
    medium = InfiniteMedium(sigma_s=0.2, sigma_f=0.08, sigma_a=0.02, nu=2.43)
    initial_bank = [SourceSite(position=0.0, direction=1, energy=2.0e6, weight=1.0) for _ in range(1000)]
    history = power_iteration(medium, initial_bank, generations=5, rng=rng)
    print("\n迭代完成，k_eff 序列:")
    for idx, value in enumerate(history, 1):
        print(f"  Gen {idx}: {value:.5f}")
