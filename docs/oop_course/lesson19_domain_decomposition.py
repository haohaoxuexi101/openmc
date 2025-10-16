"""Lesson 19: Spatial domain decomposition with recursive bisection.

OpenMC distributes particles across domains using orthogonal recursive
bisection (ORB).  This Python recreation mirrors that algorithm while
staying short: a ``BoundingBox`` value object, a decomposition orchestrator,
and a clear separation between geometry data and the partitioning strategy.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, List, Sequence, Tuple

import numpy as np


@dataclass(frozen=True)
class BoundingBox:
    """Immutable axis-aligned bounding box."""

    lower: np.ndarray
    upper: np.ndarray

    def extent(self) -> np.ndarray:
        return self.upper - self.lower

    def split(self, axis: int, value: float) -> Tuple["BoundingBox", "BoundingBox"]:
        lower_upper = self.upper.copy()
        lower_upper[axis] = value
        upper_lower = self.lower.copy()
        upper_lower[axis] = value
        return (
            BoundingBox(self.lower.copy(), lower_upper),
            BoundingBox(upper_lower, self.upper.copy()),
        )


@dataclass
class GeometryDomain:
    """Captures a geometry region and its weight (e.g. number of particles)."""

    name: str
    bounds: BoundingBox
    weight: float


class ORBPartitioner:
    """Recursive orthogonal bisection with equal-weight target."""

    def __init__(self, max_domains: int) -> None:
        self.max_domains = max_domains

    def partition(self, regions: Sequence[GeometryDomain]) -> List[List[GeometryDomain]]:
        regions = list(regions)
        result: List[List[GeometryDomain]] = []
        self._split(regions, self._total_weight(regions), result)
        return result

    def _split(
        self,
        regions: List[GeometryDomain],
        target_weight: float,
        result: List[List[GeometryDomain]],
    ) -> None:
        if len(result) >= self.max_domains - 1 or len(regions) == 1:
            result.append(regions)
            return

        axis = int(np.argmax(regions[0].bounds.extent()))
        pivot = self._find_partition_value(regions, axis, target_weight / 2.0)
        left: List[GeometryDomain] = []
        right: List[GeometryDomain] = []
        for region in regions:
            box_left, box_right = region.bounds.split(axis, pivot)
            split_weight = region.weight / 2.0
            left.append(GeometryDomain(region.name + "_L", box_left, split_weight))
            right.append(GeometryDomain(region.name + "_R", box_right, split_weight))

        self._split(left, target_weight / 2.0, result)
        self._split(right, target_weight / 2.0, result)

    @staticmethod
    def _total_weight(regions: Iterable[GeometryDomain]) -> float:
        return sum(region.weight for region in regions)

    @staticmethod
    def _find_partition_value(
        regions: Sequence[GeometryDomain], axis: int, target: float
    ) -> float:
        values = sorted({float(region.bounds.lower[axis]) for region in regions} |
                        {float(region.bounds.upper[axis]) for region in regions})
        cumulative = 0.0
        for value in values:
            cumulative += sum(
                region.weight
                for region in regions
                if region.bounds.lower[axis] <= value < region.bounds.upper[axis]
            )
            if cumulative >= target:
                return value
        return values[-1]


if __name__ == "__main__":
    regions = [
        GeometryDomain(
            name="fuel",
            bounds=BoundingBox(np.array([0.0, 0.0, 0.0]), np.array([1.0, 1.0, 1.0])),
            weight=0.6,
        ),
        GeometryDomain(
            name="moderator",
            bounds=BoundingBox(np.array([1.0, 0.0, 0.0]), np.array([2.0, 1.0, 1.0])),
            weight=0.4,
        ),
    ]
    partitioner = ORBPartitioner(max_domains=4)
    domains = partitioner.partition(regions)
    for i, domain in enumerate(domains, start=1):
        print(f"Domain {i} total weight = {sum(r.weight for r in domain):.2f}")
