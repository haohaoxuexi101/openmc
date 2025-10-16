"""Lesson 16: Unionized energy grid construction and lookup.

This script mirrors OpenMC's approach for accelerating Doppler-broadened
cross section lookups: every nuclide table contributes its resonance
structure to a common *unionized* grid so that all nuclides can be sampled
with the same energy index.  The code is intentionally compact but complete,
showing how an object-oriented design keeps responsibilities separated.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Sequence

import bisect
import numpy as np


@dataclass(frozen=True)
class CrossSectionTable:
    """Immutable cross section table for a single nuclide.

    We freeze the dataclass so instances are hashable and safe to cache.
    Each table stores monotonic energy breakpoints with *piecewise constant*
    values.  OpenMC performs more sophisticated interpolation but the lookup
    interface is identical, which is what the lesson focuses on.
    """

    nuclide: str
    energies: Sequence[float]
    total_xs: Sequence[float]

    def evaluate(self, energy: float) -> float:
        """Return the total cross section at ``energy`` using left sampling."""
        idx = bisect.bisect_right(self.energies, energy) - 1
        idx = max(0, min(idx, len(self.total_xs) - 1))
        return float(self.total_xs[idx])


class UnionizedEnergyGrid:
    """Aggregate many tables into one shared energy grid.

    The grid is responsible for two things:

    * computing the sorted union of all breakpoints;
    * precomputing table indices so lookups become O(1).

    Both steps mirror the C++ backend's ``EnergyGrid`` helpers.  Splitting
    the logic into a small class keeps the lesson short while highlighting
    why OpenMC factors this capability out of the material objects.
    """

    def __init__(self, tables: Sequence[CrossSectionTable]):
        self._tables = list(tables)
        self.grid = self._build_unionized_grid()
        self._index_map = self._precompute_indices()

    def _build_unionized_grid(self) -> np.ndarray:
        points: List[float] = []
        for table in self._tables:
            points.extend(table.energies)
        # ``np.unique`` is fast and deterministic; we keep it for clarity.
        grid = np.unique(points)
        if grid[0] > 0.0:
            # OpenMC always ensures the grid starts at zero to match
            # logarithmic interpolation.  We emulate that behaviour here.
            grid = np.insert(grid, 0, 0.0)
        return grid

    def _precompute_indices(self) -> Dict[str, np.ndarray]:
        mapping: Dict[str, np.ndarray] = {}
        for table in self._tables:
            idx = np.searchsorted(table.energies, self.grid, side="right") - 1
            idx = np.clip(idx, 0, len(table.total_xs) - 1)
            mapping[table.nuclide] = idx.astype(np.int32)
        return mapping

    def evaluate(self, energy: float) -> Dict[str, float]:
        """Return total cross sections for all tables at ``energy``.

        The method reflects OpenMC's ``Material.calculate_xs`` behaviour: the
        union grid index is located once, then reused for every nuclide.
        """

        union_idx = np.searchsorted(self.grid, energy, side="right") - 1
        union_idx = int(np.clip(union_idx, 0, len(self.grid) - 1))
        return {
            table.nuclide: float(table.total_xs[self._index_map[table.nuclide][union_idx]])
            for table in self._tables
        }


if __name__ == "__main__":
    u235 = CrossSectionTable(
        "U235",
        energies=[0.0253, 1.0, 10.0],
        total_xs=[680.0, 6.0, 2.0],
    )
    u238 = CrossSectionTable(
        "U238",
        energies=[0.0253, 0.5, 6.0],
        total_xs=[11.0, 4.0, 3.0],
    )
    grid = UnionizedEnergyGrid([u235, u238])

    energy = 0.3
    xs = grid.evaluate(energy)
    for nuclide, value in xs.items():
        print(f"{nuclide} total XS at {energy:.2f} eV = {value:.3f} barns")
