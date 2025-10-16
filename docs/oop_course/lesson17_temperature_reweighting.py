"""Lesson 17: Temperature reweighting and cache-aware material design.

OpenMC stores temperature-dependent cross sections in an on-demand cache.
Here we demonstrate the same idea with compact classes: a base interface for
reconstruction strategies, a cache that respects material state, and a tiny
example that interpolates between tabulated temperatures.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Protocol, Tuple

import numpy as np


class ReconstructionStrategy(Protocol):
    """Strategy protocol for rebuilding cross sections at a target temperature."""

    def reconstruct(self, base: np.ndarray, temperature: float) -> np.ndarray:
        ...


class LinearTemperatureInterpolator:
    """Mimics OpenMC's logarithmic temperature interpolation in spirit.

    The implementation remains linear for brevity but preserves the
    plug-in strategy pattern: materials can swap reconstruction algorithms
    without changing their public interface.
    """

    def __init__(self, tabulated: Dict[float, np.ndarray]):
        self._tabulated = dict(sorted(tabulated.items()))

    def reconstruct(self, base: np.ndarray, temperature: float) -> np.ndarray:
        lower, upper = self._bracketing_temperatures(temperature)
        if lower == upper:
            return self._tabulated[lower]
        frac = (temperature - lower) / (upper - lower)
        return (1.0 - frac) * self._tabulated[lower] + frac * self._tabulated[upper]

    def _bracketing_temperatures(self, temperature: float) -> Tuple[float, float]:
        temps = list(self._tabulated)
        if temperature <= temps[0] or temperature >= temps[-1]:
            # Clamp to available data, mirroring OpenMC's error handling.
            return temps[0], temps[-1]
        for lower, upper in zip(temps, temps[1:]):
            if lower <= temperature <= upper:
                return lower, upper
        raise RuntimeError("Failed to bracket temperature")


@dataclass
class MaterialState:
    """Stores the cross section cache keyed by temperature."""

    name: str
    strategy: ReconstructionStrategy
    base_xs: np.ndarray
    _cache: Dict[float, np.ndarray] = field(default_factory=dict)

    def xs_at(self, temperature: float) -> np.ndarray:
        """Return cross sections at ``temperature`` using cache-on-access."""
        key = round(float(temperature), 2)  # mimic OpenMC's tolerance rounding
        if key not in self._cache:
            self._cache[key] = self.strategy.reconstruct(self.base_xs, key)
        return self._cache[key]


if __name__ == "__main__":
    base = np.array([100.0, 50.0, 20.0])
    tabulated = {
        300.0: base,
        600.0: base * 0.9,
        900.0: base * 0.8,
    }
    fuel = MaterialState(
        name="Fuel",
        strategy=LinearTemperatureInterpolator(tabulated),
        base_xs=base,
    )
    print("XS @ 450 K:", fuel.xs_at(450.0))
    print("XS @ 450 K (cached):", fuel.xs_at(450.0))
    print("XS @ 1000 K (clamped):", fuel.xs_at(1000.0))
