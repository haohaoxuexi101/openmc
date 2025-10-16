"""Lesson 18: Event-based tallies with observer pattern.

The C++ kernel in OpenMC emits collision and surface events which the tally
module observes.  To keep the lesson concise we model a stripped-down event
system where tallies subscribe to a dispatcher.  The design demonstrates how
mixins and small classes cooperate to deliver a feature-rich API.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Protocol


class Event(Protocol):
    weight: float


@dataclass
class CollisionEvent:
    cell_id: int
    energy: float
    weight: float = 1.0


@dataclass
class SurfaceEvent:
    surface_id: int
    current: float
    weight: float = 1.0


class EventDispatcher:
    """Mixin that handles observer registration and broadcasting."""

    def __init__(self) -> None:
        self._observers: List[EventObserver] = []

    def register(self, observer: "EventObserver") -> None:
        self._observers.append(observer)

    def notify(self, event: Event) -> None:
        for observer in self._observers:
            observer.on_event(event)


class EventObserver(Protocol):
    def on_event(self, event: Event) -> None:
        ...


@dataclass
class CellCollisionTally:
    """Accumulates path-length tallies per cell."""

    cells: Iterable[int]
    scores: Dict[int, float] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.scores = {cell: 0.0 for cell in self.cells}

    def on_event(self, event: Event) -> None:
        if isinstance(event, CollisionEvent) and event.cell_id in self.scores:
            self.scores[event.cell_id] += event.weight


@dataclass
class SurfaceCurrentTally:
    """Records net current through surfaces, matching OpenMC's behaviour."""

    surfaces: Iterable[int]
    scores: Dict[int, float] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.scores = {surf: 0.0 for surf in self.surfaces}

    def on_event(self, event: Event) -> None:
        if isinstance(event, SurfaceEvent) and event.surface_id in self.scores:
            self.scores[event.surface_id] += event.current * event.weight


if __name__ == "__main__":
    dispatcher = EventDispatcher()
    collision_tally = CellCollisionTally(cells=[1, 2])
    surface_tally = SurfaceCurrentTally(surfaces=[10])
    dispatcher.register(collision_tally)
    dispatcher.register(surface_tally)

    dispatcher.notify(CollisionEvent(cell_id=1, energy=0.5))
    dispatcher.notify(CollisionEvent(cell_id=2, energy=1.0, weight=0.2))
    dispatcher.notify(SurfaceEvent(surface_id=10, current=1.5))

    print("Collision tally:", collision_tally.scores)
    print("Surface tally:", surface_tally.scores)
