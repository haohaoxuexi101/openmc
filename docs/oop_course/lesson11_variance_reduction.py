"""Lesson 11: 方差缩减与 Russian Roulette 算法

本示例模拟 OpenMC C++ 求解器中的 `russian_roulette` 与 `split_particle`
逻辑，演示如何通过权重调节在蒙特卡洛输运中控制计算代价与统计误差。
实现要点：

* ``Particle`` 保存权重、历史和存活标记；
* ``russian_roulette`` 根据存活概率 `p_survive` 随机决定粒子是否存活，
  并对存活粒子的权重进行放大，保持无偏性；
* ``split_particle`` 将权重过大的粒子拆分成多个等权粒子；
* ``adaptive_control`` 结合两种策略，使粒子权重保持在预设窗口内，
  仿照 OpenMC 的 weight window 策略。

通过运行示例，可观察不同初始权重如何被归一到 ``[0.25, 4.0]`` 的权
重窗口中，同时记录发生的运算。
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Tuple

import numpy as np


@dataclass
class Particle:
    weight: float
    history: List[str] = field(default_factory=list)
    alive: bool = True

    def log(self, message: str) -> None:
        self.history.append(message)


def russian_roulette(particle: Particle, p_survive: float) -> None:
    """仿照 OpenMC 的 Russian Roulette 实现。"""

    if not particle.alive:
        return

    if particle.weight <= p_survive:
        particle.log(
            f"weight={particle.weight:.3f} <= p_survive -> skipping roulette"
        )
        return

    xi = np.random.random()
    if xi < p_survive:
        particle.weight /= p_survive
        particle.log(
            f"roulette survive (xi={xi:.3f}) -> new weight={particle.weight:.3f}"
        )
    else:
        particle.alive = False
        particle.log(f"roulette killed (xi={xi:.3f})")


def split_particle(particle: Particle, max_weight: float) -> List[Particle]:
    """将权重过大的粒子拆分成多个子粒子，仿照 OpenMC 的 ``split_particle``。"""

    if not particle.alive or particle.weight <= max_weight:
        return [particle]

    n = int(np.ceil(particle.weight / max_weight))
    new_weight = particle.weight / n
    particle.log(
        f"splitting weight={particle.weight:.3f} into {n} children (w={new_weight:.3f})"
    )
    particle.alive = False

    children = [Particle(new_weight, history=particle.history + [f"child #{i}"]) for i in range(n)]
    for child in children:
        child.log("spawned from split")
    return children


def adaptive_control(
    particles: List[Particle], window: Tuple[float, float], p_survive: float
) -> List[Particle]:
    """组合 Russian Roulette 与 splitting，将粒子权重保持在窗口内。"""

    w_min, w_max = window
    managed: List[Particle] = []

    for particle in particles:
        if particle.weight < w_min:
            particle.log(f"weight below {w_min}, applying roulette")
            russian_roulette(particle, p_survive)
            if particle.alive:
                managed.append(particle)
        elif particle.weight > w_max:
            particle.log(f"weight above {w_max}, splitting")
            managed.extend(split_particle(particle, w_max))
        else:
            particle.log("weight already in window")
            managed.append(particle)

    return managed


if __name__ == "__main__":
    np.random.seed(3)

    initial = [Particle(w) for w in (0.1, 0.5, 2.0, 6.5)]
    window = (0.25, 4.0)

    active = initial
    generation = 0
    while active:
        generation += 1
        print(f"\n--- generation {generation} ---")
        active = adaptive_control(active, window, p_survive=0.4)
        for idx, particle in enumerate(active):
            status = "alive" if particle.alive else "dead"
            print(f"particle {idx}: weight={particle.weight:.3f}, status={status}")
            for record in particle.history:
                print(f"  - {record}")

        # 终止条件：所有粒子均在窗口内且上一轮未进行分裂/淘汰
        if all(p.alive and window[0] <= p.weight <= window[1] for p in active):
            break
