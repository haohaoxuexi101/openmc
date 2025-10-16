"""Lesson 13: 燃耗求解器——CRAM + 预测-校正算法

OpenMC 的燃耗模块使用常规燃耗微分方程 ``dN/dt = A(N, \phi) N``，其中
``A`` 是由吸收、裂变、衰变产生的转化矩阵。本课程复刻 OpenMC 的核心流程：

* 构建 ``DepletionChain``，记录吸收、衰变与裂变产物关系；
* 根据瞬时粒子数密度 ``N`` 与常量中子通量 ``\phi`` 计算反应速率 ``R``；
* 使用和 OpenMC 相同的 16 阶复有理近似（CRAM-16）求解矩阵指数；
* 通过预测-校正积分器迭代更新反应速率，模拟 OpenMC `operator.py`
  中的 predictor-corrector 步骤。

运行脚本即可得到 1 步燃耗后的浓度变化，并输出每一步的矩阵与收敛情况。
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Mapping

import numpy as np


# --- CRAM-16 系数（来自 OpenMC ``openmc/deplete/cram.py``） ---
ALPHA0 = 2.1248537104952237488e-16
ALPHA = np.array([
    2.2580381827439833941e-01 + 2.2580381827439833941e-01j,
    2.2580381827439833941e-01 - 2.2580381827439833941e-01j,
    1.5870243430567845411e-01 + 2.4133130934759948988e-01j,
    1.5870243430567845411e-01 - 2.4133130934759948988e-01j,
    8.4149852912854629089e-02 + 2.4640213328375597619e-01j,
    8.4149852912854629089e-02 - 2.4640213328375597619e-01j,
    1.8584588010905935585e-02 + 2.4690947684390895430e-01j,
    1.8584588010905935585e-02 - 2.4690947684390895430e-01j,
])
THETA = np.array([
    -1.0843917078696988026e+00 + 8.9015877231026243790e+00j,
    -1.0843917078696988026e+00 - 8.9015877231026243790e+00j,
    -3.5091036084149180974e+00 + 3.7132745994695882233e+00j,
    -3.5091036084149180974e+00 - 3.7132745994695882233e+00j,
    -6.4161776990994341923e+00 + 1.5489960631534797325e+00j,
    -6.4161776990994341923e+00 - 1.5489960631534797325e+00j,
    -9.0149339281534924199e+00 + 1.0481478389620078816e-01j,
    -9.0149339281534924199e+00 - 1.0481478389620078816e-01j,
])


@dataclass
class ChainNuclide:
    """表示一个链上核素的转变关系。"""

    name: str
    decay_const: float = 0.0
    decay_paths: Dict[str, float] = field(default_factory=dict)
    capture: Dict[str, float] = field(default_factory=dict)
    fission_yields: Dict[str, float] = field(default_factory=dict)


@dataclass
class DepletionChain:
    """根据反应速率构造转化矩阵。"""

    nuclides: Mapping[str, ChainNuclide]

    def matrix(self, rates: Mapping[tuple[str, str], float]) -> np.ndarray:
        nucs = list(self.nuclides)
        index = {name: idx for idx, name in enumerate(nucs)}
        size = len(nucs)
        mat = np.zeros((size, size))
        for j, name in enumerate(nucs):
            nuclide = self.nuclides[name]
            removal = 0.0

            # 吸收（不含裂变）
            absorb_rate = rates.get((name, "absorption"), 0.0)
            removal += absorb_rate
            for product, branch in nuclide.capture.items():
                i = index[product]
                mat[i, j] += branch * absorb_rate

            # 裂变移除 + 裂变产额
            fiss_rate = rates.get((name, "fission"), 0.0)
            removal += fiss_rate
            for product, yield_ in nuclide.fission_yields.items():
                i = index[product]
                mat[i, j] += yield_ * fiss_rate

            # 衰变
            if nuclide.decay_const > 0.0:
                removal += nuclide.decay_const
                for product, branch in nuclide.decay_paths.items():
                    i = index[product]
                    mat[i, j] += branch * nuclide.decay_const

            mat[j, j] -= removal
        return mat


@dataclass
class ConstantFluxField:
    """根据常量通量和微观截面计算反应速率。"""

    flux: float
    microscopic_xs: Mapping[str, Dict[str, float]]

    def rates(self, densities: Mapping[str, float]) -> Dict[tuple[str, str], float]:
        rates: Dict[tuple[str, str], float] = {}
        for name, reactions in self.microscopic_xs.items():
            n = densities.get(name, 0.0)
            for rx, xs in reactions.items():
                rates[(name, rx)] = self.flux * xs * n
        return rates


class PredictorCorrectorIntegrator:
    """复刻 OpenMC 中的 predictor-corrector 燃耗积分器。"""

    def __init__(self, chain: DepletionChain, field: ConstantFluxField, dt: float):
        self.chain = chain
        self.field = field
        self.dt = dt

    def step(self, number_densities: Mapping[str, float]) -> Dict[str, float]:
        names = list(self.chain.nuclides)
        n0 = np.array([number_densities[name] for name in names], dtype=float)

        rates0 = self.field.rates(number_densities)
        A0 = self.chain.matrix(rates0)
        predictor = self._cram_step(A0, n0)

        pred_density = {name: predictor[i] for i, name in enumerate(names)}
        rates1 = self.field.rates(pred_density)
        A1 = self.chain.matrix(rates1)
        corrected = self._cram_step(0.5 * (A0 + A1), n0)

        return {name: corrected[i] for i, name in enumerate(names)}

    def _cram_step(self, matrix: np.ndarray, n0: np.ndarray) -> np.ndarray:
        identity = np.eye(matrix.shape[0])
        mat = matrix * self.dt
        result = ALPHA0 * n0
        for alpha, theta in zip(ALPHA, THETA):
            solve = np.linalg.solve(mat - theta * identity, n0)
            result += 2.0 * np.real(alpha * solve)
        return result


if __name__ == "__main__":
    chain = DepletionChain(
        {
            "U235": ChainNuclide(
                name="U235",
                capture={"U236": 1.0},
                fission_yields={"Xe135": 0.06},
            ),
            "U236": ChainNuclide(name="U236", capture={}),
            "Xe135": ChainNuclide(
                name="Xe135",
                capture={"Xe136": 1.0},
                decay_const=np.log(2) / (9.14 * 3600.0),  # 半衰期 9.14h
                decay_paths={"Cs135": 1.0},
            ),
            "Xe136": ChainNuclide(name="Xe136"),
            "Cs135": ChainNuclide(name="Cs135"),
        }
    )

    field = ConstantFluxField(
        flux=3.0e14,
        microscopic_xs={
            "U235": {"absorption": 680e-24, "fission": 585e-24},
            "U236": {"absorption": 12e-24},
            "Xe135": {"absorption": 2.6e6 * 1e-24},
        },
    )

    initial = {"U235": 6.5e-2, "U236": 0.0, "Xe135": 0.0, "Xe136": 0.0, "Cs135": 0.0}
    integrator = PredictorCorrectorIntegrator(chain, field, dt=24 * 3600.0)
    updated = integrator.step(initial)

    print("初始数密度:")
    for name, value in initial.items():
        print(f"  {name:6s} = {value:.6e}")

    print("\n燃耗 1 天后的数密度:")
    for name, value in updated.items():
        print(f"  {name:6s} = {value:.6e}")
