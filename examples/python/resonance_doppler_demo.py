r"""用可执行代码演示多极共振截面与多普勒展宽公式.

该脚本逐步实现 OpenMC 理论手册中的窗口多极（WMP）公式，并在终端打印详细解释：

* 0 K 截面：

  .. math::
     \sigma(E, 0) = \frac{1}{E} \sum_j \operatorname{Re}\left[\frac{i r_j}{\sqrt{E} - p_j}\right]

* 含修正积分 :math:`C` 的温度展宽截面：

  .. math::
     \sigma(E, T) = \frac{1}{2 E \sqrt{\xi}} \sum_j \operatorname{Re}\left[r_j\sqrt{\pi} W_i(z)
     - \frac{r_j}{\sqrt{\pi}} C\left(\frac{p_j}{\sqrt{\xi}}, \frac{u}{2\sqrt{\xi}}\right)\right]

* 忽略 :math:`C` 项时的高能近似：

  .. math::
     \sigma(E, T) \approx \frac{1}{2 E \sqrt{\xi}} \sum_j \operatorname{Re}\left[i r_j \sqrt{\pi} W_i(z)\right]

脚本最后会：

1. 打印基础物理量（温度、多普勒宽度等）；
2. 列出 0.01–10 eV 能量网格上的 0 K、温度展宽、近似截面及相对误差；
3. 对每个能量点，逐个极点给出贡献，帮助理解多普勒展宽如何平滑共振。

运行方法：``python examples/python/resonance_doppler_demo.py``。
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Tuple

import numpy as np
from scipy.integrate import quad
from scipy.special import wofz

KB_EV_PER_K = 8.617333262145e-5
"""玻尔兹曼常数（单位 eV/K），与 OpenMC 内部一致。"""


@dataclass
class MultipoleDataset:
    """封装单个能量窗口的多极参数。"""

    poles: np.ndarray
    residues: np.ndarray

    def __post_init__(self) -> None:
        self.poles = np.asarray(self.poles, dtype=np.complex128)
        self.residues = np.asarray(self.residues, dtype=np.complex128)
        if self.poles.shape != self.residues.shape:
            raise ValueError("poles and residues must have the same shape")


def zero_kelvin_cross_section(
    energies: Iterable[float],
    dataset: MultipoleDataset,
    return_components: bool = False,
) -> np.ndarray:
    """计算 0 K 多极截面，同时可选返回逐极点贡献。"""

    energies = np.asarray(energies, dtype=float)
    if np.any(energies <= 0.0):
        raise ValueError("energies must be strictly positive")

    sqrt_e = np.sqrt(energies)[:, None]
    term = 1j * dataset.residues / (sqrt_e - dataset.poles)
    contributions = term.real / energies[:, None]
    sigma = contributions.sum(axis=1)
    if return_components:
        return sigma, contributions
    return sigma


def _faddeeva_half_plane(z: np.ndarray) -> np.ndarray:
    """计算 :math:`W_i(z)` 并保证选取正确的半平面支路。"""

    wi = np.empty_like(z, dtype=np.complex128)
    mask = np.imag(z) >= 0.0
    wi[mask] = wofz(z[mask])
    conj_z = np.conjugate(z[~mask])
    wi[~mask] = -np.conjugate(wofz(conj_z))
    return wi


def _c_integral(p: complex, u: float, xi: float) -> complex:
    """计算单个极点的修正积分 :math:`C`。"""

    def integrand(u_prime: float) -> complex:
        exponent = np.exp(-((u + u_prime) ** 2) / (4.0 * xi))
        return exponent / (p**2 - u_prime**2)

    real = quad(
        lambda val: np.real(integrand(val)),
        0.0,
        np.inf,
        limit=250,
        epsabs=1e-10,
    )[0]
    imag = quad(
        lambda val: np.imag(integrand(val)),
        0.0,
        np.inf,
        limit=250,
        epsabs=1e-10,
    )[0]
    return 2.0 * p * (real + 1j * imag)


def doppler_broadened_cross_section(
    energies: Iterable[float],
    temperature: float,
    target_mass: float,
    dataset: MultipoleDataset,
    include_correction: bool = True,
    return_components: bool = False,
) -> np.ndarray:
    """计算温度展宽后的多极截面。"""

    energies = np.asarray(energies, dtype=float)
    if np.any(energies <= 0.0):
        raise ValueError("energies must be strictly positive")
    if temperature <= 0.0:
        raise ValueError("temperature must be strictly positive")

    xi = KB_EV_PER_K * temperature / (4.0 * target_mass)
    sqrt_e = np.sqrt(energies)
    z = (sqrt_e[:, None] - dataset.poles) / (2.0 * np.sqrt(xi))
    wi = _faddeeva_half_plane(z)

    correction = None
    if include_correction:
        correction = np.empty_like(z, dtype=np.complex128)
        for idx_e, u in enumerate(sqrt_e):
            for idx_p, pole in enumerate(dataset.poles):
                correction[idx_e, idx_p] = _c_integral(pole, u, xi)

    contributions = dataset.residues * (np.sqrt(np.pi) * wi)
    if include_correction:
        contributions -= dataset.residues * correction / np.sqrt(np.pi)

    prefactor = 1.0 / (2.0 * np.sqrt(xi) * energies)
    sigma = prefactor * contributions.real.sum(axis=1)
    if return_components:
        return sigma, (prefactor[:, None] * contributions.real)
    return sigma


def debug_pole_contributions(
    energies: np.ndarray,
    temperature: float,
    target_mass: float,
    dataset: MultipoleDataset,
) -> None:
    """打印每个能量点上各极点的贡献。"""

    sigma_0k, comp_0k = zero_kelvin_cross_section(
        energies, dataset, return_components=True
    )
    sigma_t, comp_t = doppler_broadened_cross_section(
        energies,
        temperature,
        target_mass,
        dataset,
        include_correction=True,
        return_components=True,
    )

    for idx, energy in enumerate(energies):
        print("--------------------------- 极点逐项贡献 ---------------------------")
        print(f"能量 E = {energy:11.5e} eV")
        for pole_idx in range(dataset.poles.size):
            frac = comp_t[idx, pole_idx] / sigma_t[idx] if sigma_t[idx] != 0.0 else 0.0
            print(
                "  "
                f"极点 #{pole_idx}: "
                f"σ_0K = {comp_0k[idx, pole_idx]:11.4e} b, "
                f"σ_T = {comp_t[idx, pole_idx]:11.4e} b, "
                f"贡献比例 {frac:5.2f}"
            )
        print(
            "  合计: σ_0K = {0:11.4e} b, σ_T = {1:11.4e} b".format(
                sigma_0k[idx], sigma_t[idx]
            )
        )


def demo_dataset() -> MultipoleDataset:
    """提供演示用的多极参数。实际应用中请替换为真实核素数据。"""

    poles = np.array(
        [
            0.28 + 0.23j,
            0.41 + 0.12j,
            0.64 + 0.29j,
            0.92 + 0.41j,
        ]
    )
    residues = np.array(
        [
            0.19 - 0.07j,
            0.08 - 0.11j,
            0.21 - 0.05j,
            0.05 - 0.02j,
        ]
    )
    return MultipoleDataset(poles=poles, residues=residues)


def tabulate_cross_sections(
    temperature: float = 900.0,
    target_mass: float = 238.0,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """生成对比表格所需的数据。"""

    dataset = demo_dataset()
    energies = np.logspace(-2, 1.0, num=8)  # 0.01 eV 到 10 eV
    sigma_0k = zero_kelvin_cross_section(energies, dataset)
    sigma_t = doppler_broadened_cross_section(
        energies,
        temperature,
        target_mass,
        dataset,
        include_correction=True,
    )
    sigma_t_approx = doppler_broadened_cross_section(
        energies,
        temperature,
        target_mass,
        dataset,
        include_correction=False,
    )
    return energies, sigma_0k, sigma_t, sigma_t_approx


def main() -> None:
    """主函数：打印表格并执行调试输出。"""

    temperature = 900.0  # K
    target_mass = 238.0  # 近似为 U-238
    dataset = demo_dataset()

    energies, sigma_0k, sigma_t, sigma_t_approx = tabulate_cross_sections(
        temperature=temperature,
        target_mass=target_mass,
    )

    xi = KB_EV_PER_K * temperature / (4.0 * target_mass)
    print("============================== 基本参数 ==============================")
    print(
        f"温度 T = {temperature:7.3f} K, 靶核质量数 A = {target_mass:7.3f}, 多普勒宽度 ξ = {xi:8.3e} eV"
    )
    print(
        f"能量网格: {energies.size:d} 个点, 范围 [{energies.min():.3e}, {energies.max():.3e}] eV"
    )
    print()

    header = (
        "------------------------------ 截面表 -------------------------------"
    )
    print(header)
    print(
        "    E [eV]   |   σ_0K [b]   |    σ_T [b]    | σ_T (忽略 C) [b] | 相对误差"
    )
    print("-" * len(header))
    for energy, s0, sT, sTapprox in zip(energies, sigma_0k, sigma_t, sigma_t_approx):
        rel_err = 0.0
        if sT != 0.0:
            rel_err = abs(sTapprox - sT) / abs(sT)
        print(
            f" {energy:11.5e} | {s0:12.6e} | {sT:12.6e} | {sTapprox:16.6e} | {rel_err:9.2e}"
        )
    print()

    debug_energies = np.array([energies[0], energies[len(energies) // 2], energies[-1]])
    debug_pole_contributions(
        debug_energies,
        temperature=temperature,
        target_mass=target_mass,
        dataset=dataset,
    )


if __name__ == "__main__":
    main()
