#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
surface_joint_stats.py

对 HDF5 面源文件（如 surface_source.h5）的 '/source_bank' 复合数据集：
  1. 提取方向 u 和能量 E；
  2. 将方向转换为球坐标 (θ, φ)；
  3. 在 (θ, φ, E) 空间上计算联合概率密度 p(θ, φ, E)，并规范化到“每单位球面面积 per unit energy”；
  4. 保存归一化后的三维数组及边界；
  5. 示例：对固定能量层绘制 Mollweide 投影热图。
"""

import os
import argparse
import h5py
import numpy as np
import matplotlib.pyplot as plt


def load_surface_source(file_name: str, dataset: str = 'source_bank') -> np.ndarray:
    """
    读取 HDF5 文件中 '/source_bank' 复合数据集，返回结构化数组。
    """
    with h5py.File(file_name, 'r') as f:
        return f[dataset][()]


def extract_u_E(data: np.ndarray) -> tuple:
    """
    从结构化数组中提取方向 u (Nx3) 和能量 E (N,)。
    """
    ux = data['u']['x'] if 'u' in data.dtype.names else data['u.x']
    uy = data['u']['y'] if 'u' in data.dtype.names else data['u.y']
    uz = data['u']['z'] if 'u' in data.dtype.names else data['u.z']
    u = np.vstack([ux, uy, uz]).T.astype(np.float64)
    E = data['E'].astype(np.float64)
    return u, E


def compute_joint_pdf(u: np.ndarray,
                      E: np.ndarray,
                      N_theta: int,
                      N_phi: int,
                      N_E: int) -> tuple:
    """
    计算联合密度 H(θ,φ,E) 并转换到每单位球面面积 per unit energy。
    返回 H_solid, theta_centers, phi_centers, E_centers
    """
    # 球坐标
    theta = np.arccos(u[:, 2])               # [0, π]
    phi = np.arctan2(u[:, 1], u[:, 0])       # [-π, π]

    # 定义 bins
    theta_edges = np.linspace(0, np.pi, N_theta + 1)
    phi_edges = np.linspace(-np.pi, np.pi, N_phi + 1)
    E_edges = np.linspace(E.min(), E.max(), N_E + 1)

    # 三维直方图，density=True 返回每单位 (θ,φ,E) 体积的密度
    samples = np.column_stack([theta, phi, E])
    H, edges = np.histogramdd(
        samples,
        bins=[theta_edges, phi_edges, E_edges],
        density=True
    )

    # 中心坐标
    theta_centers = 0.5 * (theta_edges[:-1] + theta_edges[1:])
    phi_centers = 0.5 * (phi_edges[:-1] + phi_edges[1:])
    E_centers = 0.5 * (E_edges[:-1] + E_edges[1:])

    # 转换到 per unit solid angle per unit energy: pΩE = H / sinθ
    # sinθ 可能为 0，添加小常数避免除零
    sin_theta = np.sin(theta_centers)
    sin_theta[sin_theta == 0] = 1e-12
    H_solid = H / sin_theta[:, None, None]

    return H_solid, theta_centers, phi_centers, E_centers


def save_results(output_dir: str,
                 H_solid: np.ndarray,
                 theta_centers: np.ndarray,
                 phi_centers: np.ndarray,
                 E_centers: np.ndarray):
    """
    保存归一化后的 H_solid 和边界到 .npz
    """
    os.makedirs(output_dir, exist_ok=True)
    path = os.path.join(output_dir, 'joint_pdf.npz')
    np.savez_compressed(
        path,
        H_solid=H_solid,
        theta=theta_centers,
        phi=phi_centers,
        E=E_centers
    )
    print(f"Saved joint PDF data to {path}")
    return path


def plot_mollweide_slice(H_solid: np.ndarray,
                         theta_centers: np.ndarray,
                         phi_centers: np.ndarray,
                         E_centers: np.ndarray,
                         k: int,
                         output_dir: str):
    """
    对索引 k 对应的能量层绘制 Mollweide 投影热图。
    """
    PHI, THETA = np.meshgrid(phi_centers, theta_centers)
    lon = PHI
    lat = np.pi/2 - THETA
    H_k = H_solid[:, :, k]

    fig = plt.figure(figsize=(6,4))
    ax = fig.add_subplot(111, projection='mollweide')
    pcm = ax.pcolormesh(lon, lat, H_k, shading='auto')
    cbar = fig.colorbar(pcm, ax=ax, orientation='horizontal', pad=0.05)
    cbar.set_label(f'p(θ,φ|E≈{E_centers[k]:.3g}) [per sr per unit E]')
    ax.set_title('Joint PDF slice on sphere')

    os.makedirs(output_dir, exist_ok=True)
    path = os.path.join(output_dir, f'mollweide_slice_E{k}.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved Mollweide slice k={k} to {path}")


def main():
    parser = argparse.ArgumentParser(
        description='计算球面方向与能量的联合概率密度，并进行可视化')
    parser.add_argument('file', help='HDF5 文件路径，如 surface_source.h5')
    parser.add_argument('--dataset', default='source_bank', help='数据集名称')
    parser.add_argument('--nt', type=int, default=60, help='θ bin 数量')
    parser.add_argument('--np', type=int, default=120, help='φ bin 数量')
    parser.add_argument('--ne', type=int, default=50, help='E bin 数量')
    parser.add_argument('--outdir', default='joint_stats', help='输出目录')
    args = parser.parse_args()

    data = load_surface_source(args.file, args.dataset)
    u, E = extract_u_E(data)

    H_solid, theta_c, phi_c, E_c = compute_joint_pdf(
        u, E, args.nt, args.np, args.ne
    )
    save_results(args.outdir, H_solid, theta_c, phi_c, E_c)

    # 示例：绘制中间能量层
    k0 = len(E_c) // 2
    plot_mollweide_slice(H_solid, theta_c, phi_c, E_c, k0, args.outdir)

if __name__ == '__main__':
    main()
