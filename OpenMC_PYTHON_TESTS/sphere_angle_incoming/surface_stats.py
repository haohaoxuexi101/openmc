#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
surface_stats.py

对 HDF5 面源文件（如 surface_source.h5）的 '/source_bank' 复合数据集：
  1. 提取所有数值字段；
  2. 打印基本统计量（min, max, mean, std, count）；
  3. 绘制并保存每个字段的概率密度直方图；
  4. 在三维空间中可视化能量 E 随位置 r 的分布；
  5. 绘制更清晰的 3D 向量场；
  6. 在固定 z 平面上绘制 2D 方向切片；
  7. 在球面 Mollweide 投影上：
     a) 绘制方向密度热力图；
     b) 绘制方向平均能量热力图。
"""

import os
import argparse
import h5py
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401


def load_surface_source(file_name: str, dataset_name: str = 'source_bank') -> np.ndarray:
    """
    读取 HDF5 文件中 '/source_bank' 复合数据集，返回 NumPy 结构化数组。
    """
    with h5py.File(file_name, 'r') as f:
        return f[dataset_name][()]


def extract_numeric_fields(data: np.ndarray) -> dict:
    """
    从结构化数组中挑出所有数值字段，支持嵌套结构，返回 dict。
    key: 'r.x','r.y','r.z','u.x','u.y','u.z','E',...
    value: 1D np.ndarray
    """
    numeric = {}
    for field in data.dtype.names:
        arr = data[field]
        if arr.dtype.names is not None:
            for sub in arr.dtype.names:
                sub_arr = arr[sub]
                if np.issubdtype(sub_arr.dtype, np.number):
                    numeric[f"{field}.{sub}"] = sub_arr.astype(np.float64).flatten()
        else:
            if np.issubdtype(arr.dtype, np.number):
                numeric[field] = arr.flatten()
    return numeric


def print_stats(numeric: dict):
    for key, arr in numeric.items():
        print(f"Field '{key}': min={arr.min():.5g}, max={arr.max():.5g}, "
              f"mean={arr.mean():.5g}, std={arr.std():.5g}, count={arr.size}")


def plot_histograms(numeric: dict, bins: int, output_dir: str):
    os.makedirs(output_dir, exist_ok=True)
    for key, arr in numeric.items():
        hist, edges = np.histogram(arr, bins=bins, density=True)
        plt.figure(figsize=(8,4))
        plt.bar(edges[:-1], hist, width=np.diff(edges), alpha=0.7, edgecolor='black')
        plt.xlabel(key); plt.ylabel('Probability Density')
        plt.title(f'Distribution of {key}')
        safe = key.replace('.', '_')
        plt.savefig(os.path.join(output_dir, f'{safe}_hist.png'), dpi=150, bbox_inches='tight')
        plt.close()
        print(f"Saved histogram for {key}")


def plot_energy_vs_space(numeric: dict, output_dir: str, sample: int = 20000):
    r = np.vstack((numeric['r.x'], numeric['r.y'], numeric['r.z'])).T
    E = numeric['E']
    idx = np.random.choice(len(E), min(sample, len(E)), replace=False)
    pts, Es = r[idx], E[idx]
    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111, projection='3d')
    sc = ax.scatter(pts[:,0], pts[:,1], pts[:,2], c=Es, cmap='inferno', s=1, alpha=0.7)
    fig.colorbar(sc, ax=ax, label='Energy E')
    ax.set_xlabel('r.x'); ax.set_ylabel('r.y'); ax.set_zlabel('r.z')
    ax.set_title('Energy vs Space')
    path = os.path.join(output_dir, 'energy_space.png')
    fig.savefig(path, dpi=150, bbox_inches='tight'); plt.close(fig)
    print(f"Saved energy-space scatter to {path}")


def plot_direction_field_3d(numeric: dict, output_dir: str, sample: int = 3000,
                            length: float = 5, linewidth: float = 1.2):
    r = np.vstack((numeric['r.x'], numeric['r.y'], numeric['r.z'])).T
    u = np.vstack((numeric['u.x'], numeric['u.y'], numeric['u.z'])).T
    idx = np.random.choice(len(u), min(sample, len(u)), replace=False)
    rs, us = r[idx], u[idx]
    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111, projection='3d')
    ax.quiver(rs[:,0], rs[:,1], rs[:,2], us[:,0], us[:,1], us[:,2],
              length=length, normalize=True, linewidth=linewidth)
    ax.set_xlabel('r.x'); ax.set_ylabel('r.y'); ax.set_zlabel('r.z')
    ax.set_title('3D Direction Field')
    path = os.path.join(output_dir, 'dir_field_3d.png')
    fig.savefig(path, dpi=150, bbox_inches='tight'); plt.close(fig)
    print(f"Saved 3D direction field to {path}")


def plot_direction_slice_2d(numeric: dict, output_dir: str,
                            sample: int = 1500, z_tol: float = 105.005):
    r = np.vstack((numeric['r.x'], numeric['r.y'], numeric['r.z'])).T
    u = np.vstack((numeric['u.x'], numeric['u.y'], numeric['u.z'])).T
    zc = np.median(r[:,2])
    mask = np.abs(r[:,2] - zc) < z_tol
    rs, us = r[mask], u[mask]
    if len(rs)==0:
        print('No slice points'); return
    idx = np.random.choice(len(rs), min(sample,len(rs)), replace=False)
    rs, us = rs[idx], us[idx]
    fig, ax = plt.subplots(figsize=(6,6))
    ax.quiver(rs[:,0], rs[:,1], us[:,0], us[:,1], angles='xy', scale_units='xy', scale=0.03, width=0.003)
    ax.set_aspect('equal'); ax.set_xlabel('r.x'); ax.set_ylabel('r.y')
    ax.set_title(f'Direction Slice z≈{zc:.3f}')
    path = os.path.join(output_dir, 'dir_slice_2d.png')
    fig.savefig(path, dpi=150, bbox_inches='tight'); plt.close(fig)
    print(f"Saved 2D direction slice to {path}")

def plot_spherical_mollweide(numeric: dict, output_dir: str,
                             n_theta: int = 15, n_phi: int = 25):
    import os
    import numpy as np
    import matplotlib.pyplot as plt

    # 1) 方向向量 & 能量
    u = np.vstack((numeric['u.x'], numeric['u.y'], numeric['u.z'])).T
    E = numeric['E']

    # 2) 转成球面坐标 θ∈[0,π], φ∈[-π,π]
    theta = np.arccos(np.clip(u[:, 2], -1.0, 1.0))
    phi   = np.arctan2(u[:, 1], u[:, 0])

    # 3) 统计每个 bin 的 counts 和 sumE
    counts, th_edges, ph_edges = np.histogram2d(
        theta, phi,
        bins=[n_theta, n_phi],
        density=False
    )
    sumE, _, _ = np.histogram2d(
        theta, phi,
        bins=[th_edges, ph_edges],
        weights=E
    )
    meanE = sumE / np.maximum(counts, 1)

    # 4) 用 edges 构造网格（shape = (n_theta+1, n_phi+1)）
    TH_edges, PH_edges = np.meshgrid(th_edges, ph_edges, indexing='ij')
    # Mollweide 需要 lon=φ 边界，lat=π/2−θ 边界
    lon = PH_edges            # shape (n_theta+1, n_phi+1)
    lat = np.pi/2 - TH_edges

    # 5) 绘图
    os.makedirs(output_dir, exist_ok=True)
    fig = plt.figure(figsize=(10, 5))

    ax1 = fig.add_subplot(121, projection='mollweide')
    pcm1 = ax1.pcolormesh(lon, lat, counts, shading='auto')
    ax1.set_title('Direction Counts')
    fig.colorbar(pcm1, ax=ax1, orientation='horizontal', pad=0.05)

    ax2 = fig.add_subplot(122, projection='mollweide')
    pcm2 = ax2.pcolormesh(lon, lat, meanE, shading='auto')
    ax2.set_title('Mean Energy ⟨E⟩ per Bin')
    fig.colorbar(pcm2, ax=ax2, orientation='horizontal', pad=0.05)

    out = os.path.join(output_dir, 'spherical_stats_counts.png')
    fig.savefig(out, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved spherical counts & meanE to {out}")




# def plot_spherical_mollweide(numeric: dict, output_dir: str,
#                              n_theta: int = 15, n_phi: int = 25):
#     u = np.vstack((numeric['u.x'], numeric['u.y'], numeric['u.z'])).T
#     E = numeric['E']
#     # 球坐标
#     theta = np.arccos(u[:,2])
#     phi = np.arctan2(u[:,1], u[:,0])
#     # 密度直方图
#     H_dir, th_edges, ph_edges = np.histogram2d(theta, phi, bins=[n_theta, n_phi], density=True)
#     # 能量加权平均
#     # 累积能量
#     sumE, _, _ = np.histogram2d(theta, phi, bins=[n_theta, n_phi], weights=E)
#     meanE = sumE / (H_dir * np.sum(H_dir) * (th_edges[1]-th_edges[0]) * (ph_edges[1]-ph_edges[0]) + 1e-12)
#     # 投影网格
#     phc = (ph_edges[:-1]+ph_edges[1:])/2
#     thc = (th_edges[:-1]+th_edges[1:])/2
#     PHI, THETA = np.meshgrid(phc, thc)
#     lon = PHI; lat = np.pi/2 - THETA
#     # 绘图
#     fig = plt.figure(figsize=(10,5))
#     ax1 = fig.add_subplot(121, projection='mollweide')
#     pcm1 = ax1.pcolormesh(lon, lat, H_dir, shading='auto')
#     ax1.set_title('Direction Density')
#     fig.colorbar(pcm1, ax=ax1, orientation='horizontal', pad=0.05)
#     ax2 = fig.add_subplot(122, projection='mollweide')
#     pcm2 = ax2.pcolormesh(lon, lat, meanE, shading='auto')
#     ax2.set_title('Mean Energy per Direction')
#     fig.colorbar(pcm2, ax=ax2, orientation='horizontal', pad=0.05)
#     os.makedirs(output_dir, exist_ok=True)
#     path = os.path.join(output_dir, 'spherical_stats.png')
#     fig.savefig(path, dpi=150, bbox_inches='tight'); plt.close(fig)
#     print(f"Saved spherical stats to {path}")


def main():
    p = argparse.ArgumentParser(description="面源统计 & 球面可视化")
    p.add_argument('file'); p.add_argument('--dataset', default='source_bank')
    p.add_argument('--bins', type=int, default=100); p.add_argument('--outdir', default='stats')
    args = p.parse_args()
    data = load_surface_source(args.file, args.dataset)
    num = extract_numeric_fields(data)
    os.makedirs(args.outdir, exist_ok=True)
    print_stats(num)
    plot_histograms(num, args.bins, args.outdir)
    plot_energy_vs_space(num, args.outdir)
    plot_direction_field_3d(num, args.outdir)
    plot_direction_slice_2d(num, args.outdir)
    plot_spherical_mollweide(num, args.outdir)
    print('All visualizations saved in', args.outdir)

if __name__=='__main__':
    main()
