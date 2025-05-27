#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
surface_output_stats.py

统计单球输出特征（位置、方向、能量）并在球面位置上绘制：
  - 平均能量
  - 方向分量（径向、极向、方位分量）

方法：
  1. 提取位置 r, 方向 u, 能量 E
  2. 位置转换为球坐标 (θ_r, φ_r)
  3. 统计每个球面分区的：
     - 平均能量 meanE
     - 平均径向分量 mean_u_radial
     - 平均极向分量 mean_u_theta
     - 平均方位分量 mean_u_phi
  4. 可视化：四个 Mollweide 热力图
  5. 保存统计结果至 .npz 以供后续快速采样

默认使用较少分区：θ=30, φ=60，可通过参数调整。
"""

import os
import argparse
import h5py
import numpy as np
import matplotlib.pyplot as plt


def load_data(file_path: str, dataset: str = 'source_bank') -> np.ndarray:
    with h5py.File(file_path, 'r') as f:
        return f[dataset][()]


def extract_r_u_E(data: np.ndarray):
    def get_field(prefix, comp):
        if prefix in data.dtype.names and data[prefix].dtype.names:
            return data[prefix][comp]
        return data[f"{prefix}.{comp}"]

    rx = get_field('r','x'); ry = get_field('r','y'); rz = get_field('r','z')
    ux = get_field('u','x'); uy = get_field('u','y'); uz = get_field('u','z')
    r = np.vstack([rx, ry, rz]).T.astype(np.float64)
    u = np.vstack([ux, uy, uz]).T.astype(np.float64)
    E = data['E'].astype(np.float64)
    return r, u, E


def spherical_coords(v: np.ndarray):
    norm = np.linalg.norm(v, axis=1)
    x, y, z = v[:,0]/norm, v[:,1]/norm, v[:,2]/norm
    theta = np.arccos(np.clip(z, -1.0, 1.0))
    phi = np.arctan2(y, x)
    return theta, phi


def compute_stats(theta, phi, E, u, bins_theta, bins_phi):
    edges_t = np.linspace(0, np.pi, bins_theta+1)
    edges_p = np.linspace(-np.pi, np.pi, bins_phi+1)

    H   , _, _ = np.histogram2d(theta, phi, bins=[edges_t, edges_p])
    sumE, _, _ = np.histogram2d(theta, phi, bins=[edges_t, edges_p], weights=E)

    n = np.column_stack([np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)])
    e_theta = np.column_stack([np.cos(theta)*np.cos(phi), np.cos(theta)*np.sin(phi), -np.sin(theta)])
    e_phi   = np.column_stack([-np.sin(phi), np.cos(phi), np.zeros_like(phi)])

    u_rad = np.einsum('ij,ij->i', u, n)
    u_th  = np.einsum('ij,ij->i', u, e_theta)
    u_ph  = np.einsum('ij,ij->i', u, e_phi)

    sum_ur, _, _ = np.histogram2d(theta, phi, bins=[edges_t, edges_p], weights=u_rad)
    sum_ut, _, _ = np.histogram2d(theta, phi, bins=[edges_t, edges_p], weights=u_th)
    sum_up, _, _ = np.histogram2d(theta, phi, bins=[edges_t, edges_p], weights=u_ph)

    tc = 0.5*(edges_t[:-1]+edges_t[1:])
    pc = 0.5*(edges_p[:-1]+edges_p[1:])

    meanE    = np.divide(sumE, H, where=H>0)
    mean_ur  = np.divide(sum_ur, H, where=H>0)
    mean_ut  = np.divide(sum_ut, H, where=H>0)
    mean_up  = np.divide(sum_up, H, where=H>0)

    return H, meanE, mean_ur, mean_ut, mean_up, tc, pc


def plot_mollweide(data2d, theta_centers, phi_centers, title, fname, outdir):
    PHI, THETA = np.meshgrid(phi_centers, theta_centers)
    lon = PHI
    lat = np.pi/2 - THETA
    fig = plt.figure(figsize=(6,4))
    ax = fig.add_subplot(111, projection='mollweide')
    pcm = ax.pcolormesh(lon, lat, data2d, shading='auto')
    fig.colorbar(pcm, ax=ax, orientation='horizontal', pad=0.05)
    ax.set_title(title)
    path = os.path.join(outdir, fname)
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved {fname}")


def main():
    parser = argparse.ArgumentParser(description='球面位置输出特征统计')
    parser.add_argument('file', help='HDF5 文件路径，例如 surface_source.h5')
    parser.add_argument('--dataset', default='source_bank', help='数据集名称')
    parser.add_argument('--nt', type=int, default=30, help='θ 分区数量 (默认 30)')
    parser.add_argument('--np', type=int, default=60, help='φ 分区数量 (默认 60)')
    parser.add_argument('--outdir', default='sphere_pos_stats', help='输出目录')
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    data = load_data(args.file, args.dataset)
    r, u, E = extract_r_u_E(data)
    theta_r, phi_r = spherical_coords(r)

    H, mE, mur, mut, mup, tc, pc = compute_stats(theta_r, phi_r, E, u, args.nt, args.np)

    np.savez_compressed(os.path.join(args.outdir,'sphere_pos_stats.npz'),
                        count=H, meanE=mE, mean_u_radial=mur,
                        mean_u_theta=mut, mean_u_phi=mup,
                        theta=tc, phi=pc)

    plot_mollweide(mE,  tc, pc, 'Mean Energy on Sphere','mean_energy.png',   args.outdir)
    plot_mollweide(mur, tc, pc, 'Mean Radial Component','mean_u_radial.png', args.outdir)
    plot_mollweide(mut, tc, pc, 'Mean Theta Component','mean_u_theta.png',  args.outdir)
    plot_mollweide(mup, tc, pc, 'Mean Phi Component','mean_u_phi.png',      args.outdir)

    print('All sphere position stats saved to', args.outdir)

if __name__=='__main__':
    main()