# 章节 6：Particle 行为层与事件生命周期

本章节说明 `Particle` 类如何在 `ParticleData` 的基础上定义事件接口，并列出关键成员与方法。

## 继承关系

- `class Particle : public ParticleData`：继承使 `Particle` 拥有所有状态成员，同时增加行为方法，不引入虚函数表。【F:include/openmc/particle.h†L24-L118】

## 关键成员与方法

- `Particle()` 构造函数：初始化 `ParticleData` 字段并设置默认随机数种子。
- `void event_calculate_xs()`：调用 `physics::calculate_xs` 更新 `macro_xs` 与 `micro_xs` 缓存。【F:src/particle.cpp†L38-L175】
- `void event_advance()`：推进粒子至最近的碰撞或边界，更新 `GeometryState::surface_last`。
- `void event_cross_surface(SurfaceIndex surf)`：处理跨越表面后的边界条件与层级更新。【F:src/particle.cpp†L350-L611】
- `void event_collide()`：根据当前材料与截面采样反应，可能生成次级粒子。【F:src/particle.cpp†L176-L349】
- `void create_secondary(...)`：将次级粒子写入银行，使用 `ParticleData` 的 `n_secondary` 与 `secondary_bank_` 字段协调。【F:include/openmc/particle.h†L80-L118】
- `RNG rng()`：基于 `rng_seed` 与 `stream` 生成随机数，保证事件采样的可重现性。【F:include/openmc/particle.h†L40-L68】

## 成员交互

- `event_calculate_xs` 读取 `material` 与 `E`，写入 `macro_xs` 与 `micro_xs`。
- `event_advance` 使用 `GeometryState` 栈定位当前 `Cell`，调用几何函数求最近表面距离。
- `event_collide` 通过 `Material` 内的核素列表和 `NuclideMicroXS` 生成碰撞产物，并更新 `num_collisions`。

## 设计考量

- 行为方法保持无异常抛出路径，便于在并行环境中使用。
- 与 `ParticleData` 的 POD 布局结合，支持在 MPI 通信时直接序列化整个对象。
