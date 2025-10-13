# 章节 9：Surface 交互与边界条件

本章节解释粒子穿越表面时涉及的类与成员，展示几何边界如何与粒子事件结合。

## Surface 成员回顾

- `BoundaryType bc_`：边界条件（真空、反射、周期）。
- `bool sense_`：定义正向半空间。
- `virtual double distance(Position r, Direction u)`：计算到表面的距离。
- `virtual Direction normal(Position r)`：返回表面法向。
- `virtual void surface_interaction(Particle& p)`：部分派生类重写以处理边界特殊逻辑。【F:include/openmc/surface.h†L44-L225】【F:src/surface.cpp†L34-L379】

## Particle 与 Surface 的接口

- `Particle::event_advance` 调用 `distance_to_boundary`（`geometry.cpp` 中实现）求最近表面，并将距离最小的表面索引写入 `GeometryState::surface_last`。【F:src/particle.cpp†L214-L349】【F:src/geometry.cpp†L411-L587】
- `Particle::event_cross_surface` 根据 `Surface` 的 `bc_` 成员：
  - 若 `VACUUM`：将 `alive` 设为 `false`。
  - 若 `REFLECTIVE`：调用 `Surface::reflect` 调整 `Direction` 并更新 `GeometryState::on_boundary`。
  - 若 `PERIODIC`：调用 `Surface::periodic_translate` 修改 `Position r`。

## GeometryState 的作用

- `surface_last`：记录最后穿越的表面，供 tallies 与反射逻辑查询。
- `on_boundary`：粒子恰好位于边界时置为 `true`，避免重复计算距离。
- `level_last`：穿越表面并需要回到父宇宙时使用。【F:include/openmc/geometry_state.h†L30-L118】

## 设计要点

- `Surface` 通过虚函数实现多态，但 `Particle` 侧只需持有索引；查找具体对象由全局数组完成。
- 边界条件的状态通过 `GeometryState` 成员暴露给统计模块，实现如表面通量 tally。
