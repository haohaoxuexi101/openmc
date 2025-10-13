# 章节 13：接口扩展与设计约束总结

本章节总结对象关系与成员设计的约束，提出扩展新功能时的注意事项。

## 继承与组合策略

- 仅在需要运行时多态的 `Surface`、`Filter`、`Distribution` 等类型上使用虚函数；其他类型（如 `ParticleData`、`GeometryState`）坚持组合与 POD 布局，以降低内存开销。【F:include/openmc/surface.h†L44-L225】【F:include/openmc/tallies/filter.h†L38-L312】【F:include/openmc/particle_data.h†L20-L257】
- `Particle` 通过继承共享成员布局，但不新增虚函数，确保可以 `memcpy`。

## 成员演化注意事项

- 扩展 `ParticleData` 时，需要同步更新 `BankedParticle`、`SourceSite`、`StatePoint` 序列化逻辑，避免数据错位。【F:include/openmc/bank.h†L31-L205】【F:include/openmc/source.h†L33-L152】【F:src/state_point.cpp†L52-L410】
- 修改 `GeometryState::coord` 长度需检查 `MAX_COORD` 相关常量，以及几何遍历算法在 `geometry.cpp` 中的边界处理。【F:include/openmc/geometry_state.h†L30-L118】【F:src/geometry.cpp†L210-L587】
- 新增材料属性时，要考虑 `Material::calculate_xs` 与 `Nuclide` 接口的兼容性，确保截面缓存成员依然有效。【F:include/openmc/material.h†L170-L265】【F:include/openmc/nuclide.h†L180-L322】

## 接口扩展示例

- **新增几何类型**：继承 `Surface` 并实现 `distance` 与 `normal`；更新 `geometry.cpp` 的几何工厂以支持新类型。
- **新增材料模型**：扩展 `Material` 成员（如温度依赖表），同时修改 `ParticleData` 的缓存失效逻辑以容纳新状态。
- **新增粒子种类**：扩展 `ParticleType` 枚举、`Particle::event_collide` 分支，并更新 tallies 与输出以支持新的成员字段。【F:include/openmc/particle.h†L70-L118】【F:src/particle.cpp†L176-L412】

## 汇总

- 各章节列出的成员关系构成粒子输运仿真的基础链路：`ParticleData` -> `GeometryState` -> `Cell`/`Lattice`/`Universe` -> `Material`/`Nuclide` -> `Tally`/并行服务。
- 扩展功能时，应优先遵循现有成员命名与布局策略，确保与 tallies、银行、并行通信兼容。
