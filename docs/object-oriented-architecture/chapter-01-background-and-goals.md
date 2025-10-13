# 章节 1：文档背景与目标

本章节阐述说明书的定位与覆盖范围，明确对象是 OpenMC C++ 内核中的几何、粒子、材料及统计模块，尤其聚焦 `Particle`、`ParticleData`、`GeometryState` 与 `LocalCoord` 的数据流关系。文档旨在为架构审查、培训课程与功能拓展提供可独立研读的材料，每个章节都包含类与结构体的成员说明，帮助读者迅速映射源码结构。

## 关键类型索引

- `Particle`：粒子事件接口的核心类型，定义于 `include/openmc/particle.h`。【F:include/openmc/particle.h†L24-L118】
- `ParticleData`：粒子状态数据容器，定义于 `include/openmc/particle_data.h`。【F:include/openmc/particle_data.h†L20-L257】
- `GeometryState` 与 `LocalCoord`：负责几何层级跟踪，分别位于 `include/openmc/geometry_state.h` 与 `include/openmc/particle_data.h`。【F:include/openmc/geometry_state.h†L30-L118】【F:include/openmc/particle_data.h†L80-L175】
- 几何与材料实体：`Cell`、`Universe`、`Lattice`、`Surface`、`Material` 等，均在 `include/openmc` 目录中提供头文件接口。【F:include/openmc/cell.h†L58-L305】【F:include/openmc/universe.h†L40-L215】【F:include/openmc/lattice.h†L32-L214】【F:include/openmc/surface.h†L44-L225】【F:include/openmc/material.h†L42-L265】

## 章节目标

1. 描述面向对象层次及组合关系，帮助理解模块耦合方式。
2. 分析关键类成员的作用与数据流，强调缓存策略与性能考量。
3. 梳理与物理核、统计、并行化之间的接口设计，为后续章节做铺垫。
