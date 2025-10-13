# 章节 2：源代码组织与命名空间策略

本章节介绍 OpenMC C++ 源码的目录划分、命名空间与模块边界，帮助读者在阅读类和结构体实现时快速定位文件。

## 目录结构概览

- `include/openmc/`：面向外部组件的头文件接口，声明所有公共类与结构体。
- `src/`：与上述头文件对应的实现文件，包括几何求交、物理事件流程、并行通信等。
- `tallies/`、`physics/` 子目录：将统计与物理核逻辑进一步模块化。

## 命名空间与模块划分

- 全部核心类型均位于 `openmc` 命名空间，避免符号冲突并清晰表达模块所有权。
- 几何模块：`cell.h`、`surface.h`、`universe.h`、`lattice.h` 定义多层次几何实体。【F:include/openmc/cell.h†L58-L305】【F:include/openmc/surface.h†L44-L225】【F:include/openmc/universe.h†L40-L215】【F:include/openmc/lattice.h†L32-L214】
- 粒子与状态模块：`particle_data.h`、`particle.h`、`geometry_state.h` 提供粒子属性与几何上下文。【F:include/openmc/particle_data.h†L20-L257】【F:include/openmc/particle.h†L24-L118】【F:include/openmc/geometry_state.h†L30-L118】
- 材料与核素模块：`material.h`、`nuclide.h`、`reaction.h` 管理核数据与截面。【F:include/openmc/material.h†L42-L265】【F:include/openmc/nuclide.h†L40-L322】【F:include/openmc/reaction.h†L34-L296】

## 成员说明示例

- `openmc::Cell`：成员 `Region region_` 保存几何布尔表达式，`Fill` 联合体记录材料或子宇宙，`material_`、`temperature_` 数组缓存物性信息。【F:include/openmc/cell.h†L124-L305】
- `openmc::Universe`：成员 `std::vector<CellIndex> cells_` 建立单元集合，`bool contains` 系列方法支持层级遍历。【F:include/openmc/universe.h†L112-L215】
- `openmc::ParticleData`：成员 `Position r`、`Direction u`、`GeometryState geom`、`MaterialHandle material` 构成粒子态核心字段。【F:include/openmc/particle_data.h†L20-L257】

## 设计约束

- 模块间通信通过整数索引或轻量句柄完成，减少头文件互相包含。
- 需要跨模块访问的成员通过前向声明和 `friend` 关系控制范围，保持二进制稳定。
