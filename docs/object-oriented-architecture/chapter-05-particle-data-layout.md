# 章节 5：ParticleData 的结构化内存布局

本章节分析 `ParticleData` 的成员字段及其内存布局策略，说明其如何支持高性能粒子输运。

## 核心状态成员

- `Position r`、`Direction u`：粒子在全局坐标系中的位置和方向。
- `double E`：粒子能量。
- `double wgt`：统计权重。
- `int type`：粒子类型枚举（中子、光子等）。
- `GeometryState geom`：嵌套的几何状态。
- `MaterialHandle material`：当前材料句柄。
- `CellInstanceKey cell_instance_id`：唯一标识单元实例，用于 tallies。【F:include/openmc/particle_data.h†L20-L204】

## 截面与缓存成员

- `MacroXS macro_xs`：宏观截面缓存，包含 `total`, `absorption`, `nu_fission` 等字段。
- `MicroXS micro_xs`：微观截面缓存，维护 `nuclide` 序号与截面数组。
- `NuclideMicroXS last_nuclide_xs`：记录上一次计算的核素截面及对应能量。
- `double last_E`、`MaterialHandle last_material`：用于判断是否可复用缓存。
- `double temperature`、`int temperature_index`：当前材料温度及索引。【F:include/openmc/particle_data.h†L167-L257】

## 事件状态与控制成员

- `int generation`、`int num_collisions`：统计粒子生命周期。
- `uint64_t rng_seed`、`uint64_t stream`：随机数生成器状态。
- `bool alive`、`bool banked`：标识粒子是否仍需推进以及是否已写入粒子银行。
- `double time`、`double last_move`：时间相关变量，用于瞬态模拟。
- `int n_secondary`：生成的次级粒子数量，用于协调 `Bank` 写入。【F:include/openmc/particle_data.h†L120-L204】

## 结构化布局策略

- 成员按访问频度排列：高频访问的几何、能量、截面字段靠前，减少缓存行跨越。
- 采用 POD 风格设计，避免虚函数；`Particle` 通过继承复用内存布局而不增加虚表。
- 所有嵌套结构（如 `GeometryState`、`MacroXS`）均内嵌于 `ParticleData`，确保序列化时只需一次内存拷贝。
