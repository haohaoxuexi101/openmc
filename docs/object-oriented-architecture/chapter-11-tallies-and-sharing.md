# 章节 11：统计收集、Tallies 与粒子态共享

本章节阐述统计模块与粒子状态之间的成员交互，说明 tallies 如何利用粒子成员数据。

## Tally 成员

- `int32_t id_`、`std::string name_`：统计量标识。
- `Filters filters_`：包含 `Filter` 派生对象指针，决定粒子属性筛选条件。
- `std::vector<double> results_`：存储累积结果。
- `int estimator_`：统计估计器类型（tracklength、collision 等）。
- `bool active_`：是否启用统计。【F:include/openmc/tallies/tally.h†L40-L266】

## Filter 成员示例

- `CellFilter::bins_`：感兴趣的 `Cell` 索引集合。
- `SurfaceFilter::bins_`：用于表面通量统计的表面编号列表。
- `MaterialFilter::bins_`：材料索引集合。
- `EnergyFilter::edges_`：能量分箱边界。【F:include/openmc/tallies/filter.h†L38-L312】

## 与 ParticleData 的接口

- `score_general_tallies` 从 `ParticleData` 读取 `cell_instance_id`、`material`、`E`、`wgt`、`geom.surface_last` 等成员进行匹配与累积。【F:src/tallies/tally.cpp†L72-L491】【F:include/openmc/particle_data.h†L120-L257】
- 对于表面 tallies，`GeometryState::surface_last` 与 `on_boundary` 提供必要的上下文。
- `ParticleData::wgt` 与 `time` 可用于功率或瞬态统计。

## 粒子态共享结构

- `BankedParticle`：与 `ParticleData` 部分字段保持相同布局（位置、方向、能量、权重），用于生成源粒子列表。【F:include/openmc/bank.h†L31-L205】
- `SourceSite`：在源迭代与重启文件中保存粒子状态，同样复用 `ParticleData` 的核心字段。【F:include/openmc/source.h†L33-L152】

## 设计特点

- 通过共享成员布局，统计与源重建均可直接拷贝 `ParticleData`。
- Tallies 对粒子成员的读取不需要锁机制，利用线程局部缓存或批量归并实现并行安全。
