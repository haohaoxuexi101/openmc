# 章节 3：几何基础类型概览

本章节讲解 `Surface`、`Region`、`Cell`、`Universe`、`Lattice` 等几何实体及其成员字段，说明它们如何组合描述复杂几何。

## Surface

- `Surface` 为抽象基类，成员 `int32_t id_` 标识表面编号，`BoundaryType bc_` 表示边界条件，`bool sense_` 用于半空间判定。【F:include/openmc/surface.h†L44-L225】
- 派生类型（如 `Plane`, `XPlane`, `Sphere`）扩展几何参数成员，例如 `double x0_`、`double y0_`、`double radius_`，负责距离与法向计算。

## Region

- `Region` 由 `Surface` 的布尔表达式组成，内部成员 `std::vector<SurfaceIndex> surfaces_` 与 `std::vector<int>` 记录布尔运算树，用于 `Cell` 内部空间判定。【F:include/openmc/cell.h†L124-L217】

## Cell

- 关键成员：
  - `int32_t id_`、`std::string name_`：单元标识。
  - `Region region_`：描述几何范围。
  - `Fill fill_`：联合体，可能持有 `MaterialIndex`、`UniverseIndex` 或 `LatticeIndex`。
  - `std::array<double, MAX_TEMPS> sqrtkT_` 与 `std::array<MaterialIndex, MAX_TEMPS>`：缓存温度对应的材料指针。
  - `Position translation_`、`Rotation rotation_`：支持单元的位移与旋转操作。
  - `std::vector<int32_t> neighbor_cells_`：邻接关系，用于加速表面穿越。【F:include/openmc/cell.h†L124-L305】

## Universe

- `int32_t id_`、`std::string name_`：宇宙编号与名称。
- `std::vector<CellIndex> cells_`：包含的单元索引。
- `std::vector<double> weights_`：特定场景下用于均匀抽样单元。
- 方法 `find_cell`、`contains` 通过遍历 `cells_` 决定粒子所在单元。【F:include/openmc/universe.h†L112-L215】

## Lattice

- `int32_t id_` 与 `std::string name_` 用于识别晶格。
- `LatticeType type_` 区分正交或六角晶格。
- `std::array<int, 3> shape_` 记录晶格尺寸，`std::vector<UniverseIndex> universes_` 存储每个格点填充的宇宙。
- `Position origin_`、`std::array<double, 3> pitch_` 记录晶格原点与间距。
- 方法 `indices`, `bounds`, `find_cell` 结合 `LocalCoord` 用于粒子定位。【F:include/openmc/lattice.h†L32-L214】

## 组合关系

- `Cell` 包含 `Region`，并通过 `Fill` 关联 `Material`、`Universe` 或 `Lattice`。
- `Universe` 是 `Cell` 的集合，实现几何层级的递归。
- `Lattice` 在规则网格上引用多个 `Universe`，与 `GeometryState` 的网格索引成员对应。
