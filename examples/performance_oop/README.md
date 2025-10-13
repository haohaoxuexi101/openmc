# 面向对象设计与性能对比实验

本示例基于 OpenMC 中 `GeometryState → ParticleData → Particle` 的继承结构，解释为什么几何状态、粒子属性和传输行为被组织在同一个对象层次中，以及这样设计对性能和教学的帮助。本目录包含三种实现：

1. **`polymorphic`**：与 OpenMC 类似的继承式设计，所有状态都位于同一对象中，算法只保留一个指针即可访问几何和粒子属性。
2. **`aggregate`**：使用聚合对象分别保存几何和粒子数据，通过间接访问组合对象；每个组件单独分配在堆上，代表常见的“面向过程 + 组合”写法，存在额外的指针跳转和缓存失配。
3. **`soa`**：面向数据的结构化数组（Structure of Arrays）写法，对批量更新非常高效，用于对比继承式设计在保持封装的同时仍能接近数据导向写法的性能。

`benchmark.cpp` 将三种实现放在同一个可执行程序中，通过 `std::chrono::steady_clock` 测量 1,000,000 个粒子在 200 个传输步长下的耗时，并验证三种写法得到的物理量保持一致。

## 为什么要采用继承层次？

* **缓存局部性**：`Particle` 继承 `ParticleData` 和 `GeometryState`，所有数据连续排布。运行传输算法时只需一次缓存加载即可获得所有成员，减少聚合/指针方式产生的随机访问成本。
* **接口复用**：OpenMC 的许多算法只关心几何状态（例如判定网格单元）或粒子通量信息。有了继承层次，我们可以将算法的参数写成基类引用，实现重用并保持易读性。
* **教学友好**：通过一个对象就能演示“状态 + 行为”的封装，同时保留与真实 Monte Carlo 代码相似的接口，便于讲解多态、继承和组合的取舍。

## 编译与运行

```bash
# 在仓库根目录
c++ -std=c++17 -O3 examples/performance_oop/benchmark.cpp -o /tmp/perf_demo
/tmp/perf_demo
```

程序运行后会输出每种实现的耗时（毫秒）以及一次性校验结果，帮助教师实时演示不同写法的性能差异。

## 输出示例

```
Scenario: polymorphic  | Time: 428.24  ms | Final energy sum: 3.74e-29
Scenario: aggregate    | Time: 533.19  ms | Final energy sum: 3.74e-29
Scenario: soa          | Time: 323.87  ms | Final energy sum: 3.74e-29
All strategies produce equivalent physics results.
```

> 注：实际时间依赖于硬件与编译器版本，但相对趋势能够稳定反映继承式设计减少指针间接、聚合式设计额外开销的特点。

## 性能观察

* 继承式 `polymorphic` 写法将几何状态与粒子属性合并在同一对象内，加载和更新时只需一次指针解引用，CPU 缓存命中率最高。
* 聚合式 `aggregate` 写法每个粒子在堆上分配两个独立对象，遍历时需要经过 `AggregateParticle` → `geom/data` 的跳转，额外的缓存失效使得耗时增长约 20–30%。
* 数据导向的 `soa` 在批量操作下最为高效，但牺牲了封装性与可扩展性；通过对比可以强调 OpenMC 的继承式设计在可维护性与性能之间取得折衷。
