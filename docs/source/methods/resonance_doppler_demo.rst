.. _resonance_doppler_demo:

共振多普勒展宽公式与示例代码
==============================

本节以中文重新梳理 :ref:`methods-cross-sections` 中的多极（windowed multipole, WMP）
共振截面公式，并给出一份逐行解释的 Python 脚本，展示如何把理论表达式直接落地为
可执行程序。该脚本不仅计算 0 K 截面和温度展宽后的截面，还会把每个极点的贡献、
近似公式的误差等中间量打印出来，帮助读者更直观地理解多普勒效应对共振结构的影
响。

公式推导回顾
------------

1. **0 K 多极展开**

   对于解析共振区，一个波段内的截面可以写成极点和残量的求和：

   .. math::
      \sigma(E, 0) = \frac{1}{E} \sum_j \operatorname{Re}\left[\frac{i r_j}{\sqrt{E} - p_j}\right],

   其中 :math:`p_j` 是复平面上的极点，:math:`r_j` 是对应的残量，它们均由核数据事
   先拟合得到。推导要点如下：

   * 引入变量 :math:`u = \sqrt{E}`，将 Breit-Wigner 型的能量依赖转化为线性分式；
   * 利用部分分式展开把截面写成极点求和；
   * 通过对共轭项求实部，得到物理上可观测的实截面。

2. **自由气体近似下的多普勒展宽**

   当目标核温度升高时，其热运动导致中子在靶核参考系中的相对能量发生扰动。对自
   由气体速度分布 :math:`M(\mathbf{v}_T)` 进行平均后，截面可写成：

   .. math::
      \sigma(E, T) = \frac{1}{2 E \sqrt{\xi}} \sum_j \operatorname{Re}\left[r_j\sqrt{\pi} W_i(z)
      - \frac{r_j}{\sqrt{\pi}} C\left(\frac{p_j}{\sqrt{\xi}}, \frac{u}{2\sqrt{\xi}}\right)\right],

   其中 :math:`\xi = k_B T / (4A)` 为多普勒展宽宽度，:math:`A` 是靶核质量数，
   :math:`W_i(z)` 与 Faddeeva 函数相关。推导步骤概览：

   * 对相对速度分布积分时，令 :math:`z = (u - p_j)/(2\sqrt{\xi})`，把高斯核与极点项
     结合；
   * 通过引入复误差函数 :math:`W(z)`，把积分解析化；
   * 剩余项形成 :math:`C` 函数，即所谓的 *correction integral*，在能量远大于
     :math:`k_B T/A` 时可以忽略。

3. **高能近似**

   当 :math:`E \gg k_B T / A` 时，多数情况下可以忽略 :math:`C` 积分，得到更简洁的
   近似式：

   .. math::
      \sigma(E, T) \approx \frac{1}{2 E \sqrt{\xi}} \sum_j \operatorname{Re}\left[i r_j \sqrt{\pi} W_i(z)\right].

   在数值实现中，一般同时保留完整式和近似式，以评估忽略 :math:`C` 项带来的误差
   是否可接受。

代码示例讲解
------------

下面的脚本位于 ``examples/python/resonance_doppler_demo.py``，包含如下部分：

#. ``MultipoleDataset``：用 ``dataclass`` 封装极点和残量，并在初始化时确保形状一致；
#. ``zero_kelvin_cross_section``：矢量化实现式 (1)，直接套用 :math:`i r_j / (\sqrt{E}-p_j)`；
#. ``_faddeeva_half_plane``：处理 Faddeeva 函数在不同半平面的支路选择；
#. ``_c_integral``：对 :math:`C` 积分按实部、虚部分开使用 ``scipy.integrate.quad``；
#. ``doppler_broadened_cross_section``：组合式 (2)，可通过 ``include_correction``
   开关保留或忽略 :math:`C` 项；
#. ``debug_pole_contributions``：**新增的调试函数**，逐个极点列出其对截面的贡献与
   多普勒展宽前后的差异；
#. ``tabulate_cross_sections``：给出 0.01–10 eV 的对数能网，比较 0 K、全温度展宽、
   高能近似三者；
#. ``main``：打印详细的中文表格和极点分析结果。

运行脚本的输出示例：

.. code-block:: text

   python examples/python/resonance_doppler_demo.py

   ============================== 基本参数 ==============================
   温度 T = 900.000 K, 靶核质量数 A = 238.000, 多普勒宽度 ξ = 8.150e-05 eV
   能量网格: 8 个点, 范围 [1.000e-02, 1.000e+01] eV

   ------------------------------ 截面表 -------------------------------
        E [eV]   |   σ_0K [b]   |    σ_T [b]    | σ_T (忽略 C) [b] | 相对误差
      1.00000e-02 |  1.234567e+02 |  8.765432e+01 |      8.700000e+01 |  7.47e-03
      ...

   --------------------------- 极点逐项贡献 ---------------------------
   能量 E = 1.00000e-01 eV
     极点 #0: σ_0K = 1.2345e+02 b, σ_T = 8.9123e+01 b, 贡献比例 0.67
     ...

   以上输出展示了温度展宽如何削弱尖锐共振、各极点如何参与叠加。读者可以根据需求，
   把 ``demo_dataset`` 中的参数替换为实际核素（例如 U-238）在对应窗口的多极系数，
   以获得真实截面。

完整代码如下所示（为便于对照，保留丰富注释）：

.. literalinclude:: ../../examples/python/resonance_doppler_demo.py
   :language: python
   :lines: 1-400
   :linenos:

后续扩展建议
------------

* 若需绘制截面曲线，可在脚本末尾使用 ``matplotlib`` 绘制 ``σ(E,0)`` 与 ``σ(E,T)``；
* 通过增大 ``tabulate_cross_sections`` 中的能量点数，可观察高能区 ``C`` 项的影响；
* 对真实模拟而言，OpenMC 会在 C++ 内核中调用相同公式，并通过缓存与向量化加速计
  算，本示例旨在帮助读者理解理论与实现之间的对应关系。
