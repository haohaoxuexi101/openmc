import numpy as np
import matplotlib.pyplot as plt

# 设置参数
num_samples = 10000  # 每次模拟的样本数
sample_size = 30  # 每次抽取的样本大小

# 生成样本均值
sample_means = []
for _ in range(num_samples):
    sample = np.random.uniform(0, 1, sample_size)  # 从均匀分布抽样
    sample_means.append(np.mean(sample))

# 绘制图形
plt.figure(figsize=(10, 6))
plt.hist(sample_means, bins=100, density=True, alpha=0.6, color='g')
plt.title('Central Limit Theorem')
plt.xlabel('Sample Mean')
plt.ylabel('Density')
plt.show()
