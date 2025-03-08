import numpy as np
import matplotlib.pyplot as plt

# 设置参数
lambda_param = 1  # 指数分布的速率参数，均值为1/lambda_param = 1
sample_sizes = np.arange(1, 10001)  # 样本数量从1到10000

# 生成随机样本
samples = np.random.exponential(1/lambda_param, (10000, 1))

# 计算不同样本大小下的均值
sample_means = np.cumsum(samples) / sample_sizes

# 绘制图形
plt.figure(figsize=(10, 6))
plt.plot(sample_sizes, sample_means, label='Sample Mean')
plt.axhline(1/lambda_param, color='r', linestyle='dashed', label='True Mean (1/lambda)')
plt.xlabel('Sample Size')
plt.ylabel('Sample Mean')
plt.title('Law of Large Numbers with Exponential Distribution')
plt.legend()
plt.show()
