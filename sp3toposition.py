from datetime import datetime
import numpy as np
import bisect
import matplotlib.pyplot as plt
from astropy.coordinates import ITRS,GCRS
from astropy.time import Time
from astropy import units as u

plt.rcParams['font.sans-serif'] = ['SimHei']  # 设置中文字体
plt.rcParams['axes.unicode_minus'] = False  # 显示负号

def get_data(sp3_data, target_satellite):# 获取数据
    matched_data = []  
    position = [[], [], [], []]  # X, Y, Z, 时间

    current_time = None
    with open(sp3_data, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith("*"):#获取时间
                parts = line.split()
                year = int(parts[1])
                month = int(parts[2])
                day = int(parts[3])
                hour = int(parts[4])
                minute = int(parts[5])
                second = int(float(parts[6]))
                current_time = datetime(year, month, day, hour, minute, second)
                continue
            
            elif line.startswith(target_satellite):#获取目标卫星坐标
                if current_time is not None:
                    parts = line.split()

                    x = float(parts[1])
                    y = float(parts[2])
                    z = float(parts[3])

                    matched_data.append((x, y, z, current_time))

    for x, y, z, time in matched_data:
        position[0].append(x)
        position[1].append(y)
        position[2].append(z)
        position[3].append(time)

    print(f"时间范围: {min(position[3])} 至 {max(position[3])}")# 输出时间范围
    return position

def lagrange(position, target_time, num=8):# 拉格朗日插值法
    if not position[3]:# 如果没有时间数据
        raise ValueError("未获取时间")
    
    timestamps = np.array([t.timestamp() for t in position[3]])# 时间戳
    target_x = target_time.timestamp()
    
    idx = bisect.bisect_left(timestamps, target_x)# 查找最近的时间索引
    
    # 窗口范围
    start = max(0, idx - num//2)#大于等于0
    end = min(len(timestamps), idx + num//2)
    
    # 如果窗口不足
    if end - start < num:
        if start == 0:
            end = min(len(timestamps), start + num)
        else:
            start = max(0, end - num)
    
    window_times = timestamps[start:end]
    # 窗口数据
    window_data = [
        np.array(position[0][start:end]),
        np.array(position[1][start:end]),
        np.array(position[2][start:end])
    ]
    
    # 拉格朗日插值
    target_location = []
    for coord in window_data:
        n = len(coord)
        result = 0.0
        for j in range(n):
            numerator = 1.0
            denominator = 1.0
            for k in range(n):
                if k != j:
                    numerator *= (target_x - window_times[k])
                    denominator *= (window_times[j] - window_times[k])
            if denominator == 0:
                raise ZeroDivisionError("时间戳重复导致分母为零")
            result += coord[j] * (numerator / denominator)
        target_location.append(result)
    
    return target_location

# 运行
target_satellite = "P" + input("输入卫星编号(G02-G32:) ").strip().upper()
sp3_file = "GNSS\igv23326_12.sp3"
position = get_data(sp3_file, target_satellite)


if not position[3]:
    print("未找到匹配数据")
else:
    target_time = datetime.strptime(input("输入目标时间 (YYYY-MM-DD HH:MM:SS): "),"%Y-%m-%d %H:%M:%S"
    )
    final_result = lagrange(position, target_time)
    print(f"卫星坐标（米）: X={final_result[0]:.3f}, Y={final_result[1]:.3f}, Z={final_result[2]:.3f}")
       
    #作图
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    #原始数据
    ax.plot(position[0], position[1], position[2], 'bo', label='original data')
    
    #插值数据
    ax.scatter(final_result[0], final_result[1], final_result[2], color='g', s=100, label='predicted point')
    
    ax.set_xlabel('X Coordinate')
    ax.set_ylabel('Y Coordinate')
    ax.set_zlabel('Z Coordinate')
    ax.set_title('3D Lagrange Interpolation')
    ax.legend()
    
    plt.show()

# --------------------------
obs_time = Time(position[3], scale="utc")  # 转换为Time对象
itrs_coord = ITRS(
    x=position[0] * u.m, 
    y=position[1] * u.m, 
    z=position[2] * u.m, 
    obstime=obs_time
)

# 转换为GCRS（天球坐标系）
gcrs_coord = itrs_coord.transform_to(GCRS(obstime=obs_time))  # 注意此处传递GCRS类

# 提取转换后的坐标
scale = 1e3  # 缩放因子
x_gcrs = gcrs_coord.cartesian.x.value 
y_gcrs = gcrs_coord.cartesian.y.value 
z_gcrs = gcrs_coord.cartesian.z.value 

def generate_earth(radius=6371):  # 地球平均半径6371公里
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 50)
    x = radius * np.outer(np.cos(u), np.sin(v))
    y = radius * np.outer(np.sin(u), np.sin(v))
    z = radius * np.outer(np.ones(np.size(u)), np.cos(v))
    return x, y, z

# --------------------------
# 绘图（3D空间坐标）
# --------------------------
fig = plt.figure(figsize=(12, 10))
ax = fig.add_subplot(111, projection='3d')

# 绘制地球表面
earth_x, earth_y, earth_z = generate_earth()
ax.plot_surface(earth_x, earth_y, earth_z, 
                color='lightblue', alpha=0.3, edgecolor='grey')

# 添加经纬度网格
theta = np.linspace(0, 2*np.pi, 100)
phi = np.linspace(0, np.pi, 50)
for ph in [0, np.pi/2]:  # 赤道与经线
    x = 6371e3 * np.cos(theta) * np.sin(ph)
    y = 6371e3 * np.sin(theta) * np.sin(ph)
    z = 6371e3 * np.cos(ph) * np.ones_like(theta)
    ax.plot(x/scale, y/scale, z/scale, color='gray', linestyle='--', alpha=0.5)

# 绘制坐标点
ax.scatter(x_gcrs, y_gcrs, z_gcrs, c='blue', s=100, label='GCRS (天球系)')

# 设置观察角度与比例
ax.view_init(elev=30, azim=45)  # 仰角30°, 方位角45°
ax.set_box_aspect([1,1,1])      # 强制等比例缩放
ax.set_xlabel('X ')
ax.set_ylabel('Y ')
ax.set_zlabel('Z ')
ax.set_title('地球模型与坐标系转换可视化\nITRF vs GCRS')
ax.legend()
# 自动调整坐标范围

plt.tight_layout()
plt.show()

