import math

# 角度与弧度转换
def deg2rad(deg):
    return deg * math.pi / 180.0

def rad2deg(rad):
    return rad * 180.0 / math.pi

# BLH → XYZ
def BLH_to_XYZ(B_deg, L_deg, H, a, b):
    B = deg2rad(B_deg)
    L = deg2rad(L_deg)
    e2 = (a**2 - b**2) / a**2  # 第一偏心率平方
    W = math.sqrt(1 - e2 * math.sin(B)**2)
    N = a / W

    X = (N + H) * math.cos(B) * math.cos(L)
    Y = (N + H) * math.cos(B) * math.sin(L)
    Z = (N * (1 - e2) + H) * math.sin(B)
    return X, Y, Z

# XYZ → BLH
def XYZ_to_BLH(X, Y, Z, a, b):
    e2 = (a**2 - b**2) / a**2
    ep2 = (a**2 - b**2) / b**2

    L = math.atan2(Y, X)
    p = math.sqrt(X**2 + Y**2)
    tgBi = Z / p
    c = a * math.sqrt(1 + ep2)
    eps = 1e-12

    while True:
        tgBi1 = Z / p + (a * e2 * tgBi) / (p * math.sqrt(1 + tgBi**2) * math.sqrt(1 - e2 * tgBi**2 / (1 + tgBi**2)))
        Bi1 = math.atan(tgBi1)
        Bi = math.atan(tgBi)
        if abs(Bi1 - Bi) < eps:
            break
        tgBi = tgBi1

    B = Bi1
    N = a / math.sqrt(1 - e2 * math.sin(B)**2)
    H = p / math.cos(B) - N
    return rad2deg(B), rad2deg(L), H

# 主程序
def main():
    ellipsoids = {
        'WGS84': (6378137.0, 6356752.314245),
        'Krassovsky': (6378245.0, 6356863.0188)
    }

    choice = input("选择椭球 (WGS84/Krassovsky)：").strip()
    if choice not in ellipsoids:
        print("未知椭球，默认使用 WGS84")
        choice = 'WGS84'

    a, b = ellipsoids[choice]

    # 原始大地坐标
    B_deg = 44 + 59.0 / 60 + 59.9999 / 3600
    L_deg = 45.0
    H = 999999.9987

    # BLH → XYZ
    X, Y, Z = BLH_to_XYZ(B_deg, L_deg, H, a, b)
    print("\nXYZ 坐标:")
    print(f"X = {X:.12f}")
    print(f"Y = {Y:.12f}")
    print(f"Z = {Z:.12f}")

    # XYZ → BLH
    B1, L1, H1 = XYZ_to_BLH(X, Y, Z, a, b)
    print("\n恢复后的 BLH 坐标:")
    print(f"B = {B1:.12f}")
    print(f"L = {L1:.12f}")
    print(f"H = {H1:.12f}")

    # 差异检查
    B_diff_arcsec = abs(B1 - B_deg) * 3600
    L_diff_arcsec = abs(L1 - L_deg) * 3600
    H_diff_cm = abs(H1 - H) * 100

    print("\n误差检查:")
    print(f"B 误差（角秒）: {B_diff_arcsec:.6f}")
    print(f"L 误差（角秒）: {L_diff_arcsec:.6f}")
    print(f"H 误差（厘米）: {H_diff_cm:.6f}")

    if B_diff_arcsec < 0.1 and L_diff_arcsec < 0.1 and H_diff_cm < 1.0:
        print("转换精度满足要求。")
    else:
        print("转换精度不满足要求。")

if __name__ == "__main__":
    main()
