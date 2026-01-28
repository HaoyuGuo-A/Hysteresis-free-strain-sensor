import serial
import time
import threading
from collections import deque
import os
import csv
from datetime import datetime

# 寄存器地址
register_dict = {
    'ID': 1000,
    'clearErr': 1004,
    'angleSet': 1486,
    
}

# 每个手指的映射区间和对应角度 小拇指 无名指  中指 食指 拇指
#[15, 25, 35, 45, 55, 65, 75, 85, 95, 105]   [1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000]
finger_mappings = {
     'pinky': {
        'intervals': [87.13, 87.63, 88.14, 88.64, 89.15, 89.66, 90.16, 90.67, 91.17, 91.68, 92.19, 92.69, 93.2, 93.7, 94.21, 94.72, 95.22, 95.73, 96.23, 96.74, 97.25],
        'angles': [1000, 950, 900, 850, 800, 750, 700, 650, 600, 550, 500, 450, 400, 350, 300, 250, 200, 150, 100, 50, 0]
    },
    'ring': {
        'intervals': [74.13, 74.9, 75.68, 76.46, 77.24, 78.02, 78.8, 79.58, 80.36, 81.14, 81.92, 82.7, 83.48, 84.26, 85.04, 85.82, 86.6, 87.38, 88.16, 88.94, 89.72],
        'angles': [1000, 950, 900, 850, 800, 750, 700, 650, 600, 550, 500, 450, 400, 350, 300, 250, 200, 150, 100, 50, 0]
    },
    'middle': {
        'intervals': [90.75, 91.93, 93.1, 94.28, 95.45, 96.63, 97.81, 98.98, 100.16, 101.34, 102.51, 103.69, 104.86, 106.04, 107.22, 108.39, 109.57, 110.74, 111.92, 113.1, 114.27],
        'angles':  [1000, 950, 900, 850, 800, 750, 700, 650, 600, 550, 500, 450, 400, 350, 300, 250, 200, 150, 100, 50, 0]
    },
    'index': {
        'intervals': [81.6, 82.52, 83.45, 84.37, 85.3, 86.22, 87.14, 88.07, 88.99, 89.92, 90.84, 91.76, 92.69, 93.61, 94.54, 95.46, 96.38, 97.31, 98.23, 99.16, 100.08],
        'angles': [1000, 950, 900, 850, 800, 750, 700, 650, 600, 550, 500, 450, 400, 350, 300, 250, 200, 150, 100, 50, 0]
    },
    'thumb': {
        'intervals': [77, 77.2, 77.39, 77.59, 77.78, 77.98, 78.17, 78.37, 78.56, 78.76, 78.95, 79.15, 79.34, 79.54, 79.73, 79.93, 80.12, 80.32, 80.52, 80.71, 80.91],
        'angles': [1000, 950, 900, 850, 800, 750, 700, 650, 600, 550, 500, 450, 400, 350, 300, 250, 200, 150, 100, 50, 0]
    }
}


# ---- 电容值到角度的映射函数 ---- #
def get_angle_from_cap(cap, finger):
    mapping = finger_mappings.get(finger)
    if not mapping:
        print(f"未找到手指 {finger} 的映射配置。")
        return 1000
    intervals = mapping['intervals']
    angles = mapping['angles']
    for i, threshold in enumerate(intervals):
        if cap < threshold:
            return angles[i]
    return angles[-1]


# ---- 打开串口 ---- #
def openSerial(port, baud_rate):
    try:
        ser_conn = serial.Serial(port=port, baudrate=baud_rate, timeout=1)
        print(f"串口 {port} 已打开，波特率: {baud_rate}")
        return ser_conn
    except serial.SerialException as e:
        print(f"无法打开串口 {port}: {e}")
        exit(1)


# ---- 写寄存器 ---- #
def writeRegister(ser_conn, device_id, add, num, val):
    data_bytes = [0xEB, 0x90, device_id, num + 3, 0x12, add & 0xFF, (add >> 8) & 0xFF]
    data_bytes.extend(val)
    checksum = sum(data_bytes[2:]) & 0xFF
    data_bytes.append(checksum)
    try:
        ser_conn.write(bytes(data_bytes))
        time.sleep(0.01)
        ser_conn.read_all()
    except serial.SerialException as e:
        print(f"写寄存器失败: {e}")


# ---- 写6个寄存器 ---- #
def write6(ser_conn, device_id, reg, values):
    values = values.copy()
    values.append(1000)
    val_reg = [val & 0xFF for value in values for val in (value, value >> 8)]
    writeRegister(ser_conn, device_id, register_dict[reg], 12, val_reg)


# ---- 将字节数据转换为浮点数 ---- #
def bytes_to_float(sensor_data):
    integer_part = int.from_bytes(sensor_data[:3], byteorder='big', signed=True)
    fractional_part = int.from_bytes(sensor_data[3:], byteorder='big', signed=False)
    return round(integer_part + fractional_part / 1000000.0, 6)


# ---- 全局同步变量 ---- #
current_angles = [None, None, None, None, None]
capacitance_values = [0.0, 0.0, 0.0, 0.0, 0.0]
target_angles = [1000, 1000, 1000, 1000, 1000]
old_angles_6 = [0, 0, 0, 0, 0, 0]
WINDOW_SIZE = 5
capacitance_buffers = [deque(maxlen=WINDOW_SIZE) for _ in range(5)]
filtered_capacitance = [0.0, 0.0, 0.0, 0.0, 0.0]


# ---- 计算滑动平均的函数 ---- #
def calculate_moving_average(buffer):
    return sum(buffer) / len(buffer) if buffer else 0.0


# ---- 创建保存数据的文件夹和文件 ---- #
def setup_data_saving():
    # 获取桌面路径
    desktop_path = os.path.join(os.path.expanduser("~"), "Desktop")

    # 创建以年月日命名的文件夹
    current_date = datetime.now().strftime("%Y%m%d")
    folder_path = os.path.join(desktop_path, current_date)
    os.makedirs(folder_path, exist_ok=True)

    # 创建以时分秒命名的CSV文件
    current_time = datetime.now().strftime("%H%M%S")
    csv_filename = f"{current_time}.csv"
    csv_path = os.path.join(folder_path, csv_filename)

    # 初始化CSV文件并写入表头
    with open(csv_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Timestamp', 'Pinky', 'Ring', 'Middle', 'Index', 'Thumb'])

    return csv_path


# ---- 保存电容值到CSV ---- #
def save_capacitance_to_csv(csv_path, cap_values):
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S.%f")
    with open(csv_path, 'a', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow([timestamp] + cap_values)


# ---- 机械手操作线程 ---- #
def mechanical_hand_thread():
    global current_angles, target_angles, old_angles_6
    ser_conn = openSerial('COM3', 115200) ## ----------------修改成机械手的串口号
    writeRegister(ser_conn, 1, register_dict['clearErr'], 1, [1])
    print("机械手控制线程已启动。")
    while True:
        angles_6 = target_angles + [1000]
        if angles_6 != old_angles_6:
            write6(ser_conn, 1, 'angleSet', target_angles)
            old_angles_6 = angles_6[:]
            current_angles = target_angles[:]
            print("机械手角度已更新：", angles_6)
        time.sleep(0.05)


# ---- 数据采集与映射线程 ---- #
def data_acquisition_thread(csv_path):
    global capacitance_values, target_angles, current_angles
    ser = openSerial('COM10', 230400)       # ----------------修改成采集器的串口号
    initial_bytes = bytes([
        0xFF, 0xFF, 0x06, 0x09, 0x00, 0x00, 0x00, 0x0C,
        0xC0, 0xA8, 0x13, 0x93, 0xFF, 0xFF, 0xFF, 0xFF,
        0x00, 0x11, 0x00, 0x01
    ])
    ser.write(initial_bytes)
    time.sleep(3)
    buffer = b""
    print("数据采集线程已启动。")
    fingers = ['pinky', 'ring', 'middle', 'index', 'thumb']

    while True:
        try:
            data = ser.read(1024)
            if not data:
                continue
            buffer += data
            while len(buffer) >= 50:
                if buffer[:4] == b'\xFF\xFF\x06\x09':
                    sensor_data = buffer[20:50]
                    raw_cap_vals = [
                        bytes_to_float(sensor_data[i * 6: (i + 1) * 6])
                        for i in range(5)
                    ]
                    reordered = [
                        raw_cap_vals[4],  # 小拇指
                        raw_cap_vals[3],  # 无名指
                        raw_cap_vals[2],  # 中指
                        raw_cap_vals[1],  # 食指
                        raw_cap_vals[0],  # 大拇指
                    ]
                    capacitance_values = reordered

                    # 保存电容值到CSV
                    save_capacitance_to_csv(csv_path, capacitance_values)

                    for i, cap in enumerate(capacitance_values):
                        capacitance_buffers[i].append(cap)
                        filtered_capacitance[i] = calculate_moving_average(capacitance_buffers[i])

                    new_target = []
                    for cap, finger in zip(filtered_capacitance, fingers):
                        finger_angle = get_angle_from_cap(cap, finger)
                        new_target.append(finger_angle)

                    target_angles = new_target
                    print(f"接收到电容值: {capacitance_values}, 映射后的目标角度: {target_angles}")
                    buffer = buffer[50:]
                else:
                    buffer = buffer[1:]
        except serial.SerialException as e:
            print(f"数据采集失败: {e}")
            break


# ---- 主函数 ---- #
def main():
    # 设置数据保存路径
    csv_path = setup_data_saving()

 #加一个设置机械手速度的
   # write6(ser_conn, 1, register_dict['speedSet'], [1000, 1000, 1000, 1000, 1000, 1000])
    #write6(ser_conn, 1, register_dict['angleSet'], [1000, 1000, 1000, 1000, 1000, 0])
    #write6(ser_conn, 1, register_dict['forceSet'], [1000, 1000, 1000, 1000, 1000, 1000])

    # 启动两个线程
    t1 = threading.Thread(target=mechanical_hand_thread, daemon=True)
    t2 = threading.Thread(target=data_acquisition_thread, args=(csv_path,), daemon=True)
    t1.start()
    t2.start()
    print("主程序已启动，按 Ctrl+C 退出。")
    try:
        while True:
            time.sleep(1)
    except KeyboardInterrupt:
        print("程序已中断，退出。")
        exit(0)


if __name__ == '__main__':
    main()