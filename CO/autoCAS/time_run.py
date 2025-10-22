#!/usr/bin/env python3
import time
import subprocess
import sys

def main():
    start_time = time.time()
    start_clock = time.perf_counter()
    start_process = time.process_time()
    
    print(f"开始时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 50)
    
    try:
        # 运行原始命令
        result = subprocess.run([
            sys.executable, "-m", "scine_autocas", "-y", "full.yml"
        ], check=True)
        
    except subprocess.CalledProcessError as e:
        print(f"程序执行失败: {e}")
        return 1
    
    end_time = time.time()
    end_clock = time.perf_counter()
    end_process = time.process_time()
    
    print("=" * 50)
    print(f"结束时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    
    # 计算各种时间
    wall_clock = end_time - start_time
    cpu_clock = end_clock - start_clock
    process_time = end_process - start_process
    
    print(f"\n时间统计:")
    print(f"墙钟时间: {wall_clock:.2f} 秒")
    print(f"CPU时间: {cpu_clock:.2f} 秒")
    print(f"进程时间: {process_time:.2f} 秒")
    
    # 格式化输出
    hours = int(wall_clock // 3600)
    minutes = int((wall_clock % 3600) // 60)
    seconds = wall_clock % 60
    
    print(f"总运行时间: {hours:02d}:{minutes:02d}:{seconds:06.3f}")

if __name__ == "__main__":
    main()
