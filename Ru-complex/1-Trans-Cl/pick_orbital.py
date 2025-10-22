def main(kesi):
    # 读取one-body.txt
    one_entries = []
    with open('one-body.txt', 'r') as f:
        for line in f:
            parts = line.strip().split(',')
            orb = int(parts[0])
            corr = float(parts[2])
            one_entries.append({'orbitals': [orb], 'corr': corr})
    
    # 读取two-body.txt
    two_entries = []
    with open('two-body.txt', 'r') as f:
        for line in f:
            parts = line.strip().split(',')
            orb_i = int(parts[0])
            orb_j = int(parts[1])
            corr = float(parts[3])
            two_entries.append({'orbitals': [orb_i, orb_j], 'corr': corr})
    
    # 合并所有条目并按相关能绝对值排序
    all_entries = one_entries + two_entries
    sorted_entries = sorted(all_entries, key=lambda x: abs(x['corr']), reverse=True)
    
    # 计算总绝对和
    total = sum(abs(entry['corr']) for entry in sorted_entries)
    if total == 0:
        print("总相关能绝对和为0，无法计算占比。")
        return []
    
    # 确定满足kesi的前n个条目
    current_sum = 0.0
    target = kesi * total
    selected_entries = []
    for entry in sorted_entries:
        if current_sum >= target:
            break
        selected_entries.append(entry)
        current_sum += abs(entry['corr'])
    
    # 收集轨道并去重
    orbitals = []
    for entry in selected_entries:
        orbitals.extend(entry['orbitals'])
    unique_orbitals = sorted(list(set(orbitals)))
    
    return unique_orbitals

if __name__ == '__main__':
    kesi = 0.06  # 可设定的参数
    result = main(kesi)
    print("选中的轨道：", result)
