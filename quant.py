import pandas as pd
import numpy as np
import re

def process_file(file_path, file_name):
    # 读取文件
    df = pd.read_csv(file_path, sep='\t')
    
    # 1. 添加Pathway列（根据文件名映射）
    pathway_map = {
        'CBB': 'Calvin Benson Bassham (CBB) Cycle',
        '3-HP_4-HP_all': '3-Hydroxypropionate (3-HP) and 4-Hydroxypropionate (4-HP) Cycles'
    }
    
    base_name = file_name.split('.')[0]
    pathway_name = pathway_map.get(base_name, base_name)
    df['Pathway'] = pathway_name
    
    # 2. 从Gene列提取KO号
    def extract_ko(gene_str):
        match = re.search(r'\[(K\d{5})\]', gene_str)
        return match.group(1) if match else np.nan
    
    df['Gene_KO'] = df['Gene'].apply(extract_ko)
    
    # 3. 计算TPM值（按样本分组计算）
    tpm_values = []
    for sample, group in df.groupby('Sample'):
        total_fpkm = group['FPKM'].sum()
        tpm = (group['FPKM'] / total_fpkm) * 1e6
        tpm_values.extend(tpm)
    df['TPM'] = tpm_values
    
    # 4. 计算log10(TPM)，处理0值
    df['TPM(log10)'] = np.log10(df['TPM'].replace(0, 1e-10))
    
    # 5. 选择需要的列
    result_df = df[['Pathway', 'Sample', 'Gene_KO', 'TPM', 'TPM(log10)']]
    
    return result_df

# 处理两个文件
file1 = process_file('CBB.txt', 'CBB.txt')
file2 = process_file('3-HP_4-HP_all.txt', '3-HP_4-HP_all.txt')

# 合并结果
combined_df = pd.concat([file1, file2], ignore_index=True)

# 保存结果
combined_df.to_csv('combined_TPM_results.txt', sep='\t', index=False)

print("处理完成！结果已保存到 combined_TPM_results.txt")
