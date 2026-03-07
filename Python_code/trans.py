import pandas as pd
from Bio.Seq import Seq

# --- 配置 ---
INPUT_CSV = "final_result_with_symbols.csv"  # 刚才生成的文件
OUTPUT_ORDER = "Oligo_Order_Form.csv"


def get_reverse_complement(seq_str):
    return str(Seq(seq_str).reverse_complement())


def process_for_ordering(file_path):
    df = pd.read_csv(file_path)

    # 策略：每个 Gene_Symbol 只挑 1 个最好的 gRNA
    # 规则：优先选 Position 最小的，其次选 GC 最接近 55 的

    # 1. 转换 Position 为数值类型
    df['Position'] = pd.to_numeric(df['Position'])

    # 2. 排序：先按 Gene 分组，组内按 Position 从小到大排
    df_sorted = df.sort_values(by=['Gene_Symbol', 'Position'])

    # 3. 去重：每个基因只留第一行 (即 Position 最小的那个)
    df_unique = df_sorted.drop_duplicates(subset=['Gene_Symbol'], keep='first')

    order_list = []

    print(f"正在为 {len(df_unique)} 个基因生成合成序列...")

    for idx, row in df_unique.iterrows():
        gene = row['Gene_Symbol']
        seq = row['gRNA_Sequence']

        # --- 核心逻辑：加接头 (适配 BsmBI / lentiCRISPRv2) ---

        # 1. Top Strand (F)
        # 这里采用：如果首位不是G，则加G；如果是G，保持不变

        if seq[0].upper() == 'G':
            top_seq = "CACC" + seq
            # 对应的 Bottom 需要补配对的 sequence
            bottom_core = get_reverse_complement(seq)
            bottom_seq = "AAAC" + bottom_core
        else:
            top_seq = "CACCG" + seq  # 补一个G
            # 对应的 Bottom 核心要是 seq加了G 后的反向互补
            # 简单算法：Bottom = AAAC + RevComp(seq) + C
            bottom_core = get_reverse_complement(seq)
            bottom_seq = "AAAC" + bottom_core + "C"

        order_list.append({
            'Oligo_Name': f"{gene}_sgRNA_F",
            'Sequence': top_seq,
            'Description': 'Forward_Oligo'
        })

        order_list.append({
            'Oligo_Name': f"{gene}_sgRNA_R",
            'Sequence': bottom_seq,
            'Description': 'Reverse_Oligo'
        })

    return pd.DataFrame(order_list)


# --- 运行 ---
df_order = process_for_ordering(INPUT_CSV)
df_order.to_csv(OUTPUT_ORDER, index=False)
print(f"完成！已生成: {OUTPUT_ORDER}")