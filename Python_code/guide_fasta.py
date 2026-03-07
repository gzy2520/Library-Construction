import pandas as pd
import re
import mygene
import sys

# --- 配置 ---
FASTA_FILE = "nucleus_sequences.fasta"
OUTPUT_FILE = "final_result_with_symbols.csv"


def parse_fasta_and_extract_ids(file_path):
    """
    读取FASTA，提取序列，并从表头中抓取 RefSeq ID (NM_xxxxx)
    """
    entries = []
    current_id = None
    current_seq = []

    # 正则表达式：匹配 NM_ 开头的编号
    # 表头示例: >hg38_ncbiRefSeqCurated_NM_001005337.3 range=...
    id_pattern = re.compile(r'(NM_\d+)')

    print("1. 正在读取 FASTA 文件...")
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # 保存上一条
                if current_id and current_seq:
                    entries.append({'id': current_id, 'seq': "".join(current_seq)})

                # 解析新表头
                match = id_pattern.search(line)
                if match:
                    current_id = match.group(1)  # 提取 NM_001005337
                else:
                    current_id = "Unknown"

                current_seq = []
            else:
                current_seq.append(line)

        # 保存最后一条
        if current_id and current_seq:
            entries.append({'id': current_id, 'seq': "".join(current_seq)})

    return pd.DataFrame(entries)


def map_ids_to_symbols(df):
    """
    使用 MyGene 将 NM_ ID 批量转换为 Gene Symbol
    """
    print(f"2. 正在将 {len(df)} 个 RefSeq ID 转换为基因名...")

    mg = mygene.MyGeneInfo()
    id_list = df['id'].unique().tolist()

    # 批量查询
    results = mg.querymany(id_list, scopes='refseq', fields='symbol', species='human')

    # 建立映射字典
    id_map = {}
    for res in results:
        query_id = res['query']
        if 'symbol' in res:
            id_map[query_id] = res['symbol']
        else:
            id_map[query_id] = query_id  # 如果查不到，就保留原ID

    # 将 Symbol 映射回 DataFrame
    df['Symbol'] = df['id'].map(id_map)
    return df


def design_guides(sequence, symbol, refseq_id):
    """gRNA 设计逻辑"""
    candidates = []
    # 已经是 CDS (外显子拼接)，所以可以直接设计
    # 只取前 500bp 提高 KO 效率
    scan_seq = sequence[:500]

    for i in range(len(scan_seq) - 23):
        sub = scan_seq[i: i + 23]
        spacer = sub[0:20]
        pam = sub[20:23]

        if pam[1:] == "GG":
            gc = (spacer.count('G') + spacer.count('C')) / 20 * 100

            if "TTTT" not in spacer and 40 <= gc <= 80:
                candidates.append({
                    'Gene_Symbol': symbol,
                    'RefSeq_ID': refseq_id,
                    'gRNA_Sequence': spacer,
                    'PAM': pam,
                    'GC_Content': round(gc, 1),
                    'Position': i + 1
                })

    # 按 GC 含量排序 (最接近55%的排前面)
    candidates.sort(key=lambda x: abs(x['GC_Content'] - 55))
    return candidates[:3]  # 每个转录本只保留前3个


def main():
    # 1. 解析文件
    try:
        df_seq = parse_fasta_and_extract_ids(FASTA_FILE)
    except FileNotFoundError:
        print(f"找不到文件: {FASTA_FILE}")
        return

    if df_seq.empty:
        print("文件中没有提取到序列，请检查内容。")
        return

    # 2. 转换 ID -> Symbol
    df_mapped = map_ids_to_symbols(df_seq)

    # 3. 设计 gRNA
    print("3. 开始批量设计 gRNA...")
    all_guides = []

    for index, row in df_mapped.iterrows():
        symbol = row['Symbol']
        seq = row['seq']
        refseq_id = row['id']

        # 简单进度条
        if index % 1000 == 0:
            print(f"   已处理 {index} 条序列...")

        guides = design_guides(seq, symbol, refseq_id)
        all_guides.extend(guides)

    # 4. 保存结果
    if all_guides:
        result_df = pd.DataFrame(all_guides)
        # 顺序
        cols = ['Gene_Symbol', 'RefSeq_ID', 'gRNA_Sequence', 'PAM', 'GC_Content', 'Position']
        result_df = result_df[cols]

        result_df.to_csv(OUTPUT_FILE, index=False)
        print(f"\n成功！结果已保存至: {OUTPUT_FILE}")
        print(f"共生成 {len(result_df)} 条 gRNA。")
        print("前5行预览：")
        print(result_df.head())
    else:
        print("未生成任何 gRNA。")


if __name__ == "__main__":
    main()