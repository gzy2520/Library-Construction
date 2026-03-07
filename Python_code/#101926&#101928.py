import pandas as pd

# ==========================================
# 第二部分：处理 Addgene Guide 文件
# ==========================================
guides_df = pd.read_csv("guide.csv", header=None)


# 定义解析函数：提取第一部分 (ENSG ID) 和 第三部分 (SubLibrary)
def parse_guide_info(guide_str):
    try:
        parts = str(guide_str).split('_')
        if len(parts) >= 3:
            # 格式示例: ENSG00000158402_CDC25C_DTKP_19779.5
            gene_id = parts[0]  # 取第1段: Ensembl ID
            sub_lib = parts[2]  # 取第3段: 子库标签 (ACOC, GEEX等)
            return gene_id, sub_lib
    except:
        pass
    return None, None


# 应用解析
parsed_result = guides_df[0].apply(parse_guide_info).apply(pd.Series)
parsed_result.columns = ['ID', 'SubLib']

# 核心过滤：只保留以 'ENSG' 开头的有效 ID
valid_data = parsed_result[parsed_result['ID'].str.startswith('ENSG', na=False)]

# --- 生成 #101926 (ACOC) 列表 ---
ids_101926 = valid_data[valid_data['SubLib'] == 'ACOC']['ID'].unique()
pd.DataFrame(ids_101926, columns=['ID']).to_csv("For_Venn_101926_ID.csv", index=False)
print(f"-> [#101926 ACOC] 提取到 {len(ids_101926)} 个 ID")

# --- 生成 #101928 (GEEX) 列表 ---
ids_101928 = valid_data[valid_data['SubLib'] == 'GEEX']['ID'].unique()
pd.DataFrame(ids_101928, columns=['ID']).to_csv("For_Venn_101928_ID.csv", index=False)
print(f"-> [#101928 GEEX] 提取到 {len(ids_101928)} 个 ID")

# --- (可选) 生成 Bassik 全库 ID 列表 ---
# 如果你想看这两个子库加起来，或者整个大库的情况，可以用这个
ids_all = valid_data['ID'].unique()
pd.DataFrame(ids_all, columns=['ID']).to_csv("For_Venn_All_Bassik_ID.csv", index=False)
print(f"-> [Bassik 全库] 提取到 {len(ids_all)} 个 ID")

print("\n所有文件已生成完毕！")