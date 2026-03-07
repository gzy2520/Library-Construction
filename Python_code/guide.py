import pandas as pd
import requests
import time
import os
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

# --- 配置部分 ---
INPUT_FILE = "nucleus.tsv"  # QuickGO输入文件
OUTPUT_FILE = "nucleus_gRNAs_results.csv"  # 结果保存文件
SPECIES = "human"


class CRISPRWorkflow:
    def __init__(self):
        self.server = "https://rest.ensembl.org"
        # 配置重试机制：自动重试5次，防止网络微小波动导致程序崩溃
        self.session = requests.Session()
        retries = Retry(total=5, backoff_factor=1, status_forcelist=[500, 502, 503, 504])
        self.session.mount('https://', HTTPAdapter(max_retries=retries))

    def load_symbols(self, filepath):
        """读取并去重Symbol"""
        try:
            df = pd.read_csv(filepath, sep='\t')
            df.columns = [c.strip() for c in df.columns]  # 清理空格

            if 'SYMBOL' not in df.columns:
                print("错误：未找到 SYMBOL 列")
                return []

            all_symbols = df['SYMBOL'].unique()
            # 过滤掉非字符串和线粒体基因
            nuclear_genes = [s for s in all_symbols if isinstance(s, str) and not s.startswith('MT-')]
            return sorted(list(set(nuclear_genes)))
        except Exception as e:
            print(f"读取文件出错: {e}")
            return []

    def get_gene_seq(self, symbol):
        """获取序列 (带自动重试)"""
        try:
            # 1. Symbol -> ID
            ext = f"/xrefs/symbol/{SPECIES}/{symbol}?"
            r = self.session.get(self.server + ext, headers={"Content-Type": "application/json"}, timeout=20)
            if not r.ok: return None

            genes = [x['id'] for x in r.json() if x['type'] == 'gene']
            if not genes: return None
            gene_id = genes[0]

            # 2. ID -> Canonical Transcript
            lookup_ext = f"/lookup/id/{gene_id}?expand=1"
            r_lookup = self.session.get(self.server + lookup_ext, headers={"Content-Type": "application/json"},
                                        timeout=20)
            if not r_lookup.ok: return None
            data = r_lookup.json()

            transcript_id = None
            if 'Transcript' in data:
                for t in data['Transcript']:
                    if t.get('is_canonical'):
                        transcript_id = t['id']
                        break
                if not transcript_id and data['Transcript']:
                    transcript_id = data['Transcript'][0]['id']

            if not transcript_id: return None

            # 3. Get Sequence
            seq_ext = f"/sequence/id/{transcript_id}?type=cdna"
            r_seq = self.session.get(self.server + seq_ext, headers={"Content-Type": "text/x-fasta"}, timeout=20)
            if not r_seq.ok: return None

            # 解析FASTA
            if '\n' in r_seq.text:
                sequence = r_seq.text.split('\n', 1)[1].replace('\n', '')
                return sequence
            return None

        except Exception as e:
            print(f" [Error: {e}]", end="")
            return None

    def design_guides(self, sequence, symbol):
        """设计 gRNA"""
        candidates = []
        scan_seq = sequence[:500]  # 只设计前500bp

        for i in range(len(scan_seq) - 23):
            sub = scan_seq[i: i + 23]
            spacer = sub[0:20]
            pam = sub[20:23]

            if pam[1:] == "GG":
                gc = (spacer.count('G') + spacer.count('C')) / 20 * 100
                if "TTTT" not in spacer and 40 <= gc <= 80:
                    candidates.append({
                        'Symbol': symbol,
                        'Sequence': spacer,
                        'PAM': pam,
                        'GC': round(gc, 1)
                    })

        candidates.sort(key=lambda x: abs(x['GC'] - 55))
        return candidates[:3]


def main():
    wf = CRISPRWorkflow()

    # 1. 加载所有基因
    all_genes = wf.load_symbols(INPUT_FILE)
    if not all_genes: return
    print(f"总任务量: {len(all_genes)} 个基因")

    # 2. 读取历史进度 (断点续传)
    processed_genes = set()
    if os.path.exists(OUTPUT_FILE):
        try:
            # 这里的trick是：即使是失败的行，只要有Symbol，也被视为“已处理”
            existing_df = pd.read_csv(OUTPUT_FILE)
            if 'Symbol' in existing_df.columns:
                processed_genes = set(existing_df['Symbol'].unique())
            print(f"检测到历史进度，已完成: {len(processed_genes)} 个 (含失败项)")
        except:
            print("历史文件读取失败，将新建。")

    # 如果是第一次运行，写入表头
    if not os.path.exists(OUTPUT_FILE):
        pd.DataFrame(columns=['Symbol', 'Sequence', 'PAM', 'GC']).to_csv(OUTPUT_FILE, index=False)

    # 3. 计算剩余任务
    genes_to_do = [g for g in all_genes if g not in processed_genes]
    print(f"本次需处理: {len(genes_to_do)} 个基因\n")

    if not genes_to_do:
        print("所有基因均已处理完毕！")
        return

    # 4. 开始跑
    for idx, gene in enumerate(genes_to_do):
        print(f"[{idx + 1}/{len(genes_to_do)}] 处理 {gene} ...", end="", flush=True)

        try:
            seq = wf.get_gene_seq(gene)
            saved = False

            if seq:
                guides = wf.design_guides(seq, gene)
                if guides:
                    # 成功：写入gRNA
                    df_new = pd.DataFrame(guides)
                    df_new.to_csv(OUTPUT_FILE, mode='a', header=False, index=False)
                    print(f" √ 存入 {len(guides)} 条")
                    saved = True

            if not saved:
                # 失败：写入一行 "N/A"，防止下次重复跑
                # 这一步非常关键！
                error_row = pd.DataFrame([{
                    'Symbol': gene,
                    'Sequence': 'N/A (Failed/No gRNA)',
                    'PAM': 'N/A',
                    'GC': 0
                }])
                error_row.to_csv(OUTPUT_FILE, mode='a', header=False, index=False)
                print(" × 无结果 (已记录)")

        except Exception as e:
            print(f" 未知错误: {e}")

        # 稍微礼貌一点，不要把API刷爆
        time.sleep(0.3)

    print("\n--- 全部任务完成 ---")


if __name__ == "__main__":
    main()