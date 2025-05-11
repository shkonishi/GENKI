# GENKI
原核生物の基準株ゲノムデータベースを構築する。ゲノム情報を基準として、以下のように異なるデータベースから取得したデータを紐付けて管理する。
- **ゲノムデータ**: NCBI GenBank
- **系統情報**: NCBI Taxonomy
- **表現型データ**: BacDive

## **依存**
- bash
- taxonkit
- BacDive (R-package of API Client)
- (ncbi-genome-download)

## **ゲノムデータダウンロード**
### **assembly summaryの形式**
アセンブリサマリに記述されている情報について下記リンクのファイルから把握しておく(以下2025/05/08時点での情報)
https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/README_assembly_summary.txt

- `Column 25: group`: タクソノミーグループの情報(bacteria, archaea, その他)が記述されている。
```bash
cut -f25 assembly_summary_genbank_20250508T2015.txt | sort | uniq -c | sort -nr -k1,1
# 2657786 bacteria
#  229681 viral
#   30508 archaea
#   21468 fungi
#   17875 metagenomes
#   10393 invertebrate
#    7596 other
#    6646 plant
#    6061 vertebrate_other
#    4393 vertebrate_mammalian
#    2408 protozoa

```

- `Column 22: relation_to_type_material`列にタイプとの関連情報が記述されている
```bash
cut -f22 assembly_summary_genbank_20250508T2015.txt | sort | uniq -c | sort -nr -k1,1
# 2950577 na
#   29694 assembly from type material
#   12051 ICTV species exemplar
#    1660 ICTV additional isolate
#     607 assembly from synonym type material
#      93 assembly from pathotype material
#      90 assembly designated as neotype
#      25 assembly designated as reftype
#      18 assembly designated as clade exemplar
# type
```

- `Column 5: refseq_category` RefSeqプロジェクトの分類において、アセンブリが参照ゲノムであるかどうか
- `Column 12: assembly_level` Complete genome|Chromosome|Scaffold|Contig
- ``

      

### genki_download.shを使う場合(現時点では関数を個別に実行する仕様)  
#### **ダウンロード実行**
- NCBI FTPサーバーリンクからサマリーデータをダウンロードする

```bash
# スクリプト仕様(予定)
# genki_download.sh --dry-run -g bacteria,archae -r type,synonym > assembly_summary_select.txt
# genki_download.sh -g bacteria,archae -r type,synonym -o out_dir -i assembly_summary_select.txt
# genki_download.sh -g bacteria,archae -r type,synonym -o out_dir -i GCA_XXX,GCA_YYY,GCA_ZZZ
# genki_update.sh -i assembly_summary_new.txt -I assembly_summary_new.txt

# ユーザー定義関数
## 1. assembly_summaryをDL(genbankまたはrefseqを選択)
## get_assembly_summary 'gbk'
## 2. assembly_summaryを指定条件でフィルタリング
## filter_summary assembly_summary > filtered_assembly_summary
## geturl_assembly filtered_assembly_summary ./out_dir

```

```bash
# ユーザー定義関数をロード
source ../lib/get_assembly_processing.sh
# アセンブリサマリのダウンロード
VER=$(date +"%Y%m%dT%H%M")

## GenBank
wget -c -O "assembly_summary_genbank_${VER}.txt" https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt

# 指定条件のゲノムデータをダウンロードする
IN_SUM="assembly_summary_genbank_${VER}.txt"
OUT_DIR="out_genki_${VER}"
OUT_DAT_DIR="${OUT_DIR}/dat"
OUT_LOG_DIR="${OUT_DIR}/log"
OUT_ASMB_DIR="${OUT_DIR}/genomes"

out_sum="${OUT_DAT_DIR}/assembly_summary_genbank_type_prok_${VER}.txt"
out_md5="${OUT_ASMB_DIR}/md5_genomes"
out_md5chk="${OUT_DIR}/res_md5sum.txt"
out_md5fail="${OUT_DIR}/res_md5sum_err"

# サマリデータフィルタリング
mkdir -p "${OUT_DAT_DIR}"
filter_summary "$IN_SUM" "$out_sum"

# URLファイル作成
geturl_assembly "$out_sum" "$OUT_LOG_DIR"

# データダウンロード
get_gca_genomic "${OUT_LOG_DIR}/ftp_genome" "$OUT_ASMB_DIR" &
wait
get_md5_genomic "${OUT_LOG_DIR}/ftp_md5" "$OUT_ASMB_DIR" &
wait
# md5sumを実行
(cd "$OUT_ASMB_DIR" && md5sum -c ./md5_genomes) > "$out_md5chk" 2> "$out_md5fail"

# Log集計
# DL実行した数
wc -l ./out_genki_20250508T2015/log/ftp_genome
# DL失敗したものを確認
cat out_genki_20250508T2015/wget_gca.err
# MD5失敗したものを確認
awk -F":" 'BEGIN{ok=0; fail=0; }$2~/OK/{ok++;}$2!~/OK/{fail++;}END{print "OK:" ok "\t" "FAIL:" fail}' < "out_genki_20250508T2015/res_md5sum.txt"


```

- 追加登録のみをダウンロードする場合
- 追加登録の情報情報を抽出したassembly_addition.txtが作成される。 

```bash
VER=$(date +"%Y%m%dT%H%M")
IN_SUM="assembly_summary_genbank_${VER}.txt"
OUT_SUM="assembly_summary_genbank_type_prok_${VER}.txt"
filter_summary "$IN_SUM" "$OUT_SUM"

IN_SUM_OLD="out_genki_20250508T2015/dat/assembly_summary_genbank_type_prok_20250508T2015.txt"
IN_SUM_NEW="./assembly_summary_genbank_type_prok_20250511T1407.txt"
# genki_update.sh "$IN_SUM_OLD" "$IN_SUM_NEW"

```

## **データ構造**
```sh
out_genki_241229/
├── dat
│   └── assembly_summary_genbank_type_prok_241229.txt # フィルタしたアセンブリサマリ
├── genomes
│   ├── GCA_000003925.1_ASM392v1_genomic.fna.gz
│   ├── GCA_000006685.1_ASM668v1_genomic.fna.gz
│   ├── GCA_000006945.2_ASM694v2_genomic.fna.gz
│   ├── GCA_000006985.1_ASM698v1_genomic.fna.gz
│   ├── GCA_000007025.1_ASM702v1_genomic.fna.gz
│   ├── GCA_000007185.1_ASM718v1_genomic.fna.gz
│   └── md5_genomes
├── log
│   ├── ftp_genome # ゲノムデータURLリスト
│   └── ftp_md5    # md5ファイルURLリスト
├── res_md5sum.txt  # md5sumの結果
├── res_md5sum_err  # md5sumを失敗したファイルリスト
├── wget_gca.err    # wgetのエラーログ
└── wget_md5.err    # wgetのエラーログ


```

## **系統情報の取得**
系統情報についてはアセンブリサマリーからも確認できるが、記載内容が一貫していなかったり不足している場合がある。ここではtaxonomy-idをもとに系統情報を取得する。

### **最新版のtaxdumpを入手する**
```bash
wget -c https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
mv taxdump.tar.gz ~/.taxonkit
cd ~/.taxonkit
tar xvf taxdump.tar.gz
```
### **Taxonomy IDを使用して系統情報を取得**
- taxonkitを使ってlineageを取得&整形 (ver. 0.20.2)
- idとtaxonomyの変換表を作成しておく。
    ゲノムデータを参照データとした解析に利用する。

```bash
IN_TYPE='./dat/assembly_summary_genbank_type_prok.txt'
OUT_TAX='./dat/assembly_summary_genbank_type_prok_lineage.txt'
OUT_TAB='./dat/id_lookup.tsv'

# cut -f1,6 $IN_TYPE | taxonkit lineage -i 2 | taxonkit reformat -i 3 -f "{k};{p};{c};{o};{f};{g};{s};{t}" | cut -f1,2,4 > "$OUT_TAX"
cut -f1,6 "$IN_TYPE" | taxonkit reformat2 -I 2 -f "{domain};{phylum};{class};{order};{family};{genus};{species};{strain|subspecies|no rank}" > "$OUT_TAX"

# idルックアップデーブルの作成
# cat  | awk -F"\t" '{split($3,arr,";"); gsub(" ","_", arr[7]); print $1"\t"arr[7]"_"$1;}' > "$OUT_TAB"
cut -f1 "$OUT_TAX" | while read -r line; do
    awk -F"\t" -v gca="$line" '$1==gca{ split($8,arr," "); sub("type strain: ","",$9); sub("strain=","",$9); print $1"\t"arr[1] "_" arr[2] "_"$9 "_" $1}' "$IN_TYPE"
done > "$OUT_TAB"

```
- 以下のようにdeleteやnot foundの警告がでなければO.K.
```txt
22:28:36.509 [WARN] taxid 3238480 not found
22:28:36.509 [WARN] taxid 3133966 was deleted
```


## **表現型データの取得**
- 予めBacDiveのサイト(以下URL)からタイプストレインのリストcsvをダウンロードしておく
    https://bacdive.dsmz.de/dashboard
    1. BacDiveの公式サイトへアクセス
    2. 「Advanced search」をクリック  
    3. Choose fieldの欄のFilter条件に以下の項目を選択してSubmit
        + Type strainでyesを選択
        + Genome Sequence databaseでncbiを選択 
    4. 「Download table as CSV」をクリック

```bash
in_bdcsv='advsearch_bacdive_2024-08-02.csv'
Rscript this_script.r "$in_bdcsv" <API_ID> <API_PASSWORD>
```

## **skaniのsketchファイルを作成**
- ANI検索を実行するためのdbを作成しておく
- 工事中

