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
## get_assembly_summary <gbk|ref> <output>
## 2.1. assembly_summaryを指定条件でフィルタリング
## filter_summary -i <assembly_summary> -g <group> -y <relation_to_type> > assembly_summary_filtered.tsv
## 2.2. update されたassembly_summaryを取得(削除された)
#  update_summary -o <assembly_summary_old> -n <assembly_summary_new> [-r <removed_id>] [-a <added_id>]
## 3. genomeおよびmd5のftp-urlを取得
## geturl_assembly -i <input> [--out_genome_url <outfile>] [--out_md5_url <outfile>]
## 4. genomeおよびmd5のダウンロード
## get_gca_genomic <genomes_ftp_list> [<output_dir>] 
## get_md5_genomic <md5_ftp_list> [<output_dir>]

```

```bash
# ユーザー定義関数をロード
source ../lib/get_assembly_processing.sh

# 作業ディレクトリ構築 & Logファイル定義
VER=$(date +"%Y%m%dT%H%M")
OUT_DIR="GENKI_${VER}"
[[ ! -d "$OUT_DIR" ]] && mkdir -p "${OUT_DIR}"
LOG=${VER}_genki.log

# assembly summaryのダウンロード
OUTSUMDIR="${OUT_DIR}"/.genki
[[ ! -d "$OUTSUMDIR" ]] && mkdir -p "${OUTSUMDIR}"
out_sum="${OUTSUMDIR}"/assembly_summary_${VER}.txt

get_assembly_summary gbk "$out_sum"

# assembly summaryのフィルタリング
OUT_DAT_DIR="${OUT_DIR}/dat"
[[ ! -d "$OUT_DAT_DIR" ]] && mkdir -p "${OUT_DAT_DIR}"
out_filtered_sum="${OUT_DAT_DIR}/assembly_summary_filtered_${VER}.txt"

filter_gbsummary -i "$out_sum" > "$out_filtered_sum"

# assembly データのダウンロード
OUT_DAT_DIR="${OUT_DIR}/dat"
OUT_LOG_DIR="${OUT_DIR}/log"
OUT_ASMB_DIR="${OUT_DIR}/genomes"

out_md5="${OUT_ASMB_DIR}/md5_genomes"
out_md5chk="${OUT_LOG_DIR}/res_md5sum.txt"
out_md5fail="${OUT_LOG_DIR}/res_md5sum_err"

## URLファイル作成
geturl_assembly "$out_filtered_sum" "$OUT_LOG_DIR"

## データダウンロード
get_gca_genomic "${OUT_LOG_DIR}/ftp_genome" "$OUT_ASMB_DIR" &
wait
get_md5_genomic "${OUT_LOG_DIR}/ftp_md5" "$OUT_ASMB_DIR" &
wait
## md5sumを実行
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

```bash
# 作業ディレクトリ構築 & Logファイル定義
VER=$(date +"%Y%m%dT%H%M")
OUT_DIR=~/db/GENKI_250506 # 既に作成済みのディレクトリ
OUT_UPDIR="${OUT_DIR}/update_${VER}"
[[ ! -d "$OUT_UPDIR" ]] && mkdir -p "${OUT_UPDIR}"
LOG=${OUT_UPDIR}/${VER}_genki.log

# assembly summaryのダウンロード => フィルタリング
OUTSUMDIR="${OUT_DIR}"/.genki
[[ ! -d "$OUTSUMDIR" ]] && mkdir -p "${OUTSUMDIR}"
out_sum="${OUTSUMDIR}"/assembly_summary_${VER}.txt
get_assembly_summary gbk "$out_sum"
NEW_SUM="${OUT_DIR}/dat/assembly_summary_filtered_${VER}.txt"
filter_gbsummary -i "$out_sum" -g bacteria,archaea -y type,synonym > "$NEW_SUM" 2>> "$LOG"

# update_idを取得 =>  update urlを取得
OLD_SUM="${OUT_DIR}/dat/assembly_summary_filtered_20250508T2015.txt" # 前回作成したフィルタ済みassembly_summary
update_summary -o $OLD_SUM -n $NEW_SUM -r "${OUT_UPDIR}"/assembly_removed.tsv -a "${OUT_UPDIR}"/assembly_addition.tsv 2>> "$LOG"
geturl_assembly -i "${OUT_UPDIR}"/assembly_addition.tsv -g "${OUT_UPDIR}/ftp_genomes.url" -r "${OUT_UPDIR}/ftp_md5.url" 2>> "$LOG"


# ダウンロード
OUT_ASMB_DIR="${OUT_UPDIR}/genomes"
get_gca_genomic "${OUT_UPDIR}/ftp_genomes.url" "$OUT_ASMB_DIR" &
wait
get_md5_genomic "${OUT_UPDIR}/ftp_md5.url" "$OUT_ASMB_DIR" &
wait

## md5sumを実行
out_md5chk="${OUT_UPDIR}/res_md5sum.txt"
out_md5fail="${OUT_UPDIR}/res_md5sum.err"
(cd "$OUT_ASMB_DIR" && md5sum -c ./md5_genomes) > "$out_md5chk" 2> "$out_md5fail"

# Summary of update
# DL実行した数
num_dl=$(cat "${OUT_UPDIR}/ftp_genomes.url" | wc -l) 
# DL失敗したものを確認
cat "${OUT_UPDIR}/wget_gca.err"
# MD5失敗したものを確認
res_md5chk=$(awk -F":" 'BEGIN{ok=0; fail=0; }$2~/OK/{ok++;}$2!~/OK/{fail++;}END{print "OK:" ok "\t" "FAIL:" fail}' < "${OUT_UPDIR}/res_md5sum.txt")

echo "[INFO] Number of downloads executed: $num_dl" >&2
echo "[INFO] Check for download failures: " >&2 ; cat "${OUT_UPDIR}/wget_gca.err" >&2
echo "[INFO] Check MD5 checksum value mismatch: ${res_md5chk}" >&2 ; cat "${OUT_UPDIR}/wget_md5.err" >&2



```

## **データ構造**
```sh
~/db/GENKI_250506/
├── blast_db        # blast データベース
├── dat             # フィルタ済みアセンブルデータ
├── genomes         # ゲノムデータ
│   ├── GCA_000003925.1_ASM392v1_genomic.fna.gz
│   ├── GCA_000006685.1_ASM668v1_genomic.fna.gz
│   ├── GCA_000006945.2_ASM694v2_genomic.fna.gz
│   ├── :
│   └── md5_genomes
├── log
│   ├── ftp_genome         # ゲノムデータURLリスト
│   ├── ftp_md5            # md5ファイルURLリスト
│   ├── res_md5sum.txt     # md5sumの結果
│   └── res_md5sum.err     # md5sumを失敗したファイルリスト
├── skani_db 
├── update_20250518T1838    #
├── wget_gca.err            # wgetのエラーログ
└── wget_md5.err            # wgetのエラーログ


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
IN_TYPE="./dat/assembly_summary_filtered_20250518T0224.txt"
OUT_TAX='./dat/lineage.tsv'
OUT_TAB='./dat/id_lookup.tsv'

# Lineageの取得
cut -f1,6 "$IN_TYPE" | taxonkit reformat2 -I 2 -f "{domain};{phylum};{class};{order};{family};{genus};{species};{strain|subspecies|no rank}" > "$OUT_TAX"

# ID変換表の作成
cut -f1,8,9 "$IN_TYPE" \
| awk -F"\t" '{ gsub(/^ +/, "", $2); sub("\\[","-",$2); sub("\\]","-",$2); sub("strain=","",$3); gsub(/ +/,"_", $3); \
if (match($2, /^(Candidatus +)?[A-Za-z\[\]-]+ +[a-z\[\]]+( +subsp\. +[a-z0-9\-]+)?/, arr)) { \
  name = arr[0] ; gsub(/ +/, "_", name) ; print $1 "\t" name "_" $3 "_" $1
    } else { print $1 "\tUNKNOWN_" $3 "_" $1 ; } }' > "$OUT_TAB"

```

- 以下のようにdeleteやnot foundの警告がでなければO.K.
```txt
22:28:36.509 [WARN] taxid 3238480 not found
22:28:36.509 [WARN] taxid 3133966 was deleted
```
## blast db / skani-db
```bash
find ./genomes -name "./genomes*.fna.gz" \
| while read -r gca ; do pfx=$(basename $gca | cut -f1,2 -d"_"); gunzip -c $gca ; done \
| awk -v pfx="$pfx" '/^>/{sub(">",">" pfx "_",$1);print}!/^>/{print}' | gzip > ./blast_db/genki_250508.fna.gz

gunzip -c genki_250508.fna.gz | makeblastdb -in - -out genki250508 -title genki_250508 -dbtype nucl -parse_seqids -hash_index &

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

