# GENKI
原核生物の基準株ゲノムデータベースを構築する。ゲノム情報を基準として、以下のように異なるデータベースから取得したデータを紐付けて管理する。
- **ゲノムデータ**: NCBI GenBank
- **系統情報**: NCBI Taxonomy
- **表現型データ**: BacDive

## **依存**
- ncbi-genome-download
- taxonkit
- BacDive (R-package of API Client)

## **ゲノムデータダウンロード**
### ncbi-genome-downloadを使う場合
```bash
# dry-runを実行してIDを取得しておく
VER='241229'
ncbi-genome-download -s genbank -F fasta -M type --dry-run bacteria,archaea > gca_dl_${VER}.txt

# 新規にDLする場合
ncbi-genome-download -s genbank -F fasta -M type -o "GENVI_${VER}" -m METADATA_TABLE bacteria,archaea > gca_dl_${VER}.txt

# 追加登録されたもののみをDLする場合
ADD='250105'
ncbi-genome-download -s genbank -F fasta -M type --dry-run bacteria,archaea > gca_dl_${ADD}.txt

# 除外された登録及び追加された登録を確認

# 追加された登録をDL

```

### genki_download.shを使う場合(現時点では関数を個別に実行する仕様)  
#### **ダウンロード実行**
- NCBI FTPサーバーリンクからサマリーデータをダウンロードする
- アセンブリサマリに記述されている情報について下記リンクのファイルから把握しておく
    https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/README_assembly_summary.txt
- `relation_to_type_material`列にタイプとの関連情報が、`group`列にタクソノミーグループの情報(bacteria, archaea, その他)が記述されている。
- `assembly from type material`または`assembly from synonym type material`を選択する。ただし、それらの登録の中にはメタゲノム由来も含まれる事に注意。


```bash
# スクリプト仕様(予定)
# genki_download.sh 
#-g bacteria,archae \
#-r type,synonym \
#-o out_dir "$IN_SUM" \

# genki_download.sh -g bacteria,archae -r type,synonym -o out_dir "$IN_SUM"

# ユーザー定義関数をロード
source ../lib/get_assembly_processing.sh
# アセンブリサマリのダウンロード
VER='241229'
## GenBank
wget -c -O "assembly_summary_genbank_${VER}.txt" https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt
## RefSeq
wget -c -O "assembly_summary_refseq_${VER}.txt" https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt

# 指定条件のゲノムデータをダウンロードする
IN_SUM="assembly_summary_genbank_${VER}.txt"
OUT_DIR="out_genki_${VER}"
OUT_DAT_DIR="${OUT_DIR}/dat"
OUT_LOG_DIR="${OUT_DIR}/log"
OUT_ASMB_DIR="${OUT_DIR}/genomes"
out_md5="${OUT_ASMB_DIR}/md5_genomes"
out_md5chk="${OUT_DIR}/res_md5sum.txt"
out_md5fail="${OUT_DIR}/res_md5sum_err"

# サマリデータフィルタリング
mkdir -p "${OUT_DAT_DIR}"
OUT_SUM="${OUT_DAT_DIR}/assembly_summary_genbank_type_prok_${VER}.txt"
filter_summary "$IN_SUM" "$OUT_SUM"

# URLファイル作成
geturl_assembly "$OUT_SUM" "$OUT_LOG_DIR"

# データダウンロード
get_gca_genomic "${OUT_LOG_DIR}/ftp_genome" "$OUT_ASMB_DIR" &
wait
get_md5_genomic "${OUT_LOG_DIR}/ftp_md5" "$OUT_ASMB_DIR" &
wait
# md5sumを実行
(cd "$OUT_ASMB_DIR" && md5sum -c ./md5_genomes) > "$out_md5chk" 2> "$out_md5fail"

# DL失敗したものを確認

```

- 追加登録のみをダウンロードする場合
- 追加登録の情報情報を抽出したassembly_addition.txtが作成される。 

```bash
IN_SUM_OLD="./assembly_summary_genbank_241229.txt"
IN_SUM_NEW="./assembly_summary_genbank_250104.txt"

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
- taxonkitを使ってlineageを取得&整形
- idとtaxonomyの変換表を作成しておく。
    ゲノムデータを参照データとした解析(ANI,コアジーン系統解析, cgMLST)に利用する。

```bash
IN_TYPE='./dat/assembly_summary_genbank_type_prok.txt'
OUT_TAX='./dat/assembly_summary_genbank_type_prok_lineage.txt'
OUT_TAB='./dat/id_lookup.tsv'

cut -f1,6 $IN_TYPE | taxonkit lineage -i 2 | taxonkit reformat -i 3 -f "{k};{p};{c};{o};{f};{g};{s};{t}" | cut -f1,2,4 > "$OUT_TAX"

# idルックアップデーブルの作成
cat "$OUT_TAX" | awk -F"\t" '{split($3,arr,";"); gsub(" ","_", arr[7]); print $1"\t"arr[7]"_"$1;}' > "$OUT_TAB"
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