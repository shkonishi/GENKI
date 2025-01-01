# GENKI
原核生物の基準株ゲノムデータベースを構築する。ゲノム情報を基準として、以下のように異なるデータベースから取得したデータを紐付けて管理する。
- **ゲノムデータ**: NCBI GenBank
- **系統情報**: NCBI Taxonomy
- **表現型データ**: BacDive

## **依存**
- taxonkit
- BacDive (R-package of API Client)

## **データ構造**
```sh
.
├── dat
│   ├── assembly_summary_genbank_type_prok_241229.txt # 指定条件(原核生物のみ、type material)でフィルタしたアセンブリサマリ
│   ├── assembly_summary_genbank_type_prok_lineage_241229.txt # 系統情報
│   └── id_lookup.tsv # IDと種名の対応表
├── genomes
│   ├── GCA_000003925.1_ASM392v1_genomic.fna.gz
│   ├── GCA_000006685.1_ASM668v1_genomic.fna.gz
│   ├── GCA_000006945.2_ASM694v2_genomic.fna.gz
│   ├── GCA_000006985.1_ASM698v1_genomic.fna.gz
│   ├── GCA_000007025.1_ASM702v1_genomic.fna.gz
│   ├── GCA_000007185.1_ASM718v1_genomic.fna.gz
│   └── md5_genomes
└── log
    ├── error_gca.log  # wgetのエラーログ
    ├── error_md5.log  # wgetのエラーログ
    ├── ftp_genome     # ゲノムデータURLリスト
    ├── ftp_md5        # md5ファイルURLリスト
    ├── res_md5sum.txt # md5sumの結果 
    ├── res_md5sum_err.txt  # md5sumのエラー
    └── res_md5sum_fail.txt # ハッシュ値が不一致のファイル名

```

## **データダウンロード**
### **ユーザー定義関数をロード**
- 以下の一連の工程は、スクリプトにする予定だが、現時点では関数を個別に実行する仕様
```bash
source ../lib/get_assembly_processing.sh
```

### **アセンブリサマリのダウンロード**
- アセンブリサマリに記述されている情報について下記リンクのファイルから把握しておく
    https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/README_assembly_summary.txt

```bash
wget -c https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt
```

### **細菌・古細菌のみ、かつタイプのみの情報を抽出**
- `relation_to_type_material`列にタイプとの関連情報が、`group`列にタクソノミーグループの情報(bacteria, archaea, その他)が記述されている。
- `assembly from type material`または`assembly from synonym type material`を選択する。ただし、それらの登録の中にはメタゲノム由来も含まれる事に注意。
```bash
VER='241229'
IN_SUM="assembly_summary_genbank.txt"
OUT_SUM="./dat/assembly_summary_genbank_type_prok_${VER}.txt"

mkdir -p ./dat
filter_summary "$IN_SUM" "$OUT_SUM"
```
- 追加登録のみをダウンロードする場合
- 追加登録の情報情報を抽出したassembly_addition.txtが作成される。 

```bash
IN_SUM_OLD="./dat/assembly_summary_genbank_type_prok_XXYYZZ.txt"
IN_SUM_NEW="./dat/assembly_summary_genbank_type_prok_${VER}.txt"
compare_assembly "$SUM_TYPE_OLD" "$SUM_TYPE_NEW"
```
### **URLリンクを作成**
```bash
# 目的のデータURLファイル作成
geturl_assembly "$OUT_SUM"
# 追加登録のみのURLファイル作成
geturl_assembly assembly_addition.txt
```
### **データダウンロード**
- 3.で作成されたURLリストファイルを用いてデータをダウンロードする。
- error_gca.log内にダウンロードを失敗した登録が記録されるので、再ダウンロード

```bash
get_gca_genomic ./log/ftp_gca ./genomes &
get_md5_genomic ./log/ftp_md5 ./genomes/md5_genomes &

# md5sumを実行
cd ./genomes
md5sum -c md5_genomes > ../res_md5sum.txt 2> ../res_md5sum_err.txt
# DL失敗したものを確認
cd ../
grep -v 'OK' ./res_md5sum.txt > ./res_md5sum_fail.txt
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