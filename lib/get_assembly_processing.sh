#!/bin/bash

###########################################################################################
# filter_summary : assembly summaryデータを指定条件でフィルタリングする
# - ここでは細菌・古細菌のタイプマテリアル由来のデータを抽出する条件を指定している
# - filter_summary <input.txt> [<output.txt>] 
###########################################################################################
function filter_summary () {
    if [ $# -eq 0 ]; then echo "Usage: filter_assembly <input> [<output>]" >&2 ; return 1 ; fi
    local sum_gca=$1 
    local sum_gca_type=${2:-"assembly_summary_genbank_type_prok.txt"} 
    # カラム指定
    group_col=25 ; type_col=22
    # アセンブリサマリーの処理
    awk -F'\t' -v gc=$group_col -v rc=$type_col '$gc ~ /^(bacteria|archaea)$/ && $rc ~ /^assembly from (synonym )?type material$/ ' "$sum_gca" > "$sum_gca_type"

}

###########################################################################################
# compare_assembly : 過去データと追加データの比較
# - ここでは細菌・古細菌のタイプマテリアル由来のデータを抽出する条件を指定している
# - filter_summary <input.txt> [<output.txt>] 
###########################################################################################
function compare_assembly () {
    local asm_old=$1
    local asm_new=$2

    # unique array
    function uqarray () { 
        # usage: uqarray array1 array2 mode
        # mode: 1: Only1, 2: Only2, 3: Both
        local -n A1_REF=$1 ; local -n A2_REF=$2 ; local MODE=$3
        local BOTH ONLY1 ONLY2

        mapfile -t BOTH < <(printf "%s\n%s\n" "${A1_REF[@]}" "${A2_REF[@]}" | sort | uniq -d)
        mapfile -t ONLY1 < <(printf "%s\n%s\n" "${A1_REF[@]}" "${BOTH[@]}" | sort | uniq -u)
        mapfile -t ONLY2 < <(printf "%s\n%s\n" "${A2_REF[@]}" "${BOTH[@]}" | sort | uniq -u)

        if [[ $MODE == 1 ]]; then 
            printf "%s\n" "${ONLY1[@]}"
        elif [[ $MODE == 2 ]]; then
            printf "%s\n" "${ONLY2[@]}"
        elif [[ $MODE == 3 ]]; then
            printf "%s\n" "${BOTH[@]}"
        else
            echo "[ERROR] Invalid mode: $MODE" >&2
        fi
    } 

    # Extract matching lines for given IDs
    function get_id_stream() {
        local ids_file="$1"
        local input_file="$2"

        if [[ ! -f "$ids_file" || ! -f "$input_file" ]]; then
            echo "[ERROR] File(s) not found: $ids_file or $input_file" >&2
            return 1
        fi
        awk 'NR==FNR {ids[$1]; next} $1 in ids' "$ids_file" "$input_file"
    }

    # Check if input files exist
    if [[ ! -f "$asm_old" || ! -f "$asm_new" ]]; then
        echo "[ERROR] Input files not found: $asm_old or $asm_new" >&2
        return 1
    fi

    # Extract IDs
    # shellcheck disable=SC2034
    mapfile -t old_id < <(cut -f1 "$asm_old")
    # shellcheck disable=SC2034
    mapfile -t new_id < <(cut -f1 "$asm_new")

    # Temporary files for IDs
    tmp_old_ids=$(mktemp)
    tmp_new_ids=$(mktemp)

    # Get removed and added assemblies
    uqarray old_id new_id 1 > "$tmp_old_ids"
    uqarray old_id new_id 2 > "$tmp_new_ids"

    get_id_stream "$tmp_old_ids" "$asm_old" > assembly_removed.tsv
    get_id_stream "$tmp_new_ids" "$asm_new" > assembly_addition.tsv

    # Clean up temporary files
    rm -f "$tmp_old_ids" "$tmp_new_ids"
}

###########################################################################################
# geturl_assembly : アセンブリサマリからアセンブリ及びmd5のurlリンクを取得する
###########################################################################################
function geturl_assembly () {
    if [ $# -eq 0 ]; then echo "Usage: geturl_assembly <input> [<output>]" >&2 ; return 1 ; fi
    local sum_gca=$1
    local out_dir=${2:-'log'}
    local ftp_genome="${out_dir}/ftp_genome"    # ダウンロード用のゲノムファイルリンク
    local ftp_md5="${out_dir}/ftp_md5"          # MD5チェックサムリンク

    if [[ -d "$out_dir" ]] ; then echo "[ERROR] $out_dir already exists" ; return 1; else mkdir -p "$out_dir" ; fi
    
    # カラム指定
    group_col=25 ; type_col=22 ; ftp_col=20

    # アセンブリサマリーの処理
    awk -F'\t' -v gc=$group_col -v rc=$type_col -v ftp=$ftp_col '$gc ~ /^(bacteria|archaea)$/ && $rc ~ /^assembly from (synonym )?type material$/ \
    { split($ftp, arr, "/"); asm_name=arr[length(arr)]; glnk=$ftp "/" asm_name "_genomic.fna.gz"; md5lnk=$ftp "/md5checksums.txt"; print glnk > "'$ftp_genome'"; print md5lnk > "'$ftp_md5'"; }' "$sum_gca"
    return 0
}

###########################################################################################
# get_gca_genomic: 
# - urlリンクファイルをもとにしてゲノムデータをダウンロードする
# - ゲノムデータがリポジトリに存在しない場合、wgetでエラーを拾う
# get_md5_genomic
# - urlリンクファイルをもとにmd5データをダウンロードした後、ゲノムデータのみ抽出
# - 登録がsubmitterにより取り消された等の場合でもmd5checksums.txtには記述が存在するため、
#   md5checksums.txtをダンロードした後にその内容をチェックする
###########################################################################################
function get_gca_genomic () {
    # Usage: get_gca_genomic ./ftp_genome ./genomes
    local FTP_GCA=$1
    local OUT_DIR=${2:-'genomes'}

    if [[ ! -f "$FTP_GCA" ]] ; then echo "[ERROR] $FTP_GCA does not exists" ; return 1; fi
    if [[ -d "$OUT_DIR" ]] ; then echo "[ERROR] $OUT_DIR already exists" ; return 1; else mkdir -p "$OUT_DIR" ; fi
    
    # エラーログをクリア
    local error_gca ; error_gca=$(dirname "${OUT_DIR}")/wget_gca.err
    : > "$error_gca"

    # ファイルごとにダウンロードを試行し、エラーを記録
    while IFS= read -r url; do
        if wget -q -P "$OUT_DIR" "$url"; then
            :
        else
            echo "FAIL: $url" >> "$error_gca"
        fi
    done < "$FTP_GCA"

}

function get_md5_genomic () {
    # Usage: get_md5_genomic <ftp_md5_list> <merged_md5>
    # Error: get_md5_genomic ./ftp_md5 ./genomes/md5_genomes
    local FTP_MD5=$1
    local OUT_DIR=${2:-'genomes'}

    # 入出力チェック
    if [[ ! -f "$FTP_MD5" ]] ; then echo "[ERROR] $FTP_MD5 does not exists" ; return 1; fi
    if [[ ! -d "$OUT_DIR" ]] ; then mkdir -p "$OUT_DIR" ; fi
    local out_md5="${OUT_DIR}/md5_genomes"
    local error_md5 ; error_md5=$(dirname "$OUT_DIR")/wget_md5.err 

    # エラーログをクリア
    : > "$error_md5"
    
    # MD5リンクの内容を取得し、該当ファイルのみ抽出
    while IFS= read -r url; do
        if wget -q -O - "$url" > tmp_md5_output; then
            if ! grep "_genomic.fna.gz" tmp_md5_output | grep -v -e "_rna_" -e "_cds_" >> "$out_md5"; then
                echo "NO MATCH: $url" >> "$error_md5"
            fi
        else
            echo "FAIL: $url" >> "$error_md5"
        fi
        rm tmp_md5_output 
    done < "$FTP_MD5" > "$out_md5"    
}
