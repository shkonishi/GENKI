#!/usr/bin/Rscript

### Contents ###
# 1. Input arguments
#  1.1. Check arguments
#  1.2. Print arguments
# 2. Main routine
#

### NOTE
# Advanced searchでダウンロードしたcsvからidを取得する
##  e.g. type=yes & Genome database=ncbi
# in_f <- "~/projects2/GENKI_data/dat/advsearch_bacdive_Genome_Sequence_data_2024-12-30.csv"
# dat <- read.csv(in_f, skip = 3, header = T, comment.char="")



### Description ###
# BacDive-APIにアクセスして、タイプのゲノム関連のメタデータを取得する

###################################################################
# ユーザー定義関数:オブジェクトの存在を確認して、存在しなけれ終了
###################################################################
exists_obj <- function(obj_name) {
  if (!exists(obj_name, envir = .GlobalEnv)) {
    stop(paste("Error: Object", obj_name, "does not exist."))
  }
}
###################################################################
# ユーザー定義関数:指定ベクトルを任意の数で分割したリストを作成
###################################################################
chunkr <- function (x, n, r = FALSE) 
{
  if (r) {
    nchnk <- ceiling(seq_along(x)/n)
    split(x, nchnk)
  }
  else {
    nchnk <- ceiling(length(x)/n)
    split(x, cut(seq_along(x), nchnk, labels = FALSE))
  }
}
###################################################################


# 1. Input arguments
argv <- commandArgs(TRUE)

#  1.1. Check arguments
if (length(argv) != 3) {
  msg_u <- "# [USAGE]: Rscript this_script.r <input> <API_ID> <API_PASS>"
  msg_i <- "# [INFO]: Rscript this_script.r advsearch_bacdive_2024-08-02.csv XXXXX@gmail.com PASSWORD123"  
  message(paste(msg_u, msg_i, sep = "\n"))
  q()
}
if (!file.exists(argv[1])) {
  msg <- paste("[ERROR]", argv[1], "does not exists.")
  message(msg) ; q()
} else {
  in_f <- argv[1]
}

if (is.na(argv[2])) {
  msg <- paste("[ERROR]", argv[2], "must be user-name for BacDive-API.")
  message(msg1) ; q()
} else {
  UNAME <- argv[2]
}

if (is.na(argv[3])) {
  msg1 <- paste("[ERROR]", argv[3], "must be password for BacDive-API.")
  message(msg) ; q()
} else {
  PASWD <- argv[3]
}

out_f1 <- paste0("bacdive_summary_", Sys.Date(), ".tsv")
out_f2 <- paste0("bacdive_genbank_", Sys.Date(), ".tsv")

#  1.2. Print arguments
msg1 <- message(paste("input	[", in_f, "]"))
msg2 <- message(paste("USER 	[", UNAME, "]"))
msg3 <- message(paste("PASSWORD [", "XXXXX", "]"))
msg4 <- message(paste("output1  [", out_f1, "]"))
msg5 <- message(paste("output2  [", out_f2, "]"))
message(paste(msg1, msg2, msg3, msg4, sep = "\n"))

# 2. Main : BacDive
## 2.2. BacDiveのリストを読込
# BacDive -> statistics -> Type strains -> Taxonomy -> Download table as CSV
dsmz_dat <- read.csv(in_f, TRUE, skip = 2)
chnk <- chunkr(dsmz_dat$ID, n = 100)

# Open BacDive
bacdive <- BacDive::open_bacdive(username = UNAME, password = PASWD)
exists_obj("bacdive")
## 指定idのデータをfetch & データフレームに変換(idの並びは崩れる)
## 時間がかかる、また途中で失敗する場合がある(curlのエラー)
bd_dats <- lapply(seq_along(chnk), function(i){
  bd <- BacDive::fetch(bacdive, chnk[[i]])
  as.data.frame(bd)
})
exists_obj("bd_dats")
bd_dat <- Reduce(rbind, bd_dats)

## カテゴリごとのリストを作成(とりあえず以下の要素のみ)
dat_gen <- bd_dat$General
dat_tax <- bd_dat$`Name and taxonomic classification`
dat_seq <- bd_dat$`Sequence information`
dat_safe <- bd_dat$`Safety information`

#bd_dat$Morphology # Colony color, forms multicellular complex, incubation period, gram stain etc
#dat_growth <- bd_dat$`Culture and growth conditions` 
#bd_dat$`Physiology and metabolism` 
#bd_dat$`Isolation, sampling and environmental information`
#bd_dat$`Genome-based predictions`
#bd_dat$`External links"
#bd_dat$Reference

########################
# General
########################
bacid_gen <- sapply(dat_gen, function(x)x$`BacDive-ID`) ## BacDive-ID
desc_gen <- sapply(dat_gen, function(x)x$description) ## Description

########################
# Name and taxonomic classification
########################
## LPSN lineage
lineage <- as.data.frame(Reduce(rbind, lapply(dat_tax, function(x){
    lpsn <- x$LPSN
    domain <- ifelse(is.null(lpsn$domain), NA, lpsn$domain)
    phylum <- ifelse(is.null(lpsn$phylum), NA, lpsn$phylum)
    class <- ifelse(is.null(lpsn$class), NA, lpsn$class)
    order <- ifelse(is.null(lpsn$order), NA, lpsn$order)
    family <- ifelse(is.null(lpsn$family), NA, lpsn$family)
    genus <- ifelse(is.null(lpsn$genus), NA, lpsn$genus)
    species <- ifelse(is.null(lpsn$species), NA, lpsn$species)
    full_sfn <- ifelse(is.null(lpsn$`full scientific name`), NA, lpsn$`full scientific name`)

    setNames(c(domain, phylum, class, order, family, genus, species, full_sfn), 
        c("domain","phylum","class","order","family","genus","species","full_scientific_name"))
})))
## synonym
synonym <- sapply(dat_tax, function(x) {
  lpsn <- x$LPSN
  if(is.null(lpsn$synonyms)){
    NA
  }else{
    syn <- unlist(lpsn$synonyms)
    paste(syn[names(syn) %in% "synonym"], collapse = ";")
  }
})

## strain
strain <- sapply(dat_tax, function(x){ 
    v <- x$`strain designation`
    ifelse(is.null(v),NA,v)
})
########################
# Genome 
# データ構造確認 sapply(dat_gseq, mode) %>% table() # logical かリスト
# dat_seq[[224]][[1]]$accession
# dat_seq[[116]]$accession
# dat_seq[[1]]$`Genome sequences`$accession
# dat_seq[[2]]$`Genome sequences`[[1]]$accession
# unlist(dat_seq[[2]])[grep("Genome sequences.accession", names(unlist(dat_seq[[2]])))]
########################
dat_gseq <- lapply(dat_seq, function(x){
  if(!is.null(x$accession)){
    x
  } else if (!is.null(x$`Genome sequences`)) {
    x$`Genome sequences`
  } else if (is.null(x$accession) && is.null(x$`Genome sequences`) && !is.null(x[[1]]$accession)) {
    x
  } else {
    NA
  }
})
gdat <- as.data.frame(Reduce(rbind, lapply(seq_along(dat_gseq), function(i){
    x<- dat_gseq[[i]] ; # print(i)
    if(all(is.na(x)) || is.null(x)){
        accs_gs <- NA ; desc_gs <- NA ; lev_gs <- NA ; db_gs <- NA
    } else if(!is.null(x) && !is.null(x$accession)){
      accs_gs <- x$accession
      desc_gs <- x$description
      lev_gs  <- ifelse(is.null(x$`assembly level`),NA, x$`assembly level`)
      db_gs <- x$database
    } else if(!is.null(x) && is.null(x$accession)){
      accs_gs <- paste(sapply(x, function(v)v$accession), collapse = ";")
      desc_gs <- paste(sapply(x, function(v)v$description), collapse = ";")
      lev_gs <- paste(sapply(x, function(v)v$`assembly level`), collapse = ";")
      db_gs <- paste(sapply(x, function(v)v$database), collapse = ";")
    }
    setNames(c(accs_gs, desc_gs, lev_gs, db_gs),c("accession", "genome_description", "assembly_level", "database" ))
})))
######################
# Biosafety level
######################
biosafety_lv <- sapply(seq_along(dat_safe), function(i){
    x <- dat_safe[[i]]
    if(all(is.na(x))){
        NA
    } else {
        paste(unlist(x)[names(unlist(x)) %in% "biosafety level"], collapse=";")
    }
})

# 2.5. 要約データを保存しておく
dat <- data.frame(bacdive_id=bacid_gen, description=desc_gen, lineage, synonym, strain, biosafety_lv, gdat)
dat <- dat[order(dat$description),]
write.table(dat, out_f1, sep="\t", quote=F, row.names=F, col.names=T)


## GenBankアクセッションを取り出す(17284)
gca <- sub("\\..*$","", unlist(sapply(strsplit(dat$accession, ";"), function(x)grep("GCA", x, value=T))))
## BacDiveに記載されたGenBankアクセッションのなかで、GenBankのほうにデータが存在しない登録を確認( "GCA_900182945" )
#gca[which(is.na(match(gca, sub("\\..*$","",asmb$X.assembly_accession))))]

## genomeデータが存在するbacdiveの登録
genome_bacid <- setNames(
  stack(setNames(
    lapply(strsplit(dat$accession, ";"), function(x){
      v <- sub("\\..*$","", grep("GCA", x, value=T))
      if(length(v)==0){ NA }else{ v }
      }), 
    dat$bacdive_id
    )),
  c("gca", "bacdive_id")
  )
## Genbankのgenomeデータが存在する13258種、17284登録 (同種異登録含む)
bacid_gca <- genome_bacid[!is.na(genome_bacid$gca),]

## ユニークなGCAアクセッションでBacDiveデータを整形する
dat_gca <-merge(bacid_gca, dat, by = "bacdive_id", all.x = TRUE)
dat_gca <- dat_gca[order(dat_gca$description), c(2,1,3:ncol(dat_gca))]
write.table(dat_gca, out_f2, sep="\t", quote=F, row.names=F, col.names=T)


