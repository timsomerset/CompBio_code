library("dplyr")

filter <- function(out){
  data <- data.frame(
    GenEx = double(), Type = character(),
    S = character(), begin = double(),
    end = double(), length = double(),
    read_frame = double(), netphase = double(),
    init_accept_score = double(), donor_term_score = double(),
    cod_reg_score = double(), p_exon = double(),
    exon_score = double()
  )
  seqid <- unlist(strsplit(out[3], " "))[2]
  if ("NO EXONS/GENES PREDICTED IN SEQUENCE" %in% out){
    data <- c(seqid, rep(NA, 15))
  }else{
    zone1 <- c(grep("-----", out), grep("Predicted peptide sequence", out))
    # extract annotation statistics --------------------------------------------
    for (i in (zone1[1]+1):(zone1[2]-1)){
      # only include non-empty rows
      if (nzchar(out[i])){
        out_vec <- unlist(strsplit(out[i], " "))
        # remove spaces and NA's
        out_vec <- out_vec[-which(out_vec == ""|is.na(out_vec))]
        # input NA for entries that are missing, making all vectors equal length
        if (length(out_vec) != 13){
          out_vec <- c(out_vec[1:6], rep(NA, 6), out_vec[7])
          out_vec <- as.list(out_vec)
          out_vec[c(1,4:6,13)] <- as.numeric(unlist(out_vec[c(1,4:6,13)]))
        }else{
          out_vec <- as.list(out_vec)
          out_vec[c(1,4:13)] <- as.numeric(unlist(out_vec[c(1,4:13)]))
        }
        data[nrow(data)+1,] <- out_vec
      }
    }
    # extract annotation regions (peptide and nucleotide sequences)-------------
    zone_pep <- c(which(grepl("predicted_peptide", out)))
    zone_nucl <- which(grepl("predicted_CDS", out))
    peptides <- c()
    nucleotides <- c()
    for (i in 1:length(zone_pep)){
      pep <- paste(out[(zone_pep[i]+1):(zone_nucl[i]-2)], collapse = "")
      peptides[i] <- pep
      if(i < length(zone_nucl)){
        nucl <- paste(out[(zone_nucl[i]+1):(zone_pep[i+1]-2)], collapse = "")
        nucleotides[i] <- nucl
      }else{
        nucleotides[length(nucleotides)+1] <- paste(out[(max(zone_nucl)+1):length(out)],collapse = "")
      }
    }
    genes <- floor(data[,1])
    pep_seq <- peptides[genes]
    nucl_seq <- nucleotides[genes]
    scaf <- rep(seqid, nrow(data))
    data <- cbind(scaf, data, pep_seq, nucl_seq)
  }
  data
}
genscanner <- function(direc_in, direc_out, file_out){
  directory <- direc_in
  dir <- dir(directory)
  data <- data.frame(
    GenEx = double(), Type = character(),
    S = character(), begin = double(),
    end = double(), length = double(),
    read_frame = double(), netphase = double(),
    init_accept_score = double(), donor_term_score = double(),
    cod_reg_score = double(), p_exon = double(),
    exon_score = double())
  # ***input must be directory name, with path from work directory ***
  # only get names of fasta files, ignore subdirectories
  # dir_fasta <- dir[grep("\\.fasta", dir)]
  # make database to track progress
  track <- data.frame("Scaffold" =  dir, "Genscanned" = integer(length(dir)))

  for (i in 1:length(dir)){
    file <- dir[i]
    size <- file.info(paste0(directory, "/", file))$size
    if (size < 2000000){
      # scaf_index used to save file with correct index in name
      scaf_index <- unlist(strsplit(unlist(strsplit(file, "_"))[2], "\\."))[1]
      direc <- paste0(directory, "/", file)
      system(paste0(
        "/local/data/public/eadc2/Genomics_1/programs/genscan/genscan ",
        "/local/data/public/eadc2/Genomics_1/programs/genscan/HumanIso.smat ",
        direc, " -cds",
        " > /local/data/mphilcompbio/2023/ts932/g1_assignment2/", direc_out,"/scaf_",
        scaf_index, ".out"))

      out <- readLines(paste0(direc_out, "/scaf_", scaf_index, ".out"))

      data <- rbind(data, filter(out))
      write.table(data, paste0(file_out, ".tsv"))
      track[i,2] <- 1
      write.table(track, "genscan_tracker.tsv")

    }else{
      # *** for files larger than 2mb, we split and put into new directory ***
      scaf <- readLines(paste0(directory, "/", file))
      head <- scaf[1]
      scaf_id <- unlist(strsplit(file, "\\."))[1]
      # create new subdirectory for specific scaffold
      new_directory <- paste0(directory, "/", scaf_id, "_big")
      dir.create(new_directory)

      n_pieces <- ceiling(size/2000000) + 1
      l_piece <- ceiling(length(scaf)/n_pieces)
      piece_divide <- c(seq(2,length(scaf), l_piece), length(scaf))
      sapply(seq(1:(length(piece_divide)-1)),
             function(x) write(c(head, scaf[piece_divide[x]:piece_divide[x+1]-1]), paste0(new_directory, "/", scaf_id, "-", x, ".fasta")))


      for (j in 1:length(dir(new_directory))){
        # *** This is the same genscan process, but helps with tracking ***
        file <- dir(new_directory)[j]
        scaf_index <- unlist(strsplit(unlist(strsplit(file, "_"))[2], "\\."))[1]
        direc <- paste0(new_directory, "/", file)
        system(paste0(
          "/local/data/public/eadc2/Genomics_1/programs/genscan/genscan ",
          "/local/data/public/eadc2/Genomics_1/programs/genscan/HumanIso.smat ",
          direc, " -cds",
          " > /local/data/mphilcompbio/2023/ts932/g1_assignment2/", direc_out, "/scaf_",
          scaf_index, ".out"))
        out <- readLines(paste0(direc_out, "/scaf_", scaf_index, ".out"))

        data <- rbind(data, filter(out))
        write.table(data, paste0(file_out, ".tsv"))
      }
      track[i,2] <- 1
      write.table(track, "genscan_tracker.tsv")
    }
  }
  data
}
extract_scafid <- function(fasta_path){
  file <- readLines(fasta_path)
  header <- file[1]
  unlist(strsplit(unlist(strsplit(file, " "))[1], ">"))[2]
}
genchecker <- function(direc_in, genscanner_out){
  direc_files <- dir(direc_in)
  direc_fasta <- direc_files[grep("\\.fasta", direc_files)]
  fasta_paths <- sapply(direc_fasta, function(x) paste0(direc_in, "/", x))
  all_id <- sapply(fasta_paths, extract_scafid)
  in_genscanner_out <- as.numeric(all_id %in% genscanner_out$scaf)
  cbind(all_id, in_genscanner_out)
}

genpost <- function(genout){
  # *** Input is a file with genscan output ***
  genout <- read.table(genout)
  genout$GenEx <- floor(genout$GenEx)
  # database for all exons across all scaffolds
  simple_out <- data.frame( scaf = character(), geneid = double(),
                             S = character(), begin = double(),
                             end = double(), length = double(),
                             read_frame = double(), netphase = double(),
                             init_accept_score = double(), donor_term_score = double(),
                             cod_reg_score = double(), p_gene = double(),
                             gene_score = double()
  )
  # for each unique scafid in genout
  for (i in unique(genout$scaf)){
    long_scaf <- genout[which(genout$scaf == i),]
    # a database that'll contain ALL exons for a single scaffold (one exon per line)
    simple_scaf <- data.frame( scaf = character(), geneid = double(),
                               S = character(), begin = double(),
                               end = double(), length = double(),
                               read_frame = double(), netphase = double(),
                               init_accept_score = double(), donor_term_score = double(),
                               cod_reg_score = double(), p_gene = double(),
                               gene_score = double(), extra = character()
    )

    # THIS IS WHERE THE WHILE LOOP BEGINS
    # j is the row number on which the first instance of our exon is seen
    j <- 1
    counter <- 1
    # for each unique exonid
    while (j < nrow(long_scaf)){
      geneid <- long_scaf$GenEx[j]
      rows_same_gene<- which(long_scaf$GenEx[j:nrow(long_scaf)] == geneid) + j-1
      closest_same_gene <- c(j)
      # we take the first consecutive sequence
      for (k in 2:length(rows_same_gene)){
        if (rows_same_gene[k] == (closest_same_gene[k-1]+1)){
          closest_same_gene[k] <- rows_same_gene[k]
        } else{
          break
        }
      }

      # Process data for adding to simple_scaf dataset
      data_gene <- long_scaf[closest_same_gene,]
      structural <- which(data_gene$Type == "PlyA" | data_gene$Type == "Prom")
      if (length(structural)==0){
        simple_gene <- c(scaf = i,
                         geneid = counter,
                         S = data_gene$S[1],
                         begin = data_gene$begin[1],
                         end = data_gene$end[nrow(data_gene)],
                         length = sum(data_gene$length),
                         read_frame = mean(data_gene$read_frame),
                         netphase = mean(data_gene$netphase),
                         init_accept_score = min(data_gene$init_accept_score),
                         donor_term_score = min(data_gene$donor_term_score),
                         cod_reg_score = min(data_gene$cod_reg_score),
                         p_gene = prod(data_gene$p_exon),
                         gene_score = min(data_gene$exon_score),
                         extra = data_gene$nucl_seq[1]
        )
      }else{
        simple_gene <- c(scaf = i,
                         geneid = counter,
                         S = data_gene$S[1],
                         begin = data_gene$begin[1],
                         end = data_gene$end[nrow(data_gene)],
                         length = sum(data_gene$length),
                         read_frame = mean(data_gene$read_frame[-structural]),
                         netphase = mean(data_gene$netphase[-structural]),
                         init_accept_score = min(data_gene$init_accept_score[-structural]),
                         donor_term_score = min(data_gene$donor_term_score[-structural]),
                         cod_reg_score = min(data_gene$cod_reg_score[-structural]),
                         p_gene = prod(data_gene$p_exon[-structural]),
                         gene_score = min(data_gene$exon_score[-structural]),
                         extra = data_gene$nucl_seq[1]
        )
      }

      simple_gene <- as.list(simple_gene)
      simple_gene[c(2,4:13)] <- as.numeric(unlist(simple_gene[c(2,4:13)]))
      simple_scaf[nrow(simple_scaf)+1,] <- simple_gene
      j <- max(closest_same_gene)+1
      counter <- counter +1
    }
    simple_out <- rbind(simple_out, simple_scaf)
  }
  simple_out
}

gffinator <- function(genout, gff_out){
  genpost <- genpost(genout)
  gff <- cbind(seqname = genpost$scaf,
               source = rep("Genscan", nrow(genpost)),
               feature = rep("coding gene", nrow(genpost)),
               start = genpost$begin,
               end = genpost$end,
               score = genpost$cod_reg_score,
               strand = genpost$S,
               frame = genpost$read_frame,
               attribute = genpost$extra
  )
  directory <- paste0(gff_out, ".gff")
  write.table(gff, directory)
  gff
}






