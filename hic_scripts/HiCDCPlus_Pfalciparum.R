#remotes::install_github("thephilross/bsgenome.pf.v24")

library(dplyr)
get_chrs_v2 <- function(gen = "Pfalciparum", gen_ver = "v24") {
    genome <- paste("BSgenome.", gen, '.PlasmoDB.', gen_ver, sep = "")
    library(genome, character.only = TRUE)
    chrs <- GenomeInfoDb::seqnames(get(gen))
    chrs <- chrs[!chrs %in% c("PFC10_API_IRAB","M76611")]
    return(chrs)
}

get_chr_sizes_v2 <- function(gen = "Pfalciparum", gen_ver = "v24", chrs = NULL) {
    genome <- paste("BSgenome.", gen, '.PlasmoDB.', gen_ver, sep = "")
    library(genome, character.only = TRUE)
    if (is.null(chrs)) {
        chrs <- get_chrs_v2(gen, gen_ver)
    }
    sizes <- GenomeInfoDb::seqlengths(get(gen))[chrs]
    return(sizes)
}

get_enzyme_cutsites_v2 <- function(sig, gen = "Pfalciparum", gen_ver = "v24", chrs = NULL) {
    genome <- paste("BSgenome.", gen, '.PlasmoDB.', gen_ver, sep = "")
    library(genome, character.only = TRUE)
    if (is.null(chrs)) {
        chrs <- get_chrs(gen, gen_ver)
        chrs <- chrs[!chrs %in% c("PFC10_API_IRAB", "M76611")]
    }
    genome.chromSizes <- get_chr_sizes_v2(gen, gen_ver, chrs)
    while (any(grepl("N", sig))) {
        if (length(sig) <= 1) {
            sig <- unique(c(sub("N", "A", sig), sub("N", "C", sig), sub("N", "G", sig), sub("N", "T", sig)))
        } else {
            rng <- length(sig)
            sig_new <- NULL
            for (i in seq(rng)) {
                if (grepl("N", sig[i])) {
                sig_new_i <- unique(c(sub("N", "A", sig[i]), sub("N", "C", sig[i]), 
                                      sub("N", "G", sig[i]), sub("N", "T", sig[i])))
                sig_new <- c(sig_new, sig_new_i)
                } else {
                sig_new <- unique(c(sig_new, sig[i]))
                }
            }
            sig <- unique(sig_new)
        }
    }
    if (length(sig) <= 1) {
        enzymeCuts <- lapply(chrs, function(x) {
            range_chr <- IRanges::ranges(Biostrings::matchPattern(sig, get(gen)[[x]]))
            if (length(IRanges::start(range_chr)) <= 0) {
                warning(paste0("chromosome ", x, " does not have any cutsites with pattern ", sig, 
                               " .Using full chromosome instead"))
                range_chr <- IRanges::IRanges(start = 1, end = genome.chromSizes[x])
            }
            sort(GenomicRanges::GRanges(seqnames = x, ranges = range_chr))
        })
    } else {
        patternobj <- Biostrings::DNAStringSet(sig)
        enzymeCuts <- lapply(chrs, function(x) {
            range_chr <- IRanges::ranges(Biostrings::extractAllMatches(get(gen)[[x]], 
                                                                       Biostrings::matchPDict(patternobj, 
                                                                                              get(gen)[[x]])))
            if (length(IRanges::start(range_chr)) <= 0) {
                warning(paste0("chromosome ", x, " does not have any cutsites with pattern set ", 
                               paste(sig, collapse = ","), " .Using full chromosome instead"))
                range_chr <- IRanges::IRanges(start = 1, end = genome.chromSizes[x])
            }
            sort(GenomicRanges::GRanges(seqnames = x, ranges = range_chr))
        })
    }
    enzymeCuts <- suppressWarnings(sort(unlist(methods::as(enzymeCuts, "GRangesList"))))
    return(enzymeCuts)
}

construct_features_chr_v2 <- function(chrom, gen = "Pfalciparum", gen_ver = "v24", sig = "GATC", 
                                      bin_type = "Bins-uniform", binsize = 10000, wg_file = NULL, 
                                      feature_type = "RE-based") {
    getGC <- function(regions) {
        genome <- paste("BSgenome.", gen, '.PlasmoDB.', gen_ver, sep = "")
        library(genome, character.only = TRUE)
        gc <- rowSums(Biostrings::alphabetFrequency(BSgenome::Views(BiocGenerics::get(genome), regions), 
                                                    as.prob = TRUE)[, seq(2,3,1)])
        return(gc)
        }
    getMap <- function(regions, data) {
        hits <- GenomicRanges::findOverlaps(data, regions, type = "within")
        DT <- data.frame(queryHits = as.data.frame(hits)[[1]], subjectHits = as.data.frame(hits)[[2]])
        DT <- DT %>% dplyr::mutate(map = data$score[.data$queryHits], 
                                   len = GenomicRanges::end(data)[.data$queryHits] - 
                                       GenomicRanges::start(data)[.data$queryHits] + 1) %>% 
            dplyr::group_by(.data$subjectHits) %>% dplyr::summarize(avgmap = stats::weighted.mean(map, 
                                                                                                  w = len))
        map <- rep(0, length(regions))
        map[DT$subjectHits] <- DT$avgmap
        return(map)
        }
    genome.chromSize <- get_chr_sizes_v2(gen, gen_ver, chrom)[chrom]
    genome.chromGR <- GenomicRanges::GRanges(seqnames = chrom, 
                                             ranges = IRanges::IRanges(start = 1, end = genome.chromSize))
    if (feature_type=="RE-based"){
        msg <- paste0("Using ", paste(chrom, collapse = " "), "and cut patterns ", paste(sig, collapse = " "))
        message(msg)
        enzymeCutsites <- get_enzyme_cutsites_v2(sig, gen, gen_ver, chrom)
        RE_sites <- enzymeCutsites
        newends <- dplyr::lead(BiocGenerics::start(RE_sites)) - 1
        newends <- c(newends[-length(newends)], genome.chromSize)
        BiocGenerics::end(RE_sites) <- newends
        minstart <- min(BiocGenerics::start(RE_sites))
    if (minstart > 1) {
        RE_sites <- GenomeInfoDb::sortSeqlevels(c(GenomicRanges::GRanges(seqnames = chrom, 
                                                                         ranges = IRanges::IRanges(
                                                                             start = 1, 
                                                                             end = minstart)), RE_sites))
        RE_sites <- sort(RE_sites)
    }
    } else {
        if (bin_type!="Bins-uniform") stop("feature_type=RE-agnostic requires a fixed binsize.")
        msg <- paste0("Using ", paste(chrom, collapse = " "), " RE agnostic")
        message(msg)
    }
    if (!(is.null(wg_file))) {
        wgdata <- rtracklayer::import.bw(wg_file, which = genome.chromGR, as = "GRanges")
    } else {
        wgdata <- NULL
    }
    if (bin_type == "Bins-RE-sites") {
        endL <- GenomicRanges::GRanges(seqnames = as.character(GenomicRanges::seqnames(RE_sites)), 
                                       IRanges::IRanges(start = GenomicRanges::start(RE_sites), width = 200))
        BiocGenerics::end(endL) <- pmin(BiocGenerics::end(endL), genome.chromSize)
        BiocGenerics::start(endL) <- pmax(BiocGenerics::start(endL), 1)
        endR <- GenomicRanges::GRanges(seqnames = as.character(GenomicRanges::seqnames(RE_sites)), 
                                       IRanges::IRanges(end = GenomicRanges::end(RE_sites), width = 200))
        BiocGenerics::end(endR) <- pmin(BiocGenerics::end(endR), genome.chromSize)
        BiocGenerics::start(endR) <- pmax(BiocGenerics::start(endR), 1)
        gcL <- getGC(regions = endL)
        gcR <- getGC(regions = endR)
        gc <- (gcL + gcR) / 2
        RE_sites$gc <- gc
        if (!is.null(wgdata)) {
            endL <- GenomicRanges::GRanges(seqnames = as.character(GenomicRanges::seqnames(RE_sites)), 
                                           IRanges::IRanges(start = GenomicRanges::start(RE_sites), 
                                                            width = 500))
            BiocGenerics::end(endL) <- pmin(BiocGenerics::end(endL), genome.chromSize)
            BiocGenerics::start(endL) <- pmax(BiocGenerics::start(endL), 1)
            endR <- GenomicRanges::GRanges(seqnames = as.character(GenomicRanges::seqnames(RE_sites)), 
                                           IRanges::IRanges(end = GenomicRanges::end(RE_sites), width = 500))
            BiocGenerics::end(endR) <- pmin(BiocGenerics::end(endR), genome.chromSize)
            BiocGenerics::start(endR) <- pmax(BiocGenerics::start(endR), 1)
            mapL <- getMap(regions = endL, data = wgdata)
            mapR <- getMap(regions = endR, data = wgdata)
            map <- (mapL + mapR) / 2
            RE_sites$map <- map
        } else {
            RE_sites$map <- 0
        }
        bintolen <- as.data.frame(RE_sites, stringsAsFactors = FALSE)
        bintolen <- bintolen %>% dplyr::mutate(RE_id = dplyr::row_number())
        if (binsize > 1) {
            bintolen <- bintolen %>% dplyr::mutate(binNumb = floor((dplyr::row_number() - 1) / binsize) + 1)
            bintolen <- bintolen %>% dplyr::group_by(.data$seqnames, .data$binNumb) %>% 
                dplyr::summarize(start = min(start), end = max(end), gc = mean(gc), 
                                 map = mean(map), RE_id = min(.data$RE_id)) %>% 
                dplyr::mutate(width = .data$end - .data$start + 1) %>% 
                dplyr::rename(chr = "seqnames") %>% 
                dplyr::mutate(bins = paste(.data$chr, .data$start, .data$end, sep = "-"))
        } else {
            bintolen <- bintolen %>% dplyr::rename(chr = "seqnames") %>% 
                dplyr::mutate(bins = paste(.data$chr, .data$start, .data$end, sep = "-"))
        }
        bintolen <- bintolen %>% dplyr::select(.data$bins, .data$gc, .data$map, .data$width, .data$RE_id)
        bintolen <- as.data.frame(bintolen, stringsAsFactors = FALSE)
        rownames(bintolen) <- bintolen$bins
        return(bintolen)
    }
    if (bin_type == "Bins-uniform") {
        bins.chrom <- function(chrom, binsize) {
            cuts_start <- seq(1, genome.chromSize, by = binsize)
            cuts_end <- seq(binsize, genome.chromSize - (genome.chromSize%%binsize) + binsize, by = binsize)
            cuts_end <- c(cuts_end[-length(cuts_end)], genome.chromSize)
            bins <- GenomicRanges::GRanges(chrom, IRanges::IRanges(start = cuts_start, end = cuts_end))
            return(bins)
        }
        binsGR <- bins.chrom(chrom, binsize)
        names(binsGR) <- paste(as.character(GenomicRanges::seqnames(binsGR)), GenomicRanges::start(binsGR), 
                               GenomicRanges::end(binsGR), sep = "-")
        if (!is.null(wgdata)) {
            wgdata <- wgdata[GenomicRanges::seqnames(wgdata) == chrom]
        } else {
            wgdata <- NULL
        }
        if (feature_type=="RE-based"){
        medians <- (GenomicRanges::start(enzymeCutsites) + GenomicRanges::end(enzymeCutsites)) / 2
        FragmentendsL <- GenomicRanges::GRanges(seqnames = chrom, ranges = IRanges::restrict(
            IRanges::IRanges(end = medians - 1, width = 500), start = 1, end = genome.chromSize))
        FragmentendsR <- GenomicRanges::GRanges(seqnames = chrom, ranges = IRanges::restrict(
            IRanges::IRanges(start = medians, width = 500), start = 1, end = genome.chromSize))
        ends <- c(FragmentendsL, FragmentendsR)
        hits <- as.data.frame(GenomicRanges::findOverlaps(ends, binsGR, type = "within", select = "all"))
        LR <- data.frame(bins = names(binsGR[hits[[2]]]), start = GenomicRanges::start(ends[hits[[1]]]), 
                         end = GenomicRanges::end(ends[hits[[1]]]), stringsAsFactors = FALSE)
        } else {
        LR <- data.frame(bins=names(binsGR), start=GenomicRanges::start(binsGR), 
                         end=GenomicRanges::end(binsGR), stringsAsFactors = FALSE)
        }
        LRgr <- GenomicRanges::GRanges(chrom, IRanges::IRanges(start = LR$start, end = LR$end))
        if (!is.null(wgdata)) {
            LR$map <- getMap(LRgr, wgdata)
        } else {
            LR$map <- 0
        }
        map <- LR %>% dplyr::group_by(.data$bins) %>% dplyr::summarize(map = mean(map))
        if (feature_type=="RE-based"){
        ir <- IRanges::IRanges(LR$start, LR$end)
        LR$group <- as.data.frame(GenomicRanges::findOverlaps(ir, IRanges::reduce(ir)))[[2]]
        len <- LR %>% dplyr::group_by(.data$group, .data$bins) %>% 
            dplyr::summarize(start = min(start), end = max(end)) %>% 
            dplyr::mutate(width = .data$end - .data$start + 1) %>% 
            dplyr::group_by(.data$bins) %>% 
            dplyr::summarize(len = sum(.data$width))
        } else {
        len <- NULL
        }
        if (feature_type=="RE-based"){    
        FragmentendsL <- GenomicRanges::GRanges(seqnames = chrom, 
                                                ranges = IRanges::IRanges(end = medians - 1, width = 200))
        FragmentendsR <- GenomicRanges::GRanges(seqnames = chrom, 
                                                ranges = IRanges::IRanges(start = medians, width = 200))
        ends <- c(FragmentendsL, FragmentendsR)
        hits <- as.data.frame(GenomicRanges::findOverlaps(ends, binsGR, type = "within", select = "all"))
        LR2 <- data.frame(bins = names(binsGR[hits[[2]]]), start = GenomicRanges::start(ends[hits[[1]]]), 
                          end = GenomicRanges::end(ends[hits[[1]]]), stringsAsFactors = FALSE)
        } else {
        LR2 <- LR%>%dplyr::select(.data$bins, .data$start, .data$end)
        }
        LR2gr <- GenomicRanges::GRanges(chrom, IRanges::IRanges(start = LR2$start, end = LR2$end))
        LR2$gc <- getGC(LR2gr)
        gc <- LR2 %>% dplyr::group_by(.data$bins) %>% dplyr::summarize(gc = mean(gc))
        if (!is.null(len)){
        LR <- dplyr::left_join(dplyr::left_join(gc, map), len)
        } else {
        LR <- dplyr::left_join(gc, map)
        }
        bintolen <- as.data.frame(LR)
        if (sum(bintolen$map)==0){
            bintolen<-bintolen%>%dplyr::select(-.data$map)
        }
        rownames(bintolen) <- bintolen$bins
        allbins <- data.frame(bins = names(binsGR), stringsAsFactors = FALSE)
        bintolen <- suppressWarnings(dplyr::left_join(allbins, bintolen) %>% 
                                         tidyr::replace_na(list(gc = 0, map = 0, len = 0)))
        return(bintolen)
    }
}

construct_features_v2 <- function(output_path, gen = "Pfalciparum", gen_ver = "v24", sig = "GATC", 
                                  bin_type = "Bins-uniform", binsize = 10000, wg_file = NULL, 
                                  chrs = NULL, feature_type = "RE-based") {
    if (is.null(chrs)) {
        chrs <- get_chrs_v2(gen, gen_ver)
    }
    bintolen <- lapply(chrs,function(x) data.frame(stringsAsFactors = FALSE))
    names(bintolen) <- chrs
    for (chrom in names(bintolen)){
    bintolen[[chrom]] <- construct_features_chr_v2(chrom = chrom, gen = gen, gen_ver = gen_ver, sig = sig, 
                                                   bin_type = bin_type, binsize = binsize, wg_file = wg_file, 
                                                   feature_type = feature_type)
    }
    bintolen<-suppressWarnings(dplyr::bind_rows(bintolen))
    bintolenoutput <- path.expand(paste0(output_path, "_bintolen.txt.gz"))
    bintolenoutputdir <- gsub("/[^/]+$", "", bintolenoutput)
    if (bintolenoutputdir == bintolenoutput){
        bintolenoutputdir <- gsub("\\[^\\]+$", "", bintolenoutput)
    }
    if (bintolenoutputdir == bintolenoutput){
        bintolenoutputdir <- gsub("\\\\[^\\\\]+$", "", bintolenoutput)
    }
    if (!bintolenoutputdir == bintolenoutput&!dir.exists(bintolenoutputdir)){
        dir.create(bintolenoutputdir, showWarnings = FALSE, recursive = TRUE, mode = "0777")
    }
    data.table::fwrite(bintolen, bintolenoutput, row.names = FALSE, quote = FALSE, sep = "\t")
    return(bintolenoutput)
}

generate_bintolen_gi_list_v2 <- function(bintolen_path, chrs = NULL, Dthreshold = 2e+06, binned = TRUE, 
                                         binsize = NULL, gen = "Pfalciparum", gen_ver = "v24") {
    input.file.read <- function(filepath) {
        if (grepl("\\.txt.gz$", filepath) | grepl("\\.txt$", filepath)) {
            return(data.table::fread(filepath))
        } else if (grepl("\\.rds$", filepath)) {
            return(readRDS(filepath))
        } else {
            stop("Can only read in paths ending with .txt,.txt.gz, or .rds")
        }
    }
    if (binned & is.null(binsize) & !is.null(bintolen_path)) {
        bintolen <- input.file.read(bintolen_path)
        bintolen <- bintolen %>% tidyr::separate(.data$bins, c("chr", "start", "end"), "-") %>% 
            dplyr::mutate(start = floor(as.numeric(.data$start) / 1000) * 1000, end = as.numeric(.data$end))
        binsize <- round(mean(abs(bintolen$end - bintolen$start)) / 1000) * 1000
    } else if (binned & is.null(binsize) & is.null(bintolen_path)) {
        stop("If binned, need to specify at least one of binsize and bintolen_path")
    }
    if (!is.null(bintolen_path) & !exists("bintolen")) {
        bintolen <- input.file.read(bintolen_path)
        if (binned) {
            bintolen <- bintolen %>% tidyr::separate(.data$bins, c("chr", "start", "end"), "-") %>% 
                dplyr::mutate(start = floor(as.numeric(.data$start)/1000) * 1000, end = as.numeric(.data$end))
        } else {
            bintolen <- bintolen %>% tidyr::separate(.data$bins, c("chr", "start", "end"), "-") %>% 
                dplyr::mutate(start = as.numeric(.data$start), end = as.numeric(.data$end))
        }
    }
    if (is.null(chrs)) {
        chrs <- sort(unique(bintolen$chr))
    }
    if (binned) {
        gi_list <- generate_binned_gi_list_v2(binsize, chrs, Dthreshold, gen, gen_ver)
    } else {
        gi_list <- generate_df_gi_list(bintolen, chrs, Dthreshold, gen, gen_ver)
    }
    gi_list <- add_1D_features(gi_list, bintolen, chrs)
    return(gi_list)
}

generate_binned_gi_list_v2 <- function(binsize, chrs = NULL, Dthreshold = 2e+07, 
                                       gen = "Pfalciparum", gen_ver = "v24") {
    gi_list <- list()
    if (is.null(chrs)) {
        chrs <- get_chrs_v2(gen, gen_ver)
    } else {
        chrs <- chrs[chrs %in% get_chrs_v2(gen, gen_ver)]
    }
    for (chrom in chrs) {
        seqlen <- get_chr_sizes_v2(gen, gen_ver)[chrom]
        numbins <- ceiling(seqlen/binsize)
        maxbins <- min(round(Dthreshold/binsize), numbins)
        all.regions <- GenomicRanges::GRanges(
            chrom, IRanges::IRanges(start = seq(0, (numbins - 1), 1) * binsize, 
                                    end = pmin(seq(1, numbins, 1) * binsize, seqlen)))
        index1 <- unlist(lapply(seq(1, numbins, 1), function(x) rep(x, min(maxbins + 1, numbins - x + 1))))
        index2 <- unlist(lapply(seq(1, numbins, 1), function(x) seq(x, min(x + maxbins, numbins), 1)))
        gi_list[[chrom]] <- InteractionSet::GInteractions(index1, index2, all.regions)
        mcols(gi_list[[chrom]])$D <- InteractionSet::pairdist(gi_list[[chrom]])
        gi_list[[chrom]] <- gi_list[[chrom]][mcols(gi_list[[chrom]])$D <= Dthreshold]
    }
    return(gi_list)
}

add_hicpro_matrix_counts_v2 <- function(gi_list, absfile_path, matrixfile_path, chrs = NULL, add_inter = FALSE) {
    gi_list_validate(gi_list)
    if (is.null(chrs)) {
        chrs <- sort(names(gi_list))
    } else {
        chrs <- chrs[chrs %in% names(gi_list)]
        if (length(chrs) == 0) {
            stop("None of the chromosomes specified exists in the gi.list object")
        }
    }
    Dthreshold <- gi_list_Dthreshold.detect(gi_list)
    absI <- data.table::fread(absfile_path, sep = "\t", header = FALSE) %>% dplyr::rename(chrI = "V1", startI = "V2", end = "V3", 
                                                                                          first = "V4") %>% dplyr::select(-.data$end)
    absJ <- data.table::fread(absfile_path, sep = "\t", header = FALSE) %>% dplyr::rename(chrJ = "V1", startJ = "V2", end = "V3", 
                                                                                          second = "V4") %>% dplyr::select(-.data$end)
    count_matrix <- data.table::fread(matrixfile_path, sep = "\t", stringsAsFactors = FALSE) %>% dplyr::rename(first = "V1", 
                                                                                                               second = "V2", counts = "V3")
    count_matrix <- dplyr::left_join((dplyr::left_join(as.data.frame(count_matrix), absI) %>% dplyr::filter(!is.na(.data$chrI))), 
                                     absJ) %>% dplyr::filter(!is.na(.data$chrJ)) %>% dplyr::filter(.data$chrI == .data$chrJ & abs(.data$startJ - .data$startI) <= 
                                                                                                       Dthreshold) %>% dplyr::select(.data$chrI, .data$startI, .data$startJ, .data$counts) %>% dplyr::rename(chr = "chrI")
    rm(absI, absJ)
    for (chrom in chrs) {
        count_matrix_chr <- count_matrix %>% dplyr::filter(.data$chr == chrom)
        if (nrow(count_matrix_chr)==0){
            msg<-paste0(chrom, "does not have any counts in this file. Dropping from gi_list.")
            warning(msg)
            gi_list[[chrom]]<-NULL
            next
        }
        gi_list[[chrom]] <- add_2D_features(gi_list[[chrom]], count_matrix_chr)
        msg<-paste0("Chromosome ",chrom," intrachromosomal counts processed.")
        message(msg)
    }
    rm(count_matrix_chr, count_matrix)
    # intercounts--hicpro files
    if (add_inter) {
        absI <- data.table::fread(absfile_path, sep = "\t", header = FALSE) %>% dplyr::rename(chrI = "V1", startI = "V2", 
                                                                                              end = "V3", first = "V4") %>% dplyr::select(-.data$end)
        absJ <- data.table::fread(absfile_path, sep = "\t", header = FALSE) %>% dplyr::rename(chrJ = "V1", startJ = "V2", 
                                                                                              end = "V3", second = "V4") %>% dplyr::select(-.data$end)
        count_matrix <- data.table::fread(matrixfile_path, sep = "\t", stringsAsFactors = FALSE) %>% dplyr::rename(first = "V1", 
                                                                                                                   second = "V2", counts = "V3")
        count_matrix <- dplyr::left_join((dplyr::left_join(as.data.frame(count_matrix), absI) %>% dplyr::filter(!is.na(.data$chrI))), 
                                         absJ) %>% dplyr::filter(!is.na(.data$chrJ)) %>% dplyr::filter(!.data$chrI == .data$chrJ)
        count_matrix <- dplyr::bind_rows(count_matrix %>% dplyr::group_by(.data$chrI, .data$startI) %>% dplyr::summarize(counts = sum(.data$counts)) %>% 
                                             dplyr::rename(chr = "chrI", start = "startI"), count_matrix %>% dplyr::group_by(.data$chrJ, .data$startJ) %>% 
                                             dplyr::summarize(counts = sum(.data$counts)) %>% dplyr::rename(chr = "chrJ", start = "startJ"))
        count_matrix <- count_matrix %>% dplyr::group_by(.data$chr, .data$start) %>% dplyr::summarize_all(sum) %>% dplyr::ungroup() %>% 
            dplyr::rename(inter = "counts")
        rm(absI, absJ)
        gi_list <- add_1D_features(gi_list, count_matrix, chrs, agg = sum)
        rm(count_matrix)
    }
    return(gi_list)
}