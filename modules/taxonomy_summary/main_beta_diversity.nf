
process BETA_DIVERSITY {

    tag 'Bray-Curtis PCoA'

    publishDir "${projectDir}/results/taxonomy_summary", mode: 'copy'

    cpus 1

    input:
    path taxonomy_table

    output:
    path 'bray_curtis_distance.tsv', emit: bray_curtis
    path 'pcoa_coordinates.tsv', emit: coordinates
    path 'pcoa_bray_curtis.png', emit: pcoa_plot

    script:
    """
    Rscript - "${taxonomy_table}" <<'RSCRIPT'
    args <- commandArgs(trailingOnly = TRUE)
    input_file <- args[1]

    if (!requireNamespace("vegan", quietly = TRUE)) {
        stop("R package 'vegan' is required.")
    }
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("R package 'ggplot2' is required.")
    }

    # Read genus-level abundance table
    tax <- read.delim(
        input_file,
        check.names = FALSE,
        stringsAsFactors = FALSE
    )

    if (!"Genus" %in% colnames(tax)) {
        stop("Input table must contain a 'Genus' column.")
    }

    rownames(tax) <- tax\$Genus
    tax\$Genus <- NULL

    # Convert to numeric and transpose:
    # rows = samples, columns = genera
    tax[] <- lapply(tax, function(x) as.numeric(as.character(x)))
    community <- t(as.matrix(tax))
    storage.mode(community) <- "numeric"

    if (anyNA(community) || any(community < 0)) {
        stop("Abundances contain missing, non-numeric, or negative values.")
    }
    if (any(rowSums(community) == 0)) {
        stop("At least one sample has zero total abundance.")
    }

    # Normalize each sample to relative abundance
    community <- community / rowSums(community)

    # Bray-Curtis distances
    bray <- vegan::vegdist(community, method = "bray")
    write.table(
        as.matrix(bray),
        "bray_curtis_distance.tsv",
        sep = "\\t",
        quote = FALSE,
        col.names = NA
    )

    # PCoA
    pcoa <- cmdscale(bray, k = 2, eig = TRUE)
    positive_eig <- pcoa\$eig[pcoa\$eig > 0]
    variance <- pcoa\$eig[1:2] / sum(positive_eig) * 100

    coords <- data.frame(
        Sample = rownames(pcoa\$points),
        PCoA1 = pcoa\$points[, 1],
        PCoA2 = pcoa\$points[, 2],
        stringsAsFactors = FALSE
    )

    metadata <- data.frame(
        Sample = c("FBI016_W0", "FBI016_W16",
                   "FBI029_W0", "FBI029_W16"),
        Participant = c("FBI016", "FBI016", "FBI029", "FBI029"),
        Timepoint = c("W0", "W16", "W0", "W16")
    )

    coords <- merge(coords, metadata, by = "Sample", sort = FALSE)
    coords <- coords[match(metadata\$Sample, coords\$Sample), ]

    write.table(
        coords,
        "pcoa_coordinates.tsv",
        sep = "\\t",
        quote = FALSE,
        row.names = FALSE
    )

    # Plot samples and connect paired timepoints
    library(ggplot2)

    p <- ggplot(coords, aes(PCoA1, PCoA2)) +
        geom_path(
            aes(group = Participant),
            color = "grey60",
            linewidth = 0.7
        ) +
        geom_point(
            aes(color = Timepoint, shape = Participant),
            size = 4
        ) +
        geom_text(
            aes(label = Sample),
            vjust = -0.8,
            size = 3,
            show.legend = FALSE
        ) +
        labs(
            title = "Genus-level beta diversity",
            subtitle = "Bray-Curtis dissimilarity; exploratory paired comparison",
            x = sprintf("PCoA1 (%.1f%%)", variance[1]),
            y = sprintf("PCoA2 (%.1f%%)", variance[2]),
            color = "Timepoint",
            shape = "Participant"
        ) +
        theme_classic()

    ggsave(
        "pcoa_bray_curtis.png",
        p,
        width = 8,
        height = 6,
        dpi = 300
    )
    RSCRIPT
    """

    stub:
    """
    touch bray_curtis_distance.tsv
    touch pcoa_coordinates.tsv
    touch pcoa_bray_curtis.png
    """
}


