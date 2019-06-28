dataFiles <- list(
    TCGA=list(
        COAD="data/coad_b2c.rda",
        LIHC="data/lihc_b2c.rda",
        PRAD="data/prad_b2c.rda",
        LUAD="data/luad_b2c.rda",
        BRCA="data/brca_b2c.rda"
    ),
    PHLab=list(
        ASCL2_KD_3I4="data/ascl2_kd_3i4_b2c.rda",
        WiNTRLINC1_KD_76="data/wintrlinc1_kd_76_b2c.rda",
        lincIGSF9_KD="data/lincigsf9_kd_b2c.rda"
    ),
    ITLab=list(
        Mouse_Liver_Development="data/mouse_liver_development_b2c.rda"
    )
)

#config=read.delim("./config/metadata.txt")
#rownames(config) <- as.character(config$sample_id)
#buildTrackList(
#   config=config,
#   annoPath="/media/raid/tracks/hybridsuite/reference",
#   appBase="/media/raid/software/bigseqcvis",
#   urlBase="http://epigenomics.fleming.gr/tracks"
#)
