pdf("four_GBM_LGG_ionchannels_of_interest.pdf")

get_km_plot(get_ensg_pcg("CATSPER1"), "GBM")
get_km_plot(get_ensg_pcg("CATSPER1"), "LGG")

get_km_plot(get_ensg_pcg("SCN9A"), "GBM")
get_km_plot(get_ensg_pcg("SCN9A"), "LGG")

get_km_plot(get_ensg_pcg("AQP9"), "GBM")
get_km_plot(get_ensg_pcg("AQP9"), "LGG")

get_km_plot(get_ensg_pcg("KCNN4"), "GBM")
get_km_plot(get_ensg_pcg("KCNN4"), "LGG")

dev.off()