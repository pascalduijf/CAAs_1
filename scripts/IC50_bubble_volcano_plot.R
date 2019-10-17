c=c("#000000", "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
    "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
    "#5A0007", "#809693", "#0089A3", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80")
df2 %>%
    arrange(desc(N_FEATURE_pos)) %>%
    
    ggplot(aes(x=FEATURE_IC50_effect_size2, y=-log10(FEATURE_pval), size=N_FEATURE_pos, color=DRUG_NAME))  +
    geom_point(alpha=0.4) +
    scale_size(range = c(2,16), name="No. of altered cell lines") + geom_vline(xintercept = 0) + scale_y_continuous(limits = c(2,9),breaks = c(3,5,9))  +  scale_color_manual(values = c,name="TCGA Cell Line Code")+ theme_minimal() + labs(y="-log10(Pvalue)",x="IC50 Effect Size") + scale_x_continuous(limits = c(-11,11)) 
