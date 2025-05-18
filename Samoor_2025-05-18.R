# loading -----------------------------------------------------------------
export <- F # replace to TRUE in the final script exectution
library(indicspecies)
library(tidyverse)
library(dendextend)
library(parallel)
# cl <- makeCluster(detectCores()-1)
theme_set(theme_bw() + theme(legend.position = "bottom"))

long <- readxl::read_excel("Caspian data_2024-07-21_SA.xlsx", sheet = "main") %>% 
    select(O, sp, contains("Sm", ignore.case = FALSE)) %>% 
    # mutate(apply(.[,3:82], 1, function(a){length(unique(a))}) == 0}), .after = sp) %>% pull(3)
    pivot_longer(
        names_to = "id",
        values_to = "v", 
        -c("O", "sp"), 
        values_transform = as.character) %>%
    separate(v, into = c("ad", "jv"), sep = "\\+", convert = TRUE) %>% # 
    mutate(
        id = substr(id, 3, 7),
        abu = ad+case_when(is.na(jv) ~ 0, TRUE ~ jv), 
        .keep = "unused") %>% 
    group_by(sp) %>% 
    mutate(total = sum(abu)) %>% 
    ungroup %>% 
    filter(total > 0) %>% 
    select(-total)

orwi <- long %>% 
    filter(O == "Oribatida") %>% 
    select(-O) %>%
    pivot_wider(names_from = id, values_from = abu, values_fn = sum, values_fill = 0)
orws <- long %>% 
    filter(O == "Oribatida") %>% 
    mutate(seria = substr(id, 1, nchar(id)-1)) %>% 
    select(-O, -id) %>% 
    pivot_wider(names_from = seria, values_from = abu, values_fn = sum, values_fill = 0)
mswi <- long %>% 
    filter(O == "Mesostigmata") %>% 
    select(-O) %>%
    pivot_wider(names_from = id, values_from = abu, values_fn = sum, values_fill = 0)
msws <- long %>% 
    filter(O == "Mesostigmata") %>% 
    mutate(seria = substr(id, 1, nchar(id)-1)) %>% 
    select(-O, -id) %>%
    pivot_wider(names_from = seria, values_from = abu, values_fn = sum, values_fill = 0)


labs <- readxl::read_excel("Caspian data_2024-07-21_SA.xlsx", sheet = "samples") %>% 
    filter(distr == "Samoor") %>% 
    transmute(id = substr(id, 3, 7), 
              seria = substr(id, 1, nchar(id)-1), 
              plants.d, 
              coast = case_when(coast == "sandy dunes" ~ "dunes", TRUE ~ coast),
              N = `N, %`, 
              C = `C, %`, 
              CN = `C/N`,
              RH) %>% 
    separate(plants.d, sep = " ", into = "plants.d", extra = "drop")
taxa <- readxl::read_excel("Caspian data_2024-07-21_SA.xlsx", sheet = "taxa") %>% 
    mutate(
        order = factor(order, levels = c("Mesostigmata", "Oribatida", "Astigmata")),
        species = case_when(
            is.na(author) ~ paste0(genus, " ", species),
            is.na(brackets) ~ paste0(genus, " ", species, " ", author, ", ", year), 
            brackets == 0 ~ paste0(genus, " ", species, " ", author, ", ", year), 
            TRUE ~ paste0(genus, " ", species, " (", author, ", ", year, ")"), 
        ), 
        .keep = "unused") %>% 
    filter(sp %in% long$sp) %>% 
    select(-cf) %>% 
    arrange(order, species)

tables <- list()

# basic counts -------------------------------------------------------------
long %>% 
    filter(O != "total") %>% 
    select(O, sp) %>% 
    split(.$O) %>% 
    map(~nrow(distinct(.x)))

taxa %>% 
    split(.$order) %>% 
    map(~ .x %>% 
            select(order, family) %>% 
            distinct)

taxa %>% 
    split(.$order) %>% 
    map(~.x %>% 
            pull(new_record) %>% 
            toupper) %>% 
    lapply(function(a){a[is.na(a)] <- "NO"; table(a)})

taxa %>% 
    split(.$order) %>% 
    map(~.x %>% 
            count(range) %>% 
            filter(range != "unknown") %>% 
            mutate(p = round(n/sum(n)*100, 1)))

tables$tab.taxa_counts1 <- taxa %>% 
    filter(order != "Astigmata") %>% 
    group_by(order, family) %>% 
    count() %>% 
    ungroup %>% 
    arrange(order, desc(n))

tables$tab.taxa_counts2 <- taxa %>% 
    filter(order != "Astigmata") %>% 
    separate(sp, into = c("genus", "sp"), sep = " ", extra = "merge") %>% 
    group_by(order, family, genus) %>% 
    count() %>% 
    ungroup %>% 
    arrange(order, desc(n))

# + Рис. 2, таб. 1. Корреляция обилия групп друг с другом и RH%, C%. ---------------
# BY SAMPLES
# Сделано по пробам (не по сериям). Это корректнее, т.к. у каждой серии 
# не только свой набор видов, но и свои значения факторов среды
cor.data1 <- long %>% 
    filter(O == "total", str_detect(sp, "adult|juven", negate = TRUE)) %>% 
    transmute(sp = str_replace(sp, "_total", ""), 
              id, 
              abu) %>% 
    pivot_wider(names_from = sp, values_from = abu) %>% 
    left_join(select(labs, id, N:RH), by = "id") %>% 
    select(1, 4, 6, 5, 3, 2, 8, 7, 9:10)

cor.val1 <- cor(cor.data1[,2:ncol(cor.data1)], method = "spearman")
cor.pval1 <- expand_grid(v1 = 2:ncol(cor.data1), v2 = 2:ncol(cor.data1)) %>% 
    split(1:nrow(.)) %>% 
    lapply(function(a){
        p = cor.test(as_vector(cor.data1[,a$v1]), as_vector(cor.data1[,a$v2]), 
                     method = "spearman")
        data.frame(p = p$p.value, 
                   v1 = colnames(cor.data1)[a$v1], 
                   v2 = colnames(cor.data1)[a$v2])
    }) %>% 
    map_df(tibble) %>% 
    mutate(p = p.adjust(p, method = "BY")) %>% 
    pivot_wider(names_from = v2, values_from = p) %>% 
    column_to_rownames("v1") %>% 
    as.matrix()

if(export){
    pdf(
        paste0("export/Fig 2. Correlation_by.samples ", Sys.Date(), ".pdf"), 
        width = 6, height = 4
    )
}
corrplot::corrplot(
    corr = cor.val1[1:5,], 
    p.mat = cor.pval1[1:5,], 
    type="upper", 
    order = "original",
    diag = FALSE,
    col =  corrplot::COL2('RdYlBu', 10)[10:1], 
    sig.level = 0.05)
if(export){
    dev.off()    
}

tables$tab1_cor_id.samples <- paste(
    "ρ=", round(cor.val1, 2), "; p=", round(cor.pval1, 2), sep = "") %>% 
    matrix(ncol = 9, byrow = TRUE) %>% 
    `colnames<-`(colnames(cor.val1)) %>% 
    `rownames<-`(rownames(cor.val1)) %>% 
    as.data.frame() %>% 
    rownames_to_column("taxa") 
tables$tab1_cor_id.samples

# + Табл. 2/3. Таксономический состав панцирных и мезостигматических -----------
tables$tab.2_all_species <- long %>% 
    filter(O != "total", abu > 0) %>% 
    left_join(labs, by = "id") %>% 
    group_by(O, sp, coast) %>% 
    summarise(abu = sum(abu), .groups = "drop") %>% 
    left_join(taxa, by = "sp") %>% 
    pivot_wider(names_from = coast, values_from = abu, values_fill = 0) %>% 
    select(sp, Order = O, Family = family, Species = species, 
           habitat = habit, distribution:pebbly) %>% 
    mutate(new_record = case_when(
        new_record == "science" ~ "***",
        new_record == "Russia" ~ "**",
        new_record == "Dagestan" ~ "*", 
        TRUE ~ ""
    )) %>% 
    rownames_to_column("No.") 

tables$tab.3_dom.species <- long %>% 
    group_by(O, sp) %>% 
    mutate(
        dom = max(abu), 
        sp = case_when(
            O == "Mesostigmata" & dom < 5 ~ "other Mesostigmata", 
            O == "Oribatida" & dom < 5 ~ "other Oribatida", 
            TRUE ~ sp)) %>% 
    group_by(O, sp, id) %>% 
    summarise(abu = sum(abu), .groups = "drop") %>% 
    left_join(select(labs, id, plants.d, coast), by = "id") %>% 
    filter(abu > 0, O == "Mesostigmata" | O == "Oribatida") %>% 
    mutate(seria = substr(id, 1, nchar(id)-1), .keep = "unused") %>% 
    group_by(seria) %>% 
    mutate(
        plants.d = paste0(sort(unique(plants.d)), collapse = "/"), 
        coast =    paste0(sort(unique(coast)),    collapse = "/")
    ) %>% 
    unite(seria, seria, coast, plants.d, sep = " - ") %>% 
    arrange(seria, O, sp) %>% 
    pivot_wider(names_from = seria, values_from = abu, values_fn = sum) %>% 
    arrange(O, sp)

# + Табл. 4. Виды-специалисты (индикаторы) ----------------------------------
tmp <- data.frame(id = colnames(select(orwi, -sp, -Se1:-Se5))) %>% 
    left_join(distinct(select(labs, id, coast)), by = "id") %>% 
    pull(coast)

set.seed(3); iv.i <- orwi %>% 
    rbind(mswi) %>%
    select(-Se1:-Se5) %>% 
    mutate(total = apply(.[,-1], 1, sum), .before = 2) %>% 
    filter(total >= 5) %>% 
    select(-total) %>% 
    column_to_rownames("sp") %>% 
    t %>%
    indicspecies::multipatt(
        tmp,
        # pull(arrange(distinct(select(labs, id, coast)), id) , coast), 
        control = how(nperm=999),
        max.order = 4,
        func = "indval",
        duleg = FALSE
    )
iv.i <- iv.i$str %>% 
    as.data.frame() %>% 
    rownames_to_column("sp") %>% 
    as_tibble() %>% 
    pivot_longer(names_to = "biotop", values_to = "iv", -sp) %>% 
    group_by(sp) %>% 
    filter(iv == max(iv)) %>% 
    ungroup() %>% 
    left_join(select(rownames_to_column(iv.i$sign, "sp"), sp, p.value), 
              by = "sp") %>% 
    mutate(`iv, %` = round(iv*100), 
           p.value = round(p.value, 4),
           sign = case_when(p.value <= 0.001 ~ "***", 
                            p.value <= 0.01 ~ "**", 
                            p.value <= 0.05 ~ "*", 
                            TRUE ~ ""), 
           .keep = "unused", 
           .after = 2) %>% 
    left_join(select(taxa, sp, order), by = "sp") %>% 
    select(order, sp:`iv, %`, p.value, sign) 


tables$tab.4_indicator_species <- filter(iv.i, p.value <= 0.05)

if(export){
tables$tab.4_indicator_species %>% 
    DT::datatable(., 
        filter = 'top', 
        extensions = c('FixedColumns',"FixedHeader"),
        options = list(
            scrollX = TRUE, 
            paging=FALSE,
            fixedHeader=TRUE)) %>% 
    DT::formatStyle('sp', fontStyle = list(fontStyle = 'italic')) %>% 
    htmlwidgets::saveWidget(
        file = "indicator_species_tmp.html", 
        selfcontained = TRUE)
    
zip(
    zipfile = paste0("export/indicator_species_", Sys.Date()),
    files = "indicator_species_tmp.html"
); file.remove("indicator_species_tmp.html")
} else {
    tables$tab.4_indicator_species    
}

# + Рис. 4. Обилие микроартропод по типам берега и парцеллам  ---------------
abundance <- long %>% 
    filter(O == "total", str_detect(sp, "adult|juven", negate = TRUE)) %>% 
    left_join(labs, by = "id") %>%
    mutate(
        sp = str_replace_all(sp, "_total", ""),
        sp = factor(sp, levels = c("Collembola", "Astigmata", 
                                   "Mesostigmata", "Oribatida", "Prostigmata")))

xaxis <- c("SdJj", "SdEq", "SdTu", "SdEc", "SdTa", "SdJm", "SdFn", 
           "PbAe", "PbDe", "PbPo", "PbTu", "PbTl", "DuCJ", "DuCS", "RsFd", "Se") 

p4a <- abundance %>% 
    group_by(seria, coast, sp) %>% 
    summarise(abu = sum(abu), 
              plants.e = paste0(unique(plants.d), collapse = " / "), 
              .groups = "drop") %>% 
    ggplot(aes(x = seria, y = abu, fill = sp)) + 
    geom_col(width = 0.68) + # position = "dodge"
    geom_text(mapping = aes(y = 3000, label = plants.e), angle = 90, color = "black", fontface = "italic") +
    geom_text(mapping = aes(y = 5200, label = coast), angle = 90, color = "black") +
    scale_fill_manual(values = c("#ABA300", "#C77CFF", "#00B8E7", "#F8766D", "#00c19A"), drop = TRUE) + 
    scale_x_discrete(limits = xaxis) +
    labs(y = NULL, x = NULL, fill = NULL) + 
    theme(axis.text.x = element_text(angle = 90, hjust = -0.5), 
          legend.position = "top")
p4b <- abundance %>% 
    filter(sp %in% c("Oribatida", "Mesostigmata")) %>% 
    group_by(seria, coast, sp) %>%  
    summarise(abu = sum(abu), 
              plants.e = paste0(unique(plants.d), collapse = " / "), 
              .groups = "drop") %>% 
    ggplot(aes(x = seria, y = abu, fill = sp)) + 
    geom_col(width = 0.68) + # position = "dodge"
    geom_text(mapping = aes(y = 1200, label = plants.e), angle = 90, color = "black", fontface = "italic") +
    geom_text(mapping = aes(y = 2000, label = coast), angle = 90, color = "black") +
    scale_x_discrete(limits = xaxis) +
    scale_y_reverse() + 
    scale_fill_manual(values = c("#00B8E7", "#F8766D")) + 
    guides(fill="none") +
    labs(y = NULL, x = NULL) + 
    theme(axis.text.x = element_text(angle = 90, hjust = -0.5))
p4 <- gridExtra::grid.arrange(p4a, p4b, ncol = 1, 
                              left = "Total abundance, individuals")
if(export){
    ggsave(
        paste0("export/Fig 4. General abundances ", Sys.Date(), ".pdf"), 
        plot = p4, 
        width = 297/25, height = 210/25)
}



# + Рис. 5. Куммуляты по типам берега ---------------------------------
s <- Sys.time()
rar <- long %>% 
    filter(O == "Oribatida" | O == "Mesostigmata") %>% 
    left_join(select(labs, coast, id), by = "id") %>% 
    group_by(O, sp, coast) %>% 
    summarise(abu = sum(abu), .groups = "drop") %>% 
    filter(abu > 0) %>% 
    unite("type", O, coast) %>% 
    split(.$type) %>% 
    map(~mutate(.x, s = case_when(str_detect(type, "Meso") ~ 800, TRUE ~ 2000))) %>% 
    lapply(function(a){
        list(abu = sort(a$abu), 
             s = max(a$s))
    }) %>% 
    mclapply(mc.cores = parallel::detectCores(), FUN = function(
        a, 
        nb = switch(as.character(export), 
                    "TRUE" = 999, "FALSE" = 9)
        ){
        iNEXT::iNEXT( 
            a$abu,
            q=0, 
            nboot = nb,
            datatype="abundance", 
            se = TRUE,
            size = seq(0, a$s, by = 5)) |>
            purrr::pluck("iNextEst", "size_based") |>
            dplyr::select(m, Method, qD, qD.LCL, qD.UCL)
        
    }) %>% 
    map_dfr(rbind, .id = "type") %>% 
    distinct() %>% 
    separate(type, into = c("taxa", "coast"), sep = "_") %>% 
    as_tibble()
Sys.time() - s # 36 mins
rar %>% 
    group_by(taxa, coast) %>% 
    count()

p5o <- rar %>% # основа
    ggplot(aes(x = m, y = qD, group = coast, fill = coast, color = coast)) + 
    labs(x = "individuals", y = "Number of species") + 
    facet_wrap(~taxa, scales = "free")
p5b <- p5o + # подложка из доверительных областей
    geom_ribbon(aes(ymin = qD.LCL, ymax = qD.UCL), alpha = 0.3, color = "transparent")
p5o <- p5o +
    geom_line(data = filter(rar, Method != "Extrapolation")) +
    geom_line(data = filter(rar, Method == "Extrapolation"), linetype = "dashed") + 
    geom_point(data = filter(rar, Method == "Observed"), 
               shape = 15, size = 2)
p5o # without confidence areas
if(export){
    ggsave(paste0("export/Fig 5. Rarefication ", Sys.Date(), ".pdf"), 
           width = 297*0.6, height = 150*0.6, units = "mm")
}

p5b + # with confidence areas
    geom_line(data = filter(rar, Method != "Extrapolation")) +
    geom_line(data = filter(rar, Method == "Extrapolation"), linetype = "dashed") + 
    geom_point(data = filter(rar, Method == "Observed"), 
               shape = 15, size = 2)

# + Рис. 6. Кладограмма по фаунистическим спискам отдельных мероценозов --------
dis <- list()
# dissimilarity
dis$or.bin <- orwi %>% 
    column_to_rownames("sp") %>%
    select_if(function(a){sum(a)>0}) %>%
    t %>% 
    as.data.frame() %>% 
    vegan::vegdist(method = "jaccard", binary = TRUE)
dis$or.num <- orwi %>% 
    column_to_rownames("sp") %>% 
    select_if(function(a){sum(a)>0}) %>% 
    t %>% 
    as.data.frame() %>% 
    vegan::vegdist(method = "bray", binary = FALSE)
dis$ms.bin <- mswi %>% 
    column_to_rownames("sp") %>% 
    select_if(function(a){sum(a)>0}) %>% 
    t %>% 
    as.data.frame() %>% 
    vegan::vegdist(method = "jaccard", binary = TRUE)
dis$ms.num <- mswi %>% 
    column_to_rownames("sp") %>% 
    select_if(function(a){sum(a)>0}) %>% 
    t %>% 
    as.data.frame() %>% 
    vegan::vegdist(method = "bray", binary = FALSE)

# https://cran.r-project.org/web/packages/dendextend/vignettes/dendextend.html

dend <- lapply(dis, function(a){ 
    dend <- a %>% 
        hclust(method = 'ward.D2') %>% ### METHOD
        as.dendrogram()
    L <- labels(dend) %>% 
        tibble(id = .) %>% 
        left_join(labs, by = "id") %>% 
        mutate(l = factor(coast), 
               l = as.numeric(l)) 
    dend %>% 
        set("labels_col", L$l) %>% 
        set("labels_cex", 0.5)
})

# par(mfrow = c(2,2))
if(export){
    pdf(paste0("export/Fig 6. Dendrogramms ", Sys.Date(), ".pdf"), 
        width = 7, height = 7)
}
plot(dend[[1]],
     horiz = TRUE, 
     main = "Order = Oribatida\n Data = binary (Jaccard)\n Method = Ward")
# plot(dend[[2]],
#      horiz = TRUE, 
#      main = "Order = Oribatida\n Data = numeric (Bray-Curtis)\n Method = Ward")
plot(dend[[3]],
     horiz = TRUE, 
     main = "Order = Mesostigmata\n Data = binary (Jaccard)\n Method = Ward")
# plot(dend[[4]],
#      horiz = TRUE, 
#      main = "Order = Mesostigmata\n Data = numeric (Bray-Curtis)\n Method = Ward")
if(export){
    dev.off()
}

# + Рис. 7. Видовое богатство по сериям  ------------------------------------
# в среднем в пробе этой серии
div <- orwi %>% 
    select(-sp) %>% 
    as.list() %>% 
    lapply(function(a){data.frame(Oribatida = length(a[a>0]))}) %>% 
    map_df(rbind, .id = "id")
div <- mswi %>% 
    select(-sp) %>% 
    as.list() %>% 
    lapply(function(a){data.frame(Mesostigmata = length(a[a>0]))}) %>% 
    map_df(rbind, .id = "id") %>% 
    left_join(div, by = "id") %>% 
    inner_join(select(labs, -RH, -C, -N, -CN), by = "id") %>% 
    as_tibble()

p7 <- div %>% 
    select(-id) %>% 
    pivot_longer(names_to = "taxa", values_to = "nsp", -c("seria", "coast", "plants.d")) %>% 
    group_by(seria, taxa) %>% 
    summarise(
        coast = paste0(unique(coast), collapse = " / "), 
        plants.e = paste0(unique(plants.d), collapse = " / "), 
        nsp_mean = mean(nsp),
        nsp_sd   = sd(nsp), 
        .groups = "drop") %>% 
    mutate(ymax = nsp_mean + nsp_sd, 
           ymin = nsp_mean - nsp_sd, 
           ymin = case_when(ymin <= 0 ~ 0, TRUE ~ ymin)) %>% 
    ggplot(aes(x = seria, y = nsp_mean, 
               fill = taxa, 
               ymin = ymin, ymax = ymax)) + 
    geom_col(position = "dodge", width = 0.68) + 
    geom_errorbar(position = "dodge", color = "black", alpha = 0.5) + 
    geom_text(mapping = aes(y = 10, label = plants.e), angle = 90, color = "black") +
    geom_text(mapping = aes(y = 4, label = coast), angle = 90, color = "black") +
    scale_x_discrete(limits = xaxis) +
    scale_fill_manual(values = c("#00B8E7", "#F8766D")) + 
    theme(axis.text.x = element_text(angle = 90)) +
    labs(x = NULL, fill = NULL, #y = "Среднее количество видов в серии ± SD")
         y = "Average number of species in seria ± SD")
p7
if(export){
    ggsave(paste0("export/Fig 7. Diversity by series ", Sys.Date(), ".pdf"), 
           p7, width = 297/25, height = 210/25)
}

tables$tab.for.fig7 <- div %>% 
    select(-id) %>% 
    pivot_longer(names_to = "taxa", values_to = "nsp", -c("seria", "coast", "plants.d")) %>% 
    group_by(seria, taxa) %>% 
    summarise(
        coast = paste0(unique(coast), collapse = " / "), 
        plants.e = paste0(unique(plants.d), collapse = " / "), 
        nsp_mean = paste0(mean(nsp), " ± ", round(sd(nsp), 1)), 
        .groups = "drop") %>% 
    pivot_wider(names_from = taxa, values_from = nsp_mean)

# + Рис. 8. Ареалогический состав по сериям: качественный и количест --------
p8 <- rbind(mutate(mswi, Order = "Mesostigmata", .before = 1), 
            mutate(orwi, Order = "Oribatida", .before = 1)
            # mutate(rbind(mswi, orwi), .before = 1,
            #        Order = "Oribatida & Mesostigmata")
) %>% 
    pivot_longer(names_to = "id", values_to = "abu", -Order:-sp) %>% 
    left_join(taxa, by = "sp") %>% 
    filter(range != "unknown" | Order != "Oribatida & Mesostigmata", 
           abu > 0) %>% 
    transmute(
        seria = substr(id, 1, nchar(id)-1), 
        range = factor(range, ordered = TRUE, levels = 
                           c("(Semi)Cosmopolitan", "Holarctic", "Palaearctic", 
                             "European-Caucasian", "Mediterranean-Caucasian",
                             "Caspian and/or Caucasian", "unknown")), 
        id, Order, sp, abu, 
        fau = case_when(abu > 0 ~ 1, TRUE ~ 0)) %>% 
    group_by(seria, Order, range) %>% 
    summarise_if(is.numeric, sum) %>% 
    mutate(Community = abu/sum(abu)*100, Fauna = fau/sum(fau)*100, .keep = "unused") %>% 
    ungroup() %>% 
    pivot_longer(names_to = "VAR", values_to = "VAL", -1:-3) %>% 
    # filter(order  == "Oribatida & Mesostigmata") %>% 
    ggplot(aes(x = seria, y = VAL, fill = range)) + 
    geom_col(width = 0.68) + 
    facet_grid(cols = vars(VAR), rows = vars(Order)) + 
    scale_fill_manual(values = c(
        "#0D0786", 
        "#3FA9F5", 
        "#79D151", 
        "#8F2773", 
        "#FB9009", 
        "#FCE724", 
        alpha("white", 0))) +
    guides(fill = guide_legend(override.aes = list(alpha = 
                                                       c(1, 1, 1, 1, 1, 1, 0)))) + 
    scale_y_reverse() +
    scale_x_discrete(limits = xaxis) +
    labs(y = "Ratio, %", x = NULL, fill = NULL) + 
    theme(axis.text.x = element_text(angle = 90))
p8

if(export){
    ggsave(paste0("export/Fig 8. Range compound ", Sys.Date(), ".pdf"), 
           p8, height = 210/40, width = 297/40)
}


# + Рис. 9. ССА распределения видов -----------------------------------------
# Факторы: тип берега (? 1-4), 
# тип растительности (рогозы, злаки, ситники, хвощ, лох, тростники), 
# расстояние до моря (м), RH%, C%  

"https://uw.pressbooks.pub/appliedmultivariatestatistics/chapter/ca-dca-and-cca/"
"https://gist.github.com/perrygeo/7572735"
"https://pmassicotte.github.io/stats-denmark-2019/07_rda.html#/redundancy-analysis"
limit = 10
comm <- labs %>% 
    filter(seria != "Se") %>% 
    group_by(seria) %>% 
    summarise(
        plants_ = paste0(unique(plants.d), collapse = " / "), 
        coast_ = unique(coast), 
        RH = mean(RH), 
        CN = mean(CN),
        C = mean(C), 
        N = mean(N)) %>% 
    column_to_rownames("seria")

comm <- list(orws = select(orws, -Se, -PbTl, -PbPo),
             msws = select(msws, -Se, -DuCJ, -DuCS, -RsFd)
) %>%
    map(~ .x %>% 
            select_if(~ !is.numeric(.) || sum(.) > 1) %>% 
            mutate(total = apply(.[,-1], 1, sum)) %>% 
            filter(total > limit) %>% ###
            left_join(select(tables$tab.2_all_species, sp, No.), by = "sp")
    ) %>% 
    map(~ .x %>% 
            column_to_rownames("No.") %>%
            # column_to_rownames("sp") %>% 
            select_if(is.numeric) %>% 
            select(-total) %>% 
            t %>% 
            as.data.frame() 
    ) %>% 
    lapply(function(a){
        df0 <- comm[match(rownames(a), rownames(comm)),]
        if(prod(rownames(df0) == rownames(a)) == 1){ 
            vegan::cca(a, df0)
        } else {
            "Somethng is wrong"
        }
    })

if(export){
    pdf(paste0("export/Fig 9. Canonycal analysis_", Sys.Date(), ".pdf"), height = 7, width = 10)    
}
plot(comm[[1]], main = paste0("Oribatida, excl.: Se, PbTl, PbPo; limit = ", limit))
plot(comm[[2]], main = paste0("Mesostigmata, excl.: Se, DuCJ, DuCS, RsFd; limit = ", limit))
if(export){
    dev.off()
}

summary(comm[[1]])
summary(comm[[2]])


# + Рис. 10, таб. 5 Ординация мероценозов Mesostigmata и Oribatida  ----------------
# pcoa
PCOA <- dis %>% 
    lapply(function(a){
        p <- ape::pcoa(a)
        e <- p$values$Eigenvalues
        if(min(e) < 0){
            e <- e + abs(min(e))
            e <- round(e/sum(e)*100, 1)
        } else { 
            e <- round(e/sum(e)*100, 1)
        }
        p <- tibble::tibble(id = rownames(p$vectors), 
                            axis1 = p$vectors[,1], 
                            axis2 = p$vectors[,2]) 
        list(eig = e, pc = p)
    }) %>% 
    purrr::transpose()

M2 <- PCOA %>% 
    pluck("pc") %>% 
    map_df(rbind, .id = "D") %>% 
    filter(str_detect(id, "Se", negate = TRUE)) %>% 
    separate(D, into = c("taxa", "type")) %>% 
    left_join(labs, by = "id") %>% 
    select(taxa:coast) %>% 
    mutate(
        # axis1 = case_when(taxa == "ms" & type == "bin" ~ axis1*-1, TRUE ~ axis1),
        axis1 = case_when(taxa == "or" & type == "bin" ~ axis1*-1, TRUE ~ axis1),
        axis2 = case_when(taxa == "or" & type == "bin" ~ axis2*-1, TRUE ~ axis2),
        taxa = case_when(taxa == "or" ~ "Oribatida", TRUE ~ "Mesostigmata"), 
        type = case_when(type == "bin" ~ "Binary data (Jaccard)", 
                         TRUE ~ "Numeric data (Bray-Curtis)")) 

eig <- PCOA %>%
    pluck("eig") %>%
    map(~data.frame(axis1 = .x[1], axis2 = .x[2])) %>%
    map_df(rbind, .id = "a") %>%
    mutate_if(is.numeric, function(a){paste0(a, " %")}) %>%
    separate(a, into = c("taxa", "type")) %>% 
    mutate(taxa = case_when(taxa == "or" ~ "Oribatida", TRUE ~ "Mesostigmata"), 
           type = case_when(type == "bin" ~ "Binary data (Jaccard)", 
                            TRUE ~ "Numeric data (Bray-Curtis)"))


p10a <- M2 %>% 
    select(-id) %>% 
    distinct() %>% 
    ggplot(aes(x = axis1, y = axis2, color = plants.d)) + 
    geom_point() + 
    stat_ellipse() +
    geom_text(aes(label = axis1, x = 0, y = -0.77), color = "black",alpha = 0.68, 
              data = eig, ) +
    geom_text(aes(label = axis2, x = -0.99, y = 0), color = "black", alpha = 0.68, 
              data = eig, angle = 90) +
    facet_grid(cols = vars(type), rows = vars(taxa)) + 
    labs(x = NULL, y = NULL, color = NULL, #subtitle = "Б. Доминантные виды растений") + 
         subtitle = "A. Dominant plant species") + 
    theme(legend.text = element_text(face = "italic"))
p10b <- M2 %>% 
    select(-id) %>% 
    distinct() %>% 
    ggplot(aes(x = axis1, y = axis2, color = coast)) + 
    geom_point() + 
    stat_ellipse() +
    geom_text(aes(label = axis1, x = 0, y = -0.63), color = "black", data = eig, alpha = 0.68) +
    geom_text(aes(label = axis2, x = -0.9, y = 0), color = "black", data = eig, alpha = 0.68, angle = 90) +
    scale_color_manual(values = c("#F8766D", "#00A9FF", "#0CB720", "#CD9600"))+
    facet_grid(cols = vars(type), rows = vars(taxa)) + 
    labs(x = NULL, y = NULL, color  = NULL, #subtitle = "A. Тип берега")
         subtitle = "B. Coast type") 
p10 <- gridExtra::grid.arrange(p10a, p10b, ncol = 1) #, left = "Суммарное обилие в серии, экз.")
if(export){
    ggsave(paste0("export/Fig 10. Ordination ", Sys.Date(), ".pdf"), 
           p10, width = 210/25, height = 297/25)
}

# PERMANOVA
PERMANOVA <- expand_grid(id = names(dis),
                         type = c("coast + plants.d", "plants.d + coast")) %>% 
    mutate(no = rep(1:4, each = 2), 
           nm = paste0(id, " ~ ", type))

P <- PERMANOVA %>% 
    split(1:nrow(.)) %>% 
    lapply(function(a){
        dis_tmp <- dis[[a$no]]
        labs_tmp <- dis_tmp %>% 
            labels %>% 
            tibble(id = .) %>% 
            left_join(labs, by = "id")
        vegan::adonis2(
            formula = as.formula(paste0("dis_tmp ~", a$type)), 
            data = labs_tmp, permutations = 999)
    }) %>% 
    `names<-`(PERMANOVA$nm)

tables$tab.5a_permanova <- P %>% 
    map(~ as.data.frame(.x)) %>% 
    map_dfr(rbind, .id = "id") %>% 
    rownames_to_column("type") %>% 
    separate(type, into = "type", sep = "\\.\\.", extra = "drop") %>% 
    transmute(
        id, type, Df, SumOfSqs = round(SumOfSqs, 1), 
        R2 = round(R2*100, 1), `F` = round(F, 1), p.value = `Pr(>F)`)

# tables$tab.5b_permanova <- tables$tab.5a_permanova %>% 
#     mutate(taxa_type = substr(id, 1, 6)) %>% 
#     split(.$taxa_type) %>% 
#     map(~.x a %>% 
#             filter(type != "Total") %>% 
#             group_by(type) %>% 
#             summarise(r = min(R2)) %>% 
#             pivot_wider(names_from = type, values_from = r) %>% 
#             transmute(
#                 coast, 
#                 plants.d, 
#                 variable = 100 - Residual - plants.d - coast, 
#                 Residual
#             )
#     ) %>% 
#     map_dfr(rbind, .id = "id")

# export tables -----------------------------------------------------------
tables %>% 
    .[str_detect(names(tables), "tab")] %>% 
    writexl::write_xlsx(paste0("export/samoor_tables_", Sys.Date(), ".xlsx"))

long %>% 
    filter(str_detect(sp, "Zercon"), abu > 0)


