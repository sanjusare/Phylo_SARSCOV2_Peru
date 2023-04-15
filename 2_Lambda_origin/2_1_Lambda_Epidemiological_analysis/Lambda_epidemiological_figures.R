# 1:set working directory####
args <- commandArgs(trailingOnly = TRUE)
# 2:load libraries and establish settings#### function for
# integer_breaks (from jhrcook github)
library(scales)
integer_breaks <- function(n = 5, ...) {
    fxn <- function(x) {
        breaks <- floor(pretty(x, n, ...))
        names(breaks) <- attr(breaks, "labels")
        breaks
    }
    return(fxn)
}
library(dplyr)
library(lubridate)
library(ggplot2)
library(tidyr)
library(ggpubr)
theme_set(theme_bw(base_size = 15))
Sys.setlocale("LC_TIME", "English")
cols <- c(Omicron = "#db8a8a", Gamma = "#dbd48a", Lambda = "#8fdb8a",
    Delta = "#8ad3db", Mu = "#a18adb", Alpha = "#db8ac0", Other = "darkgrey")
colors <- c("#cbcbf7", "#f7cbcb", "#f7f4cb", "thistle1", "#91e3db",
    "#abc4d9", "#9c9181", "aquamarine4", "#a9c7af", "#99a5d6",
    "#f4b771", "#ec74dc", "#9ed5ff", "#dedede")
count_colors <- c(Argentina = "dodgerblue", Chile = "forestgreen",
    Peru = "firebrick", Colombia = "orangered", Ecuador = "goldenrod",
    Mexico = "darkmagenta")
# 3:read files####
metadata <- read.csv(args[1], sep = "\t", header = TRUE)
# download from
# https://github.com/owid/covid-19-data/tree/master/public/data
cas_wor <- read.csv(args[2])
# 4:Plot of Lambda genomes (world vs Peru)####
lams <- metadata %>%
    dplyr::select(name_analysis, pango_lineage, variant, Country,
        date) %>%
    dplyr::filter(grepl("Lambda", variant)) %>%
    tidyr::drop_na(date) %>%
    dplyr::mutate(date = as.Date(date)) %>%
    dplyr::filter(date >= as.Date("2020-06-01")) %>%
    dplyr::filter(date <= as.Date("2021-12-31")) %>%
    dplyr::mutate(week = cut.Date(date, breaks = "1 week")) %>%
    dplyr::mutate(week = as.Date(week)) %>%
    group_by(Country, week) %>%
    summarise(count = n())
lams$Peru <- "World"
lams$Peru[lams$Country == "Peru"] <- "Peru"
lams_comp <- ggplot(lams, aes(x = as.POSIXct(week), fill = Peru)) +
    geom_bar(aes(y = count), stat = "identity") + facet_wrap(~Peru,
    ncol = 1) + theme(legend.position = "none") + scale_x_datetime(date_breaks = "1 month",
    date_labels = "%Y-%b", limits = c(as.POSIXct("2020-06-01"),
        as.POSIXct("2021-01-30"))) + theme(axis.text.x = element_text(angle = 90,
    hjust = 0, vjust = 0), axis.title.x = element_blank(), legend.position = "none") +
    ylab("Number of Genomes") + scale_fill_manual(values = c(Peru = "forestgreen",
    World = "black"))
ggsave("First_Lambda_genomes.jpg", plot = lams_comp, dpi = 500,
    height = 3, width = 3)
# 5:Determine the countries with at least two cities with
# 15 % Lambda prevalence before April 2021 and make
# plots#### 5.1:cities with at least one Lambda genome
# between june 2020 and may 2021
metadata$City <- gsub("Amazona", "Amazonas", metadata$City)
metadata$City <- gsub("Amazonass", "Amazonas", metadata$City)
metadata$City <- gsub("CiudadAutonomadeBuenosAires", "BuenosAires",
    metadata$City)
metadata$City <- gsub("CiudadAutónomadeBuenosAires", "BuenosAires",
    metadata$City)
metadata$City <- gsub("CityofBuenosAires", "BuenosAires", metadata$City)
metadata$City <- gsub("Eloro", "ElOro", metadata$City)
metadata$City <- gsub("Ñuble", "Nuble", metadata$City)
metadata$City <- gsub("RegionMetropolitanadeSantiago", "Santiago",
    metadata$City)
metadata$City <- gsub("RíoNegro", "RioNegro", metadata$City)
metadata$City <- gsub("StateofMexico", "MexicoCity", metadata$City)
metadata$City <- gsub("ValleDelCauca", "ValledelCauca", metadata$City)
lam_city <- metadata %>%
    filter(grepl("Lambda", variant)) %>%
    tidyr::drop_na(date) %>%
    dplyr::mutate(date = as.Date(date)) %>%
    dplyr::filter(date >= as.Date("2020-06-01")) %>%
    dplyr::filter(date <= as.Date("2021-12-31")) %>%
    dplyr::mutate(City_con = paste(City, Country, sep = "_")) %>%
    drop_na(City) %>%
    dplyr::select(name_analysis, pango_lineage, variant, City_con,
        date) %>%
    group_by(City_con) %>%
    summarise(count = n())
# 5.2:Count the genomes by city by week from the cities
# previously selected in step 5.1 (lam_city$City_con)
gen_num_city <- metadata %>%
    dplyr::mutate(City_con = paste(City, Country, sep = "_")) %>%
    dplyr::filter(City_con %in% lam_city$City_con) %>%
    tidyr::drop_na(date) %>%
    dplyr::mutate(date = as.Date(date)) %>%
    dplyr::filter(date >= as.Date("2020-06-01")) %>%
    dplyr::filter(date <= as.Date("2021-12-31")) %>%
    dplyr::mutate(week = cut.Date(date, breaks = "1 week")) %>%
    dplyr::group_by(City_con, week) %>%
    dplyr::summarise(tot = n()) %>%
    dplyr::ungroup() %>%
    tidyr::complete(week, nesting(City_con), fill = list(tot = 0))
# 5.3:calculate the relative prevalence of Lambda by week
# from the cities selected previously in step 5.1
# (lam_city$City_con)
rel_prev_city <- metadata %>%
    dplyr::mutate(City_con = paste(City, Country, sep = "_")) %>%
    dplyr::filter(City_con %in% lam_city$City_con) %>%
    tidyr::drop_na(date) %>%
    dplyr::mutate(date = as.Date(date)) %>%
    dplyr::filter(date >= as.Date("2020-06-01")) %>%
    dplyr::filter(date <= as.Date("2021-12-31")) %>%
    dplyr::mutate(Lambda = "yes") %>%
    dplyr::mutate(Lambda = ifelse(grepl("Lambda", variant), "yes",
        "no")) %>%
    dplyr::mutate(week = cut.Date(date, breaks = "1 week")) %>%
    dplyr::group_by(City_con, week, Lambda) %>%
    dplyr::summarise(count = n()) %>%
    dplyr::mutate(relprev = count/sum(count)) %>%
    dplyr::ungroup() %>%
    tidyr::complete(week, Lambda, City_con, fill = list(count = 0,
        relprev = 0)) %>%
    dplyr::full_join(gen_num_city) %>%
    dplyr::mutate(CI = 1.96 * (sqrt(((1 - relprev) * (relprev))/tot))) %>%
    dplyr::mutate(week = as.Date(week), City = sapply(strsplit(City_con,
        "_"), "[", 1), Country = sapply(strsplit(City_con, "_"),
        "[", 2))
# 5.4:select the countries with at least two cities with 15
# % prevalence of Lambda before April
countries_list <- list()
for (i in 1:length(levels(as.factor(rel_prev_city$Country)))) {
    print(i)
    query <- levels(as.factor(rel_prev_city$Country))[i]
    df <- rel_prev_city %>%
        dplyr::filter(week <= "2021-04-25") %>%
        dplyr::filter(Country == query)
    cities_n <- length(levels(as.factor(df$City)))
    if (cities_n > 2) {
        dftest <- df %>%
            dplyr::filter(Lambda == "yes")
        if (max(dftest$relprev) >= 0.15) {
            dftest2 <- dftest %>%
                dplyr::filter(relprev >= 0.15)
            if (length(levels(as.factor(dftest2$City))) >= 2) {
                countries_list[[query]] <- query
            }
        }
    }
}
# 5.5:Select the cities from the selected countries
# (countries_list) that reached 15 % Lambda prevalence
# before April 2021
rel_prev_city_sel <- rel_prev_city %>%
    dplyr::filter(Country %in% countries_list) %>%
    dplyr::arrange(Country)
cities_selected <- c()
for (i in levels(as.factor(rel_prev_city_sel$City_con))) {
    dt <- rel_prev_city_sel %>%
        dplyr::filter(City_con == i & Lambda == "yes") %>%
        dplyr::filter(week <= "2021-04-25")
    if (max(dt$relprev) > 0.15) {
        cities_selected <- append(cities_selected, i)
    }
}
rel_prev_city_sel_p <- rel_prev_city_sel %>%
    dplyr::filter(City_con %in% cities_selected) %>%
    dplyr::arrange(Country)
# 5.6:Calculation of the relative prevalence in each
# selected city (cities_selected)
p_countries_prev_list <- list()
for (i in levels(as.factor(rel_prev_city_sel_p$Country))) {
    print(i)
    dcon <- rel_prev_city_sel_p %>%
        dplyr::filter(Country == i & Lambda == "yes")
    p_countries_prev_list[[i]] <- ggplot(dcon, (aes(x = as.POSIXct(week)))) +
        geom_point(aes(y = relprev, color = Country)) + scale_color_manual(values = count_colors) +
        geom_errorbar(aes(ymin = pmax(relprev - CI, 0), ymax = pmin(1,
            relprev + CI))) + facet_wrap(~City) + ylim(c(0, 1)) +
        labs(y = "Relative prevalence") + scale_x_datetime(date_breaks = "1 month",
        date_labels = "%Y-%b", limits = c(as.POSIXct("2020-08-01"),
            as.POSIXct("2021-04-25"))) + theme(axis.text.x = element_text(angle = 90,
        hjust = 0, vjust = 0), axis.title.x = element_blank(),
        legend.position = "none") + ggtitle(i)
}
# 5.7: Plot of Lambda relative prevalence in the selected
# cities
p1 <- ggarrange(p_countries_prev_list[["Argentina"]], p_countries_prev_list[["Chile"]])
p2 <- ggarrange(p_countries_prev_list[["Colombia"]], p_countries_prev_list[["Ecuador"]],
    p_countries_prev_list[["Mexico"]], ncol = 3, nrow = 1)
p_fin <- ggarrange(p1, p2, p_countries_prev_list[["Peru"]], ncol = 1,
    nrow = 3, heights = c(0.38, 0.23, 0.46))
ggsave("Cities_Lambda_relative_prevalence.jpg", plot = p_fin,
    dpi = 500, height = 17.5, width = 15)
# 5.8: Extraction of number of Genomes in each selected
# city (cities_selected)
p_countries_gen_list <- list()
for (i in levels(as.factor(rel_prev_city_sel_p$Country))) {
    print(i)
    dcon <- rel_prev_city_sel_p %>%
        dplyr::filter(Country == i & Lambda == "yes")
    p_countries_gen_list[[i]] <- ggplot(dcon, (aes(x = as.POSIXct(week)))) +
        geom_bar(aes(y = tot), alpha = 0.5, stat = "identity") +
        geom_bar(aes(y = count, fill = Country), alpha = 0.5,
            stat = "identity") + scale_fill_manual(values = count_colors) +
        facet_wrap(~City, scales = "free_y") + labs(y = "Genome count") +
        scale_x_datetime(date_breaks = "1 month", date_labels = "%Y-%b",
            limits = c(as.POSIXct("2020-08-01"), as.POSIXct("2021-04-25"))) +
        scale_y_continuous(breaks = integer_breaks(n = 3)) +
        theme(axis.text.x = element_text(angle = 90, hjust = 0,
            vjust = 0), axis.title.x = element_blank(), legend.position = "none") +
        ggtitle(i)
}
# 5.9: Plot of number of Lambda genomes and total in the
# selected cities
p1_g <- ggarrange(p_countries_gen_list[["Argentina"]], p_countries_gen_list[["Chile"]])
p2_g <- ggarrange(p_countries_gen_list[["Colombia"]], p_countries_gen_list[["Ecuador"]],
    p_countries_gen_list[["Mexico"]], ncol = 3, nrow = 1)
p_fin_g <- ggarrange(p1_g, p2_g, p_countries_gen_list[["Peru"]],
    ncol = 1, nrow = 3, heights = c(0.38, 0.23, 0.46))
ggsave("Cities_Lambda_genomes.jpg", plot = p_fin_g, dpi = 500,
    height = 17.5, width = 17.5)
# 6:Adjustment of the Lambda prevalence by country using
# loess####
log_adj_c <- function(x, voc, from, to) {
    df_p <- metadata %>%
        dplyr::filter(Country == x) %>%
        drop_na(date) %>%
        dplyr::filter(date >= as.Date(from)) %>%
        dplyr::filter(date <= as.Date(to))
    df_p$month <- as.numeric(sapply(strsplit(as.character(df_p$date),
        "-"), "[", 2))
    df_p$M <- as.character(month(ymd(10101) + months(df_p$month -
        1), label = TRUE, abbr = TRUE))
    df_p$year <- sapply(strsplit(as.character(df_p$date), "-"),
        "[", 1)
    df_p$y_m <- paste(df_p$year, df_p$M, sep = "-")
    df_p$VOC <- sapply(strsplit(df_p$variant, " "), "[", 2)
    df_p$VOC[is.na(df_p$VOC)] <- "Other"
    for (i in 1:nrow(df_p)) {
        if (!(df_p$VOC[i] %in% c("Lambda", "Gamma", "Alpha",
            "Delta", "Mu", "Omicron", "Other"))) {
            df_p$VOC[i] <- "Other"
        }
    }
    df_p$date <- as.Date(df_p$date)
    df_p_week <- df_p %>%
        dplyr::mutate(week = cut.Date(date, breaks = "1 week",
            labels = FALSE), collection_date = cut.Date(date,
            breaks = "1 week")) %>%
        dplyr::arrange(date)
    df_p_cg <- df_p_week %>%
        dplyr::group_by(week, VOC) %>%
        dplyr::summarise(count = n()) %>%
        dplyr::mutate(percentage = count/sum(count))
    if (length(setdiff(1:max(df_p_cg$week), unique(df_p_cg$week))) >
        0) {
        coms_weeks <- data.frame(week = setdiff(1:max(df_p_cg$week),
            unique(df_p_cg$week)), VOC = "Other", count = 0,
            percentage = 0)
        df_p_com <- df_p_cg %>%
            dplyr::ungroup() %>%
            full_join(coms_weeks) %>%
            tidyr::complete(VOC, nesting(week), fill = list(count = 0,
                percentage = 0))
    } else {
        df_p_com <- df_p_cg %>%
            dplyr::ungroup() %>%
            tidyr::complete(VOC, nesting(week), fill = list(count = 0,
                percentage = 0))
    }
    df_p_com_fin <- df_p_com %>%
        dplyr::group_by(week) %>%
        dplyr::mutate(tot_gen_week = sum(count))
    df_p_com_fin$CI <- 1.96 * (sqrt(((1 - df_p_com_fin$percentage) *
        (df_p_com_fin$percentage))/df_p_com_fin$tot_gen_week))
    df_p_com_fin$CI[is.nan(df_p_com_fin$CI)] <- 0
    df_voc_1 <- df_p_com_fin %>%
        dplyr::filter(VOC == voc)
    weeks_index <- df_p_week %>%
        dplyr::select(week, collection_date) %>%
        unique() %>%
        tidyr::complete(collection_date) %>%
        dplyr::mutate(week = 1:nrow(df_voc_1))
    df_voc <- df_voc_1 %>%
        dplyr::full_join(weeks_index)
    fit <- loess(percentage ~ week, degree = 1, span = 21/nrow(df_voc),
        data = df_voc)
    loessfit <- data.frame(week = 1:length(fit$fitted), fitted = fit$fitted)
    f_df <- Reduce(function(x, y) dplyr::full_join(x, y, all = TRUE),
        list(df_voc, loessfit)) %>%
        mutate(collection_date = as.Date(collection_date), Country = x)
    return(f_df)
}
Peru_lambda <- log_adj_c("Peru", "Lambda", "2020-06-01", "2021-12-31")
Arge_lambda <- log_adj_c("Argentina", "Lambda", "2020-06-01",
    "2021-12-31")
Chil_lambda <- log_adj_c("Chile", "Lambda", "2020-06-01", "2021-12-31")
Colo_lambda <- log_adj_c("Colombia", "Lambda", "2020-06-01",
    "2021-12-31")
Ecua_lambda <- log_adj_c("Ecuador", "Lambda", "2020-06-01", "2021-12-31")
Mexi_lambda <- log_adj_c("Mexico", "Lambda", "2020-06-01", "2021-12-31")
all_lambda <- Reduce(function(x, y) dplyr::full_join(x, y, all = TRUE),
    list(Peru_lambda, Arge_lambda, Chil_lambda, Colo_lambda,
        Ecua_lambda, Mexi_lambda))
# 7:Calculation of the number of cases based on the
# adjusted prevalence from step 6####
cases_calc <- function(x, from, to) {
    cas_filt <- cas_wor %>%
        dplyr::filter(location == x) %>%
        dplyr::select(country = location, date, total_cases,
            new_cases) %>%
        dplyr::mutate(date = as.Date(date)) %>%
        dplyr::filter(date >= as.Date(from)) %>%
        dplyr::filter(date <= as.Date(to))
    cas_count <- cas_filt %>%
        dplyr::select(date = date, confirm = new_cases, country) %>%
        dplyr::arrange(date)
    cas_count$Date <- as.Date(cas_count$date, format = "%Y-%m-%d")
    weeks_index <- cas_count %>%
        dplyr::mutate(week = cut.Date(Date, breaks = "1 week",
            labels = FALSE)) %>%
        dplyr::mutate(collection_date = cut.Date(Date, breaks = "1 week")) %>%
        dplyr::select(week, collection_date) %>%
        unique()
    rep_week <- cas_count %>%
        dplyr::mutate(week = cut.Date(Date, breaks = "1 week",
            labels = FALSE)) %>%
        dplyr::select(country, Date, confirm, week)
    rep_week$confirm[is.na(rep_week$confirm)] <- 0
    rep_week_count_com <- rep_week %>%
        complete(Date, nesting(country), fill = list(confirm = 0)) %>%
        tidyr::fill(week)
    rep_week_count_com$week <- as.numeric(rep_week_count_com$week)
    case_roll <- rep_week_count_com %>%
        dplyr::mutate(roll_14 = zoo::rollmean(confirm, k = 14,
            fill = 0)) %>%
        dplyr::group_by(week, country) %>%
        dplyr::summarise(confirm = sum(confirm), roll_14 = sum(roll_14)) %>%
        dplyr::ungroup() %>%
        dplyr::full_join(weeks_index)
    case_roll$roll_14[case_roll$week == 1] <- case_roll$confirm[case_roll$week ==
        1]
    case_roll$roll_14[case_roll$week == nrow(case_roll)] <- case_roll$confirm[case_roll$week ==
        nrow(case_roll)]
    case_roll$roll_14[case_roll$week == nrow(case_roll) - 1] <- case_roll$confirm[case_roll$week ==
        nrow(case_roll) - 1]
    case_roll_f <- case_roll %>%
        dplyr::select(Country = country, cases = confirm, roll_14,
            collection_date) %>%
        dplyr::mutate(collection_date = as.Date(collection_date))
    return(case_roll_f)
}
Per_lambda_cases <- cases_calc("Peru", "2020-06-01", "2021-12-31")
Chi_lambda_cases <- cases_calc("Chile", "2020-06-01", "2021-12-31")
Arg_lambda_cases <- cases_calc("Argentina", "2020-06-01", "2021-12-31")
Col_lambda_cases <- cases_calc("Colombia", "2020-06-01", "2021-12-31")
Ecu_lambda_cases <- cases_calc("Ecuador", "2020-06-01", "2021-12-31")
Mex_lambda_cases <- cases_calc("Mexico", "2020-06-01", "2021-12-31")
all_cases_lambda <- Reduce(function(x, y) dplyr::full_join(x,
    y, all = TRUE), list(Per_lambda_cases, Arg_lambda_cases,
    Chi_lambda_cases, Col_lambda_cases, Ecu_lambda_cases, Mex_lambda_cases))
# 8:Calculation of sampling proportion, and smoothing of
# Lambda cases####
final_data <- full_join(all_lambda, all_cases_lambda) %>%
    dplyr::mutate(sampl_prop_cas = tot_gen_week/cases, sampl_prop_roll = tot_gen_week/roll_14,
        Lambda_cases = fitted * roll_14)
# 9:Plot of smoothed relative prevalences and estimated
# Lambda cases by countries####
p_fit_rel_prev <- ggplot(final_data, aes(x = as.POSIXct(collection_date))) +
    geom_point(aes(color = Country, y = percentage), alpha = 0.5,
        size = 1) + geom_errorbar(aes(ymin = pmax(0, percentage -
    CI), ymax = pmin(percentage + CI, 1))) + scale_x_datetime(date_breaks = "1 month",
    date_labels = "%Y-%b", limits = c(as.POSIXct("2020-06-01"),
        as.POSIXct("2021-04-30"))) + theme(axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
        size = 8), axis.text.y = element_text(size = 8), legend.position = "none",
    strip.background = element_rect(fill = "white")) + labs(y = "Relative prevalence") +
    geom_smooth(aes(color = Country, fill = Country, y = percentage),
        method = "loess", span = 21/nrow(Peru_lambda), method.args = list(degree = 1),
        alpha = 0.2, size = 0.5) + scale_color_manual(values = count_colors) +
    scale_fill_manual(values = count_colors) + facet_wrap(~Country,
    ncol = 1) + scale_y_continuous(limits = c(0, 1)) + geom_vline(xintercept = as.POSIXct("2020-11-01"),
    linetype = "dotted") + geom_vline(xintercept = as.POSIXct("2020-12-31"),
    linetype = "dotted")
p_est_lam_cas <- ggplot(final_data, aes(x = as.POSIXct(collection_date))) +
    geom_line(aes(color = Country, y = Lambda_cases)) + scale_x_datetime(date_breaks = "1 month",
    date_labels = "%Y-%b", limits = c(as.POSIXct("2020-10-01"),
        as.POSIXct("2021-01-31"))) + ylim(c(0, 15000)) + theme(axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
        size = 8), axis.text.y = element_text(size = 8), legend.position = "none",
    strip.background = element_rect(fill = "white")) + labs(y = "Estimated Lambda cases") +
    scale_color_manual(values = count_colors) + facet_wrap(~Country,
    ncol = 1) + geom_vline(xintercept = as.POSIXct("2020-11-01"),
    linetype = "dotted") + geom_vline(xintercept = as.POSIXct("2020-12-31"),
    linetype = "dotted")
p_rel_prev_lam_cas <- ggarrange(p_fit_rel_prev, p_est_lam_cas,
    nrow = 1, ncol = 2, labels = c("A", "B"))
ggsave(filename = "Fitted_relprev_lam_cases_country.jpg", plot = p_rel_prev_lam_cas,
    dpi = 500, height = 8, width = 8)
# 10:Plot of Number of total and Lambda genomes and
# sampling proprotions by countries####
p_lam_gen <- ggplot(subset(final_data, collection_date >= "2020-10-01" &
    collection_date < "2021-04-30"), aes(x = as.POSIXct(collection_date))) +
    geom_bar(aes(y = tot_gen_week), alpha = 0.5, stat = "identity") +
    geom_bar(aes(y = count, fill = Country), alpha = 0.5, stat = "identity") +
    scale_x_datetime(date_breaks = "1 month", date_labels = "%Y-%b",
        limits = c(as.POSIXct("2020-10-01"), as.POSIXct("2021-04-30"))) +
    theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45,
        hjust = 1, vjust = 1, size = 8), axis.text.y = element_text(size = 8),
        legend.position = "none", strip.background = element_rect(fill = "white")) +
    labs(y = "Number of Genomes") + scale_color_manual(values = count_colors) +
    scale_fill_manual(values = count_colors) + facet_wrap(~Country,
    ncol = 1, scales = "free_y")
p_samp_prop <- ggplot(subset(final_data, collection_date >= "2020-10-01" &
    collection_date < "2021-04-30"), aes(x = as.POSIXct(collection_date))) +
    geom_line(aes(group = Country, y = sampl_prop_roll * 100)) +
    scale_x_datetime(date_breaks = "1 month", date_labels = "%Y-%b",
        limits = c(as.POSIXct("2020-10-01"), as.POSIXct("2021-04-30"))) +
    theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45,
        hjust = 1, vjust = 1, size = 8), axis.text.y = element_text(size = 8),
        legend.position = "none", strip.background = element_rect(fill = "white")) +
    labs(y = "Sampling proportion (%)") + facet_wrap(~Country,
    ncol = 1, scales = "free_y")
p_lam_gen_sam_prop <- ggarrange(p_lam_gen, p_samp_prop, nrow = 1,
    ncol = 2, labels = c("A", "B"))
ggsave(filename = "Lambda_genomes_sampprop_country.jpg", plot = p_lam_gen_sam_prop,
    dpi = 500, height = 8, width = 8)
# 11:Save the processed data of Lambda epidemiology in
# selected countries
write.table(final_data, "Lambda_epidemiological_data.tsv", quote = FALSE,
    sep = "\t", col.names = TRUE, row.names = FALSE)
