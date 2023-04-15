# 1:set working directory####
args <- commandArgs(trailingOnly = TRUE)
# 2:load libraries and establish settings####
library(dplyr)
library(ggplot2)
library(treeio)
library(tidyr)
library(ggnewscale)
theme_set(theme_bw(base_size = 12))
Sys.setlocale("LC_TIME", "English")
colreg <- c(center = "firebrick", north = "dodgerblue", south = "forestgreen",
    `mid-east` = "darkmagenta", `north-east` = "goldenrod", `south-east` = "hotpink3")
collin <- c(SubL1 = "#58578a", SubL2 = "#8ad3db")
# 3:read files#### read the ordered data (from step 1)
tabla <- read.csv(file = args[1], sep = "\t", header = TRUE)
# download from
# https://www.datosabiertos.gob.pe/dataset/casos-positivos-por-covid-19-ministerio-de-salud-minsa
cas_city <- read.csv(args[2], sep = ";")
# read files of regions
reg <- read.csv(args[3], sep = "\t")
# read sublineages classifications
sel_lins <- list()
for (i in 1:2) {
    sel_lins[[i]] <- readLines(paste("sublineage_SubL", i, ".list",
        sep = ""))
}
# 4:Analysis of prevalences of the sublineages by
# region#### data preparations
df_p <- tabla %>%
    dplyr::filter(Country == "Peru") %>%
    dplyr::mutate(date = as.Date(date)) %>%
    drop_na(date) %>%
    dplyr::filter(date >= as.Date("2020-08-01")) %>%
    dplyr::filter(date <= as.Date("2022-01-31"))
df_p$SubLs <- "no"
for (i in 1:length(sel_lins)) {
    df_p$SubLs[df_p$name_analysis %in% sel_lins[[i]]] <- paste("SubL",
        i, sep = "")
}
df_p_week <- df_p %>%
    dplyr::mutate(week = cut.Date(date, breaks = "1 week", labels = FALSE),
        collection_date = cut.Date(date, breaks = "1 week")) %>%
    dplyr::arrange(date)
df_p_week$City[df_p_week$City == "Amazona"] <- "Amazonas"
df_p_week$City[df_p_week$City == "Amazonass"] <- "Amazonas"
df_p_week$City[df_p_week$City == "HuancavelIca"] <- "Huancavelica"
df_p_week$City[df_p_week$City == "Huánuco"] <- "Huanuco"
df_p_week$City[df_p_week$City == "Lalibertad"] <- "LaLibertad"
df_p_week$City[df_p_week$City == "MadreDeDios"] <- "MadredeDios"
df_p_reg <- left_join(df_p_week, reg) %>%
    tidyr::drop_na(City)
df_p_cg <- df_p_reg %>%
    dplyr::group_by(week, Region, SubLs) %>%
    dplyr::summarise(count = n()) %>%
    dplyr::mutate(percentage = count/sum(count))
if (length(setdiff(1:max(df_p_cg$week), unique(df_p_cg$week))) >
    0) {
    coms_weeks <- data.frame(week = setdiff(1:max(df_p_cg$week),
        unique(df_p_cg$week)), SubLs = "no", count = 0, percentage = 0)
    df_p_com <- df_p_cg %>%
        dplyr::ungroup() %>%
        full_join(coms_weeks) %>%
        tidyr::complete(SubLs, week, nesting(Region), fill = list(count = 0,
            percentage = 0))
} else {
    df_p_com <- df_p_cg %>%
        dplyr::ungroup() %>%
        tidyr::complete(SubLs, week, nesting(Region), fill = list(count = 0,
            percentage = 0))
}
df_p_com_fin <- df_p_com %>%
    dplyr::group_by(week, Region) %>%
    dplyr::mutate(tot_gen_week = sum(count))
df_p_com_fin$CI <- 1.96 * (sqrt(((1 - df_p_com_fin$percentage) *
    (df_p_com_fin$percentage))/df_p_com_fin$tot_gen_week))
df_p_com_fin$CI[is.nan(df_p_com_fin$CI)] <- 0
df_p_com_fin$est <- "No"
df_p_com_fin$est[df_p_com_fin$tot_gen_week == 0] <- "Yes"
f_df_list <- list()
for (i in levels(as.factor(df_p_com$Region))) {
    for (w in levels(as.factor(df_p_com$SubLs))) {
        if (w != "no") {
            df_voc_1 <- df_p_com_fin %>%
                dplyr::filter(SubLs == w & Region == i)
            for (z in 2:nrow(df_voc_1)) {
                if (df_voc_1$est[z] == "No") {
                  break
                }
            }
            for (x in z:nrow(df_voc_1)) {
                if (df_voc_1$est[x] == "Yes") {
                  one <- df_voc_1$percentage[x - 1]
                  for (y in x:nrow(df_voc_1)) {
                    if (df_voc_1$est[y] == "No") {
                      break
                    }
                  }
                  sec <- df_voc_1$percentage[y]
                  df_voc_1$percentage[x] <- (one + sec)/2
                }
            }
            weeks_index <- df_p_week %>%
                dplyr::select(week, collection_date) %>%
                unique() %>%
                tidyr::complete(collection_date)
            weeks_index <- weeks_index %>%
                dplyr::mutate(week = 1:nrow(weeks_index))
            df_voc <- df_voc_1[z:nrow(df_voc_1), ] %>%
                dplyr::left_join(weeks_index) %>%
                dplyr::arrange(week)
            fit <- loess(percentage ~ week, degree = 1, span = 21/nrow(df_voc),
                data = df_voc)
            fit_values <- predict(fit, data.frame(week = z:nrow(df_voc_1)),
                se = T)
            loessfit <- data.frame(week = z:nrow(df_voc_1), fitted = fit_values[["fit"]],
                se_fitted = fit_values[["se.fit"]])
            f_df_list[[paste(i, w, sep = "_")]] <- Reduce(function(x,
                y) dplyr::full_join(x, y, all = TRUE), list(df_voc,
                loessfit)) %>%
                mutate(collection_date = as.Date(collection_date))
            print(paste(i, w, sep = "_"))
        }
    }
}
f_df <- do.call(rbind, f_df_list) %>%
    mutate(fitted = pmax(0, fitted))
# plot of sublineages prevalences by regions
rel_prev_all <- ggplot(subset(f_df), aes(x = as.POSIXct(collection_date))) +
    geom_point(data = subset(f_df, est == "No" & percentage >
        0), aes(y = percentage, color = SubLs), show.legend = F) +
    geom_errorbar(data = subset(f_df, est == "No" & percentage >
        0), aes(ymin = pmax(0, percentage - CI), ymax = pmin(1,
        percentage + CI)), color = "lightgrey") + scale_color_manual(values = collin) +
    new_scale_color() + new_scale_fill() + geom_line(data = subset(f_df),
    aes(y = fitted, color = SubLs)) + geom_ribbon(data = subset(f_df),
    aes(fill = SubLs, ymin = pmax(0, fitted - se_fitted), ymax = pmin(1,
        fitted + se_fitted)), alpha = 0.5) + scale_color_manual(values = collin) +
    scale_fill_manual(values = collin) + facet_wrap(~Region) +
    ylim(c(0, 1)) + ylab("Relative prevalence") + scale_x_datetime(date_breaks = "1 month",
    date_labels = "%Y-%b", limits = c(as.POSIXct("2020-10-01"),
        as.POSIXct("2022-01-31"))) + theme(axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5),
    legend.position = "none", strip.background = element_rect(fill = "white"))
ggsave("Relprev_Lambda_sublineages_reg.jpg", plot = rel_prev_all,
    dpi = 300, height = 5, width = 7.5)
# 5:Overall prevalence of Lambda#### data preparations
df_p_lam <- df_p_reg %>%
    dplyr::mutate(Lambda = "no")
df_p_lam$Lambda[grepl("Lambda", df_p_lam$variant)] <- "yes"
df_p_lam_gr_tot <- df_p_lam %>%
    dplyr::group_by(week, Lambda) %>%
    dplyr::summarise(count = n()) %>%
    dplyr::mutate(percentage = count/sum(count))
if (length(setdiff(1:max(df_p_lam_gr_tot$week), unique(df_p_lam_gr_tot$week))) >
    0) {
    coms_weeks <- data.frame(week = setdiff(1:max(df_p_lam_gr_tot$week),
        unique(df_p_lam_gr_tot$week)), Lambda = "no", count = 0,
        percentage = 0)
    df_p_com_lam_tot <- df_p_lam_gr_tot %>%
        dplyr::ungroup() %>%
        full_join(coms_weeks) %>%
        tidyr::complete(Lambda, week, fill = list(count = 0,
            percentage = 0))
} else {
    df_p_com_lam_tot <- df_p_lam_gr_tot %>%
        dplyr::ungroup() %>%
        tidyr::complete(Lambda, week, fill = list(count = 0,
            percentage = 0))
}
df_p_com_lam_tot_fin <- df_p_com_lam_tot %>%
    dplyr::group_by(week) %>%
    dplyr::mutate(tot_gen_week = sum(count))
df_p_com_lam_tot_fin$CI <- 1.96 * (sqrt(((1 - df_p_com_lam_tot_fin$percentage) *
    (df_p_com_lam_tot_fin$percentage))/df_p_com_lam_tot_fin$tot_gen_week))
df_p_com_lam_tot_fin$CI[is.nan(df_p_com_lam_tot_fin$CI)] <- 0
df_p_com_lam_tot_fin$est <- "No"
df_p_com_lam_tot_fin$est[df_p_com_lam_tot_fin$tot_gen_week ==
    0] <- "Yes"
df_p_com_lam_tot_fin$Region <- "all"
f_df_list_lam_tot <- list()
for (i in levels(as.factor(df_p_com_lam_tot_fin$Region))) {
    for (w in levels(as.factor(df_p_com_lam_tot_fin$Lambda))) {
        if (w != "no") {
            df_voc_1 <- df_p_com_lam_tot_fin %>%
                dplyr::filter(Lambda == w & Region == i)
            for (z in 2:nrow(df_voc_1)) {
                if (df_voc_1$est[z] == "No") {
                  break
                }
            }
            for (x in z:nrow(df_voc_1)) {
                if (df_voc_1$est[x] == "Yes") {
                  one <- df_voc_1$percentage[x - 1]
                  for (y in x:nrow(df_voc_1)) {
                    if (df_voc_1$est[y] == "No") {
                      break
                    }
                  }
                  sec <- df_voc_1$percentage[y]
                  df_voc_1$percentage[x] <- (one + sec)/2
                }
            }
            weeks_index <- df_p_week %>%
                dplyr::select(week, collection_date) %>%
                unique() %>%
                tidyr::complete(collection_date)
            weeks_index <- weeks_index %>%
                dplyr::mutate(week = 1:nrow(weeks_index))
            df_voc <- df_voc_1[z:nrow(df_voc_1), ] %>%
                dplyr::left_join(weeks_index) %>%
                dplyr::arrange(week)
            fit <- loess(percentage ~ week, degree = 1, span = 21/nrow(df_voc),
                data = df_voc)
            fit_values <- predict(fit, data.frame(week = z:nrow(df_voc_1)),
                se = T)
            loessfit <- data.frame(week = z:nrow(df_voc_1), fitted = fit_values[["fit"]],
                se_fitted = fit_values[["se.fit"]])
            f_df_list_lam_tot[[paste(i, w, sep = "_")]] <- Reduce(function(x,
                y) dplyr::full_join(x, y, all = TRUE), list(df_voc,
                loessfit)) %>%
                mutate(collection_date = as.Date(collection_date))
            print(paste(i, w, sep = "_"))
        }
    }
}
f_df_lam_tot <- do.call(rbind, f_df_list_lam_tot) %>%
    mutate(fitted = pmax(0, fitted))
# plot of relative prevalence by week
rel_prev_lam_tot <- ggplot(subset(f_df_lam_tot), aes(x = as.POSIXct(collection_date))) +
    geom_point(data = subset(f_df_lam_tot, est == "No" & percentage >
        0), aes(y = percentage), color = "forestgreen", show.legend = F) +
    geom_errorbar(data = subset(f_df_lam_tot, est == "No" & percentage >
        0), aes(ymin = pmax(0, percentage - CI), ymax = pmin(1,
        percentage + CI)), color = "lightgrey") + geom_line(aes(y = fitted),
    color = "forestgreen") + geom_ribbon(aes(ymin = pmax(0, fitted -
    se_fitted), ymax = pmin(1, fitted + se_fitted)), alpha = 0.5,
    fill = "forestgreen") + ylim(c(0, 1)) + ylab("Prevalence") +
    scale_x_datetime(date_breaks = "2 month", date_labels = "%Y-%b",
        limits = c(as.POSIXct("2020-10-01"), as.POSIXct("2022-01-31"))) +
    theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90,
        hjust = 0, vjust = 0.5), strip.background = element_rect(fill = "white"))
ggsave("Relprev_Lambda_tot.jpg", plot = rel_prev_lam_tot, dpi = 300,
    height = 2, width = 3)
# 6:Estimation of cases by sublineage####
cas_city_cl <- cas_city %>%
    dplyr::select(City = DEPARTAMENTO, Date = FECHA_RESULTADO) %>%
    tidyr::drop_na() %>%
    dplyr::mutate(month = substr(as.character(Date), start = 5,
        stop = 6), year = substr(as.character(Date), start = 1,
        stop = 4), day = substr(as.character(Date), start = 7,
        stop = 8)) %>%
    dplyr::mutate(date = as.Date(paste(year, month, day, sep = "-")),
        City = gsub(" ", "", City))
regions_cases <- reg %>%
    dplyr::mutate(City = toupper(City))
cas_city_week <- cas_city_cl %>%
    dplyr::filter(date >= as.Date("2020-08-01")) %>%
    dplyr::filter(date <= as.Date("2022-01-31")) %>%
    dplyr::mutate(week = cut.Date(date, breaks = "1 week", labels = FALSE),
        collection_date = cut.Date(date, breaks = "1 week"))
cas_city_grouped <- cas_city_week %>%
    dplyr::filter(City != "") %>%
    dplyr::left_join(regions_cases) %>%
    dplyr::group_by(Region, collection_date) %>%
    dplyr::summarise(cases = n()) %>%
    tidyr::drop_na() %>%
    dplyr::ungroup() %>%
    tidyr::complete(collection_date, nesting(Region), fill = list(count = 0)) %>%
    dplyr::mutate(collection_date = as.Date(collection_date))
cases_reg_week <- f_df %>%
    dplyr::ungroup() %>%
    dplyr::select(week, collection_date) %>%
    unique() %>%
    right_join(cas_city_grouped)
cases_reg_fit_list <- list()
for (i in levels(as.factor(cases_reg_week$Region))) {
    cases_reg_week_fin <- cases_reg_week %>%
        dplyr::filter(Region == i)
    fit <- loess(cases ~ week, degree = 1, span = 21/nrow(cases_reg_week_fin),
        data = cases_reg_week_fin)
    fit_values <- predict(fit, data.frame(week = 1:length(fit$fitted)),
        se = T)
    loessfit <- data.frame(week = 1:length(fit$fitted), smooth_cases = fit_values[["fit"]],
        se_smooth_cases = fit_values[["se.fit"]])
    cases_reg_fit_list[[i]] <- Reduce(function(x, y) dplyr::full_join(x,
        y, all = TRUE), list(cases_reg_week_fin, loessfit)) %>%
        mutate(collection_date = as.Date(collection_date))
}
cases_reg_fit <- do.call(rbind, cases_reg_fit_list) %>%
    mutate(smooth_cases = pmax(0, smooth_cases))
f_df_cas_all <- left_join(f_df, cases_reg_fit) %>%
    mutate(smooth_lin_cases = fitted * smooth_cases, lin_cases = fitted *
        cases)
cas_Ls <- ggplot(subset(f_df_cas_all), aes(x = as.POSIXct(collection_date))) +
    geom_line(aes(y = smooth_lin_cases, color = SubLs)) + scale_color_manual(values = collin[names(collin) !=
    "Lambda"]) + facet_wrap(~Region) + ylab("Number of Cases") +
    scale_x_datetime(date_breaks = "1 month", date_labels = "%Y-%b",
        limits = c(as.POSIXct("2020-10-01"), as.POSIXct("2022-01-31"))) +
    theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90,
        hjust = 0, vjust = 0.5), legend.position = "bottom",
        legend.title = element_blank(), strip.background = element_rect(fill = "white"))
ggsave("Cases_Lambda_sublineages_reg.jpg", plot = cas_Ls, dpi = 300,
    height = 5, width = 7.5)
# 7:Estimation of total Lambda cases####
f_df_cas_lam_tot <- f_df_cas_all %>%
    dplyr::group_by(week, collection_date) %>%
    dplyr::summarise(totcas = sum(smooth_lin_cases))
cas_tot_fin <- ggplot(subset(f_df_cas_lam_tot), aes(x = as.POSIXct(collection_date))) +
    geom_line(aes(y = totcas), color = "forestgreen") + geom_line(data = f_df_cas_lam_tot,
    aes(y = totcas), color = "forestgreen") + ylab("Cases") +
    scale_x_datetime(date_breaks = "2 month", date_labels = "%Y-%b",
        limits = c(as.POSIXct("2020-10-01"), as.POSIXct("2022-01-31"))) +
    theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90,
        hjust = 0, vjust = 0.5), legend.position = "bottom",
        legend.title = element_blank())
ggsave("Cases_Lambda_tot.jpg", plot = cas_tot_fin, dpi = 300,
    height = 2, width = 3)
# 8:Save tables with data of prevalence and estimated cases
# by sublineage by region####
write.table(f_df_cas_all, "subLambda_epid_data.tsv", sep = "\t",
    quote = F, row.names = F)
write.table(df_p_reg, "subLambda_geno_data.tsv", sep = "\t",
    quote = F, row.names = F)
