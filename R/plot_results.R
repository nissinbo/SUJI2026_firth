source("R/project_helpers.R", local = TRUE)
set_project_locale()
library(ggplot2)

tabdir <- "results/tables"
figdir <- "results/figures"
dir.create(figdir,showWarnings=FALSE)
pairs <- read.csv(file.path(tabdir,"pair_summary.csv"),stringsAsFactors=FALSE)
scenarios <- read.csv("data/scenarios.csv",stringsAsFactors=FALSE)
covdesign <- read.csv("data/covariate_scenarios.csv",stringsAsFactors=FALSE)
exact <- merge(pairs,scenarios,by="scenario_id")
cov <- merge(pairs,covdesign,by="scenario_id")
comparisons <- c("default_logistf_sas","default_brglm2_sas","aligned_PL_logistf_sas")
labels <- c("PROC LOGISTIC Wald 対 logistf PL","PROC LOGISTIC Wald 対 brglm2 Wald","PROC LOGISTIC PL 対 logistf PL")
palette <- c("#E87722","#007C91","#4F775A")
theme_set(theme_minimal(base_size=17,base_family="Yu Gothic")+
 theme(panel.grid.minor=element_blank(),legend.position="bottom",legend.title=element_blank(),
       axis.text=element_text(colour="#23373B"),plot.title=element_text(face="bold"),plot.margin=margin(10,15,10,10)))
save_plot <- function(p,name,height=5) ggsave(file.path(figdir,paste0(name,".png")),p,width=11,height=height,dpi=180,bg="white")
z <- exact[exact$cohort=="all_outputs" & exact$comparison%in%comparisons,]
z$comparison<-factor(z$comparison,comparisons,labels);z$allocation<-paste0("曝露：非曝露 = 1:",z$ratio)
save_plot(ggplot(z,aes(n,100*ci_difference_conditional_probability,colour=comparison))+
 geom_line(linewidth=1)+geom_point(size=3)+facet_wrap(~allocation)+scale_colour_manual(values=palette)+
 guides(colour=guide_legend(ncol=1))+
 scale_x_continuous(breaks=c(20,50,100,200))+labs(x="総症例数",y="信頼区間の判定差の割合（%）"),"exact_ci_decision")

z <- cov[cov$cohort=="all_outputs" & cov$comparison%in%comparisons,]
z$comparison<-factor(z$comparison,comparisons,labels)
z$background<-factor(ifelse(z$n==40,"40例・20:20","100例・20:80"),c("40例・20:20","100例・20:80"))
z$x<-match(z$rho,c(0,.5,.9,.99))
save_plot(ggplot(z,aes(x,100*ci_difference_conditional_probability,colour=comparison,group=comparison))+
 geom_line(linewidth=.9)+geom_errorbar(aes(ymin=100*pmax(0,ci_difference_conditional_probability-1.96*ci_difference_mcse),
 ymax=100*pmin(1,ci_difference_conditional_probability+1.96*ci_difference_mcse),
 linetype="エラーバー：Monte Carlo誤差に基づく95%信頼区間"),width=.06,position=position_dodge(.12),
 show.legend=c(colour=FALSE,linetype=TRUE))+
 geom_point(size=2.8,position=position_dodge(.12))+facet_wrap(~background)+scale_colour_manual(values=palette)+
 scale_linetype_manual(values=1)+
 guides(colour=guide_legend(ncol=1,order=1),
        linetype=guide_legend(order=2,override.aes=list(colour="#23373B")))+
 theme(legend.box="vertical",legend.box.just="left",legend.spacing.y=grid::unit(2,"pt"))+
 scale_x_continuous(breaks=1:4,labels=c("0","0.5","0.9","0.99"))+
 labs(x="連続共変量間の相関",y="信頼区間の判定差\n（%）"),"covariate_ci_decision")


# Shared display order and wording for confidence interval and p-value figures.
# Keep the original set of contrasts; %FL's combined-control comparison is separate.
selected_factors <- c(
  "LF_inference",
  "SAS_algorithm_default", "LF_algorithm_default", "LF_algorithm_tight",
  "LF_maxit", "BG_maxit",
  "SAS_gradient", "LF_threshold", "BG_threshold",
  "SAS_PL_control", "LF_PL_control"
)
factor_implementations <- c(
  "logistf", "PROC LOGISTIC", "logistf", "logistf", "logistf", "brglm2",
  "PROC LOGISTIC", "logistf", "brglm2", "PROC LOGISTIC", "logistf"
)
factor_categories <- c("区間の方法", rep("計算アルゴリズム", 3),
  rep("反復回数", 2), rep("点推定の収束基準", 3), rep("区間探索の収束基準", 2))
factor_actions <- c(
  "Wald区間 → PL区間",
  "Fisher scoring → Newton法",
  "Newton法 → IRLS",
  "Newton法 → IRLS（厳しい設定）",
  "反復上限を増やす",
  "反復上限を増やす",
  "基準を厳しく（Wald）",
  "基準を厳しく",
  "基準を厳しく",
  "PL区間の探索条件を厳しく",
  "PL区間の探索条件を厳しく"
)
factor_labels <- paste0(factor_implementations, "：", factor_actions)

factor_difference_plot <- function(z, metric = c("ci", "p")) {
  metric <- match.arg(metric)
  stopifnot(nrow(z) == length(selected_factors), !anyDuplicated(z$comparison),
            setequal(z$comparison, selected_factors))
  i <- match(z$comparison, selected_factors)
  categories <- factor_categories
  if (metric == "p") categories[categories == "区間の方法"] <- "検定の方法"
  z$category <- factor(categories[i], unique(categories))
  z$change <- factor(z$comparison, rev(selected_factors))
  z$median <- pmax(z[[paste0(metric, "_delta_median")]], 1e-12)
  z$p99 <- pmax(z[[paste0(metric, "_delta_p99")]], 1e-12)
  action_labels <- setNames(factor_labels, selected_factors)
  if (metric == "p") action_labels["LF_inference"] <- "logistf：Wald検定 → PLRT"
  top_tick <- if (metric == "ci") 1 else 0.1
  ggplot(z, aes(y = change)) +
    geom_segment(aes(x = median, xend = p99, yend = change), linewidth = 2, colour = "#B4C3C6") +
    geom_point(aes(x = median, colour = "中央値"), size = 3) +
    geom_point(aes(x = p99, colour = "99パーセンタイル"), shape = 15, size = 3) +
    facet_grid(category ~ ., scales = "free_y", space = "free_y", switch = "y",
      labeller = labeller(category = c("区間の方法" = "区間の方法", "検定の方法" = "検定の方法",
        "計算アルゴリズム" = "計算\nアルゴリズム", "反復回数" = "反復回数",
        "点推定の収束基準" = "点推定の\n収束基準", "区間探索の収束基準" = "区間探索の\n収束基準"))) +
    scale_y_discrete(labels = action_labels) +
    scale_x_log10(breaks = c(1e-12, 1e-9, 1e-6, 1e-3, top_tick),
      labels = c("≤10⁻¹²", "10⁻⁹", "10⁻⁶", "10⁻³", as.character(top_tick))) +
    scale_colour_manual(values = c("中央値" = "#007C91", "99パーセンタイル" = "#E87722"),
      breaks = c("中央値", "99パーセンタイル")) +
    labs(x = if (metric == "ci") "信頼区間の数値差（対数オッズ比）" else "p値の絶対差", y = NULL) +
    theme(strip.placement = "outside", strip.background = element_blank(),
      strip.text.y.left = element_text(angle = 0, size = 14),
      panel.spacing.y = grid::unit(9, "pt"), axis.text.y = element_text(size = 14),
      legend.spacing.x = grid::unit(10, "pt"))
}

z<-exact[exact$scenario_id=="E050R1" & exact$cohort=="common_inference_available" & exact$comparison%in%selected_factors,]
save_plot(factor_difference_plot(z, "ci"), "factor_ci_difference_full", 5.5)

# Main slide: medians for the four principal comparisons.
main_factors <- c("LF_inference", "SAS_gradient", "SAS_PL_control", "BG_threshold")
main_labels <- c("区間の方法（logistf）\nWald → PL",
  "点推定の収束基準\nPROC LOGISTIC（Wald）",
  "区間探索の収束基準\nPROC LOGISTIC",
  "点推定の収束基準\nbrglm2")
main <- z[match(main_factors, z$comparison), ]
other <- z[!z$comparison %in% main_factors, ]
stopifnot(nrow(main) == 4L, !anyNA(main$comparison), nrow(other) == 7L,
  all(other$ci_delta_p99 < 1e-11))
main$value <- main$ci_delta_median
main$change <- factor(main$comparison, rev(main_factors), rev(main_labels))
main$highlight <- main$comparison == "LF_inference"
# Direct labels retain three significant figures without long strings of zeros.
main$label <- vapply(main$value, function(v) {
  if (v >= .01) return(format(signif(v, 3), scientific = FALSE, trim = TRUE))
  exponent <- floor(log10(v))
  sprintf('%.3g %%*%% 10^{%d}', v / 10^exponent, exponent)
}, character(1))
save_plot(ggplot(main, aes(value, change, colour = highlight)) +
  geom_point(size = 3.5) +
  geom_text(aes(label = label), parse = TRUE, vjust = -1.1, size = 4.1,
    show.legend = FALSE) +
  scale_colour_manual(values = c("FALSE" = "#56747A", "TRUE" = "#E87722"), guide = "none") +
  scale_x_log10(limits = c(1e-10, 10), breaks = c(1e-9, 1e-6, 1e-3, 1),
    labels = c("10⁻⁹", "10⁻⁶", "0.001", "1")) +
  scale_y_discrete(expand = expansion(add = c(.45, .8))) +
  labs(x = "信頼区間の数値差の中央値（対数オッズ比）", y = NULL) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_line(colour = "#E5EBEC"),
    axis.text.y = element_text(size = 14, lineheight = 1.05)),
  "factor_ci_difference", 3.8)
save_plot(factor_difference_plot(z, "p"), "factor_p_difference", 5.5)
