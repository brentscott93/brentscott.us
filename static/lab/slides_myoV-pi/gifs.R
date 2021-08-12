library(gganimate)
library(tidyverse)
library(transformr)
wt_trace_data <- read_csv("~/lasertrapr/project_myoV-phosphate/myoV-WT_pH-7.0_0mM-Pi/2020-06-25/obs-14/trap-data.csv")
wt_trace_events <- read_csv("~/lasertrapr/project_myoV-phosphate/myoV-WT_pH-7.0_0mM-Pi/2020-06-25/obs-14/measured-events.csv")
wt_trace <- data.frame(index = (1:nrow(wt_trace_data)/5000),
                       raw = wt_trace_data$processed_bead,
                       model = wt_trace_data$hm_overlay,
                       id = "WT 0mM-Pi") %>% 
  slice(1:25500) %>% 
  mutate(time = 1:nrow(.)/5000)

times <- seq(100, 25500, by = 100)

state_id <- list()
for(row in 1:nrow(wt_trace_events)){
  state_id[[row]] <- wt_trace_events$cp_event_start_dp[[row]]:wt_trace_events$cp_event_stop_dp[[row]]
}
state_id <- unlist(state_id)

state <- tibble(seq = full_seq(c(1, state_id), 1),
                state = ifelse(seq %in% state_id, "event", "baseline"),
                index = seq/5000)
  

overlay <- data.frame(index = (1:nrow(wt_trace_data)/5000),
                      raw = wt_trace_data$processed_bead,
                      model = wt_trace_data$hm_overlay,
                      id = "WT 0mM-Pi") %>% 
  inner_join(state)
gif_color <- RColorBrewer::brewer.pal( 8, "Dark2")[c(8, 1)]
for(i in seq_along(times)){
  overlay2 <- overlay %>% 
  slice(1:times[[i]]) %>% 
  mutate(time = 1:nrow(.)/5000)

plotter <- 
  ggplot()+
  geom_point(data = wt_trace, aes(index, raw), size = 0.5, alpha = 0.5, shape = 16)+
  geom_point(data = overlay2, aes(index, raw, color = state),size = 0.5, alpha = 0.5, shape = 16)+
    scale_color_manual(values = gif_color)+
  theme_void()+
    theme(legend.position = "none")

ggsave(paste0('images/analysis2/plot', str_pad(i, pad = 0,width = 3 , "left"), '.png'), plot = plotter )

}

########ensemble avearge animation
library(cowplot)
library(tidyverse)
library(data.table)
gif_color <- RColorBrewer::brewer.pal( 8, "Dark2")[c(8, 1)]
ee <- list.files("~/lasertrapr/project_myoV-phosphate/myoV-WT_pH-7.0_0mM-Pi",
                 pattern = 'ensemble-data.csv',
                 recursive = T,
                 full.names = T) %>% 
  map(fread) %>% 
  rbindlist() %>% 
  filter(direction == 'forward') %>% 
  mutate(id = ifelse(ensemble_index < 0 , "base", "event"))

ee_unique <- ee %>% 
  select(conditions, date, obs, event_index) %>% 
  unique %>% 
  rownames_to_column() %>% 
  mutate(rowname = as.numeric(rowname)) %>% 
  inner_join(ee)

n_event <- max(ee_unique$rowname)
size_seq <- seq(1, by = 0.005, length.out = n_event)
alpha_seq <- seq(0.4, 0.8, length.out = n_event)
for(i in seq_len(n_event)){
  data <- 
    ee_unique %>% 
    filter(rowname <= i) %>% 
    group_by(ensemble_index, id) %>% 
    summarize(mean = mean(data)) %>% 
    ggplot()+
    geom_point(aes(ensemble_index, mean, color = id), size = size_seq[[i]],
               alpha = alpha_seq[[i]])+
    #geom_text(aes(ensemble_index, mean, label = ensemble_index), color = 'black')+
    coord_cartesian(c(-50, 100), 
                    c(-5, 15)
                    )+
    ggtitle(paste0("Average of ", i,  ' events'))+
    scale_color_manual(values = gif_color)+
    xlab('Time (datapoints)')+
    ylab('Displacement (nm)')+
    theme_cowplot()+
    theme(legend.position = "none")
  
 # ggsave('images/ee-figs/ee-horizontal3.png')
  ggsave(paste0('images/ee-figs/gif3/plot', str_pad(i, pad = 0,width = 3 , "left"), '.png'), 
         plot = data)
  
}

plot_grid(ee1, ee2, ncol = 2)

library(tidyverse)
library(lasertrapr)
f <- list(project=list(path="~/lasertrapr/project_myoV-phosphate"))

rv <- list()
rv$data <- lasertrapr:::summarize_trap_data(f = f,
                               hz = 5000,
                               factor_order = c("myoV-WT_pH-7.0_0mM-Pi",
                                                "myoV-WT_pH-7.0_30mM-Pi",
                                                "myoV-S217A_pH-7.0_0mM-Pi",
                                                "myoV-S217A_pH-7.0_30mM-Pi"),
                               is_shiny = FALSE)


plot_colors <- c(RColorBrewer::brewer.pal(8, "Blues")[c(5, 8)],
                 RColorBrewer::brewer.pal(8, "Reds")[c(5, 8)])

write_csv(rv$data$event_files_filtered, "data/event_files_filtered.csv")
write_csv(rv$data$summary, "data/summary.csv")
########### step size ##################
step_histo <- ggplot(data = event_files_filtered,
                     aes(x = displacement_nm,
                         fill = conditions))+
  geom_histogram(aes(y = stat(density)),
                 binwidth = 3,
                 color = "black",
                 size = 0.9)+
  facet_wrap(~conditions, ncol = 2)+
  xlab("Step Size (nm)")+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(breaks = seq(-100, 100, by = 10))+
  scale_fill_manual(values = plot_colors)+
  #scale_fill_brewer(palette = "Dark2")+
  theme_linedraw(base_size = 11)+
  theme(panel.grid = element_blank(),
        legend.position = 'none')

step_aov <-  event_files_filtered %>%
  rstatix::anova_test(displacement_nm ~ conditions)

step_tukey <-  event_files_filtered %>%
  stats::aov(displacement_nm ~ conditions, data = .) %>%
  rstatix::tukey_hsd()

step_is_sig <- which(step_tukey$p.adj < 0.05)

step_table <- step_tukey %>%
  dplyr::select(group1, group2, p.adj) %>%
  dplyr::mutate_if(is.numeric, round, digits = 3) %>%
  ggpubr::ggtexttable (theme = ggpubr::ttheme(base_style = 'light', base_size = 9), rows = NULL)

if(any(step_is_sig))  step_table %<>% ggpubr::table_cell_bg(row = step_is_sig + 1, column = 3, fill = '#ffa9a6')

step_tukey <- rstatix::add_xy_position(step_tukey, data = event_files_filtered, formula = displacement_nm ~ conditions)

step_violin <- ggpubr::ggviolin(event_files_filtered,
                                x = "conditions",
                                y = "displacement_nm",
                                fill = "conditions",
                                palette = plot_colors,
                                add = "boxplot",
                                xlab = '',
                                ylab = 'Displacement (nm)',
                                ggtheme = ggpubr::theme_pubr(base_size=11))+
  ggpubr::stat_pvalue_manual(step_tukey, bracket.nudge.y = seq(10, by = 4, length.out = nrow(step_tukey)))+
  scale_y_continuous(breaks = seq(-50, 1000, by = 15))+
  rotate_x_text(angle = 45)+
  # rotate()+
  #yscale('log10', .format = T)+
  labs(
    subtitle = rstatix::get_test_label(step_aov, detailed = TRUE),
    caption = rstatix::get_pwc_label(step_tukey)
  ) +
  theme(legend.position = 'none')



#  step_top <- cowplot::plot_grid(step_violin, step_histo, align = "hv", axis = "bt", rel_widths = c(1,1))
step_side <- cowplot::plot_grid(step_histo, step_table, ncol = 1, nrow = 2,  rel_heights = c(0.75,0.25))
step_main <-  cowplot::plot_grid(step_violin, step_side, ncol = 2, nrow = 1,  rel_heights = c(1, 1))
# now add the title
step_title <- cowplot::ggdraw() +
  cowplot::draw_label(
    "Displacement Distributions, ANOVA, & Tukey Table",
    fontface = 'bold',
    size = 18,
    x = 0.5
  )

cowplot::plot_grid(step_title, 
                   step_main, #step_table,
                   ncol = 1, 
                   nrow =2,
                   rel_heights = c(0.05, 1))

#############################################################
library(gganimate)
library(tidyverse)
library(transformr)
wt_trace_data <- read_csv("~/lasertrapr/project_myoV-phosphate/myoV-WT_pH-7.0_0mM-Pi/2020-06-25/obs-14/trap-data.csv")
wt_trace_events <- read_csv("~/lasertrapr/project_myoV-phosphate/myoV-WT_pH-7.0_0mM-Pi/2020-06-25/obs-14/measured-events.csv")
state_id <- list()
for(row in 1:nrow(wt_trace_events)){
  state_id[[row]] <- wt_trace_events$cp_event_start_dp[[row]]:wt_trace_events$cp_event_stop_dp[[row]]
}
state_id <- unlist(state_id)

state <- tibble(seq = full_seq(c(1, state_id), 1),
                state = ifelse(seq %in% state_id, "event", "baseline"),
                index = seq/5000)

gif_color <- RColorBrewer::brewer.pal( 8, "Dark2")[c(8, 1)]
start <- seq(1, 351000, by = 1000)
stop <- seq(25000, by = 1000, along.with = start)

for(i in seq_along(start)){
overlay <- data.frame(index = (1:nrow(wt_trace_data)/5000),
                      raw = wt_trace_data$processed_bead,
                      model = wt_trace_data$hm_overlay,
                      id = "WT 0mM-Pi") %>% 
  inner_join(state)  %>%
  slice(start[[i]]:stop[[i]])

  plotter <-
    ggplot()+
    geom_point(data = overlay, aes(index, raw, color = state), size = 0.5, alpha = 0.5, shape = 16)+
    #geom_point(data = overlay2, aes(index, raw, color = state),size = 0.5, alpha = 0.5, shape = 16)+
    scale_color_manual(values = gif_color)+
    theme_void()+
    theme(legend.position = "none")
  
  ggsave(paste0('images/analysis2/plot', str_pad(i, pad = 0,width = 3 , "left"), '.png'), plot = plotter )
  
}

####################################################################################
event_files_filtered <- read_csv("data/event_files_filtered.csv") %>%
  mutate(conditions2 = case_when(
    conditions == "myoV-WT_pH-7.0_0mM-Pi" ~ "WT 0mM-Pi",
    conditions == "myoV-WT_pH-7.0_30mM-Pi" ~"WT 30mM-Pi",
    conditions ==  "myoV-S217A_pH-7.0_0mM-Pi" ~ "S217A 0mM-Pi",
    conditions ==  "myoV-S217A_pH-7.0_30mM-Pi" ~ "S217A 30mM-Pi"
  ))

event_files_filtered$conditions2 <- factor(event_files_filtered$conditions2,
                                           levels = c("WT 0mM-Pi", "WT 30mM-Pi",
                                                      "S217A 0mM-Pi", "S217A 30mM-Pi"))

plot_colors <- c(RColorBrewer::brewer.pal(8, "Blues")[c(5, 8)],
                 RColorBrewer::brewer.pal(8, "Reds")[c(5, 8)])

data2 <- filter(event_files_filtered, conditions2 == 'WT 0mM-Pi')
for(i in 1:nrow(data2)){
 
p <-  data2 %>% 
  slice(1:i) %>% 
ggplot()+
    geom_dotplot(aes(x = displacement_nm, fill = conditions2,),
                                                  binwidth = 1)+
                                     # geom_richtext(data = summary_data,
                                     #               aes(x = 3ing = grid::unit(rep(0, 4), "pt"))+
                                     # geom_richtext(data = summary_data,
                                     #               aes(x = -20, y =0.9, label = n_events),
                                     #               fill = NA,5, y =0.9, label = mean_label),
                                     #               fill = NA, label.color = NA, size = 6, # remove background and outline
                                     #               label.padd label.color = NA, size = 6, # remove background and outline
                                     #               label.padding = grid::unit(rep(0, 4), "pt"))+
                                     #facet_wrap(~data_type, ncol = 2, scales = "free_y")+
      xlab("nanometers")+
                                     ylab("")+
                                     ggtitle('Displacements')+
                                     scale_y_continuous(expand = c(0,0))+
                                     scale_x_continuous(breaks = seq(-100, 100, by = 10))+
                                     scale_fill_manual(values = gif_color[[2]])+
                                     coord_cartesian(c(-20, 30), c(0, 1))+
                                     #scale_fill_brewer(palette = "Dark2")+
                                     theme_linedraw(base_size = 15)+
                                     theme(panel.grid = element_blank(),
                                           legend.position = 'none',
                                           #axis.line.y = element_blank(),
                                           axis.text.y = element_blank(),
                                           axis.ticks.y = element_blank(),
                                           #strip.background = element_rect("white"),
                                           strip.text = element_text(face = "bold"))
                                     
          ggsave(paste0('images/ton-dotplot/plot', str_pad(i, pad = 0,width = 3 , "left"), '.png'), plot = p )
          }



<video width="400" height="300" controls style = "position: fixed; bottom: 0; right: 0;">
  <source src="images/ee-figs/gif/ee.gif" type = "video/mov">
  Your browser does not support the video tag.
</video>  
  