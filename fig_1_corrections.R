library(tidyverse) 
library(estimatr)
library(magrittr)
library(emmeans)
library(broom)
library(metafor)

t2 <- "https://github.com/thomasjwood/corrections_eeurope/raw/main/data/t2.rds" %>% 
  url %>% 
  readRDS

t3 <- t2 %>% 
  group_by(
    country, trial
  ) %>% 
  nest %>% 
  mutate(
    mods = data %>% 
      map(
        \(i)
        lm_robust(
          ans_num ~ cond, 
          data = i
        )
      )
  )

t3$emms <- t3$mods %>% 
  map(
    \(i)
    i %>% 
      emmeans(
        revpairwise ~ cond
      )
  )

t4 <- t3 %>% 
  select(
    country, trial, emms
  ) %>% 
  pmap(
    \(country, trial, emms)

    emms %>% 
      extract2("emmeans") %>% 
      tidy(conf.int = T) %>% 
      mutate(
        country, trial
      )
  ) %>% 
  list_rbind %>% 
  mutate(
    colfac = "Fitted"
  ) %>% 
  bind_rows(
    t3 %>% 
      select(
        country, trial, emms
      ) %>% 
      pmap(
        \(country, trial, emms)
        
        emms %>% 
          extract2("contrasts") %>% 
          tidy(conf.int = T) %>% 
          mutate(
            country, trial
          )
      ) %>% 
      list_rbind %>% 
      mutate(
        colfac = "Contrasts"
      ) %>% 
      rename(
        cond = contrast
      ) %>% 
      select(-term)
  ) %>% 
  mutate(
    colfac = colfac %>% 
      factor(
        c("Fitted", "Contrasts")
      )
  )

# descriptive labels

t4$dep_lab <- case_when(
  t4$country %>% 
    equals("blr") &
    t4$trial %>% 
    equals("fc_1") ~ "Ukraine is not a real country, belongs to Russia",
  t4$country %>% 
    equals("blr") &
    t4$trial %>% 
    equals("fc_2") ~ "Bucha war crime videos are fake",
  t4$country %>% 
    equals("blr") &
    t4$trial %>% 
    equals("fc_3") ~ "Ukrainian textbooks endorse Nazi ideology",
  t4$country %>% 
    equals("est") &
    t4$trial %>% 
    equals("fc_1") ~ "Estonians defend killer of Daria Dugin",
  t4$country %>% 
    equals("est") &
    t4$trial %>% 
    equals("fc_2") ~ "Estonians discriminate against Russians in Narva",
  t4$country %>% 
    equals("est") &
    t4$trial %>% 
    equals("fc_3") ~ "Estonia sabotaged Nord Stream pipeline",
  t4$country %>% 
    equals("kaz") &
    t4$trial %>% 
    equals("fc_1") ~ "Ukrainian refugees create German teacher shortage",
  t4$country %>% 
    equals("kaz") &
    t4$trial %>% 
    equals("fc_2") ~ "Finland recognizes Donetsk and Luhanks People's Republics",
  t4$country %>% 
    equals("kaz") &
    t4$trial %>% 
    equals("fc_3") ~ "Ukraine creates fake videos of Russian war crimes",
  
  t4$country %>% 
    equals("rus") &
    t4$trial %>% 
    equals("fc_1") ~ "Russia not aggressor in the Ukraine war",
  t4$country %>% 
    equals("rus") &
    t4$trial %>% 
    equals("fc_2") ~ "Ukraine grain sales used to buy weapons",
  t4$country %>% 
    equals("rus") &
    t4$trial %>% 
    equals("fc_3") ~ "Russian military not involved in Bucha massacre",
  
  t4$country %>% 
    equals("ukr") &
    t4$trial %>% 
    equals("fc_1") ~ "Ukraine killed Russians civilians in Donbas",
  t4$country %>% 
    equals("ukr") &
    t4$trial %>% 
    equals("fc_2") ~ "Ukraine raises money to buy nuclear weapons",
  t4$country %>% 
    equals("ukr") &
    t4$trial %>% 
    equals("fc_3") ~ "US builds biological weapons in Ukraine"
  )

t4$dep_lab %<>%
  factor(
    t4 %>% 
      filter(
        colfac %>% 
          equals("Fitted") %>% 
          not
      ) %>% 
      arrange(
        desc(estimate)
      ) %>% 
      use_series(dep_lab)
  )

# meta summaries

t5 <- t4 %>% 
  select(
    colfac, cond, estimate, std.error
  ) %>% 
  group_by(colfac, cond) %>% 
  nest %>% 
  mutate(
    country = "Overall"
  ) %>% 
  # binding for country specific meta analytic summaries
  bind_rows(
    t4 %>% 
      select(
        colfac, cond, estimate, std.error, country
      ) %>% 
      group_by(colfac, cond, country) %>% 
      nest
  )


t5$mods <- t5$data %>% 
  map(
    \(i)
    rma(
      yi = i$estimate,
      sei = i$std.error
    )
  )

t6 <- t4 %>% 
  bind_rows(
    t5 %>% 
      pmap(
        \(colfac, cond, data, mods, country)
        
        mods %>% 
          tidy(
            conf.int = T, 
            conf.level = 0.95
          ) %>% 
          mutate(
            colfac, cond, country
          )
      ) %>% 
      list_rbind %>% 
      mutate(
        dep_lab = "Meta analytic summary",
      ) %>% 
      select(-term, -type)
  ) %>% 
  mutate(
    colfac = colfac %>% 
      case_match(
        "Fitted" ~ "Conditional means\n(5pt scale of agreement)",
        "Contrasts" ~ "Correction effects"
      ),
    country = country %>% 
      case_match(
        "blr" ~ "Belarus",
        "est" ~ "Estonia",
        "kaz" ~ "Kazakhstan",
        "rus" ~ "Russia",
        "ukr" ~ "Ukraine",
        .default = country 
      )
  )

t6$country %<>%
  factor(
    t6 %>% 
      filter(
        colfac %>%
          str_detect("Conditional ")
      ) %>% 
      group_by(
        country
      ) %>% 
      summarize(
        mu = estimate %>% mean
      ) %>% 
      arrange(mu) %>% 
      use_series(country) %>% 
      str_remove("Overall") %>% 
      stringi::stri_remove_empty() %>% 
      c("Overall")
  )

t6$dep_lab <- t6$dep_lab %>% 
  factor(
    t6$dep_lab %>% 
      fct_reorder(
        t6$estimate, .desc = T
      ) %>% 
      levels %>% 
      str_remove("Meta analytic summary") %>% 
      stringi::stri_remove_empty() %>% 
      c("Meta analytic summary")
  )

t6$lab <- t6$estimate %>% 
  round(2) %>% 
  as.character %>% 
  str_replace_all(
    "0\\.", "."
  )


t6$lab %<>% 
  equals("0") %>% 
  ifelse(
    t6$estimate %>% 
      round(3) %>% 
      as.character %>% 
      str_replace_all(
        "0\\.", "."
      ),
    t6$lab
  )

t6$lab %<>% 
  str_detect("-") %>% 
  ifelse(
    t6$lab %>% 
      str_pad(width = 4, side = "right", pad = "0"),
    t6$lab %>% 
      str_pad(width = 3, side = "right", pad = "0")
  )

t6$lab <- t6$colfac  %>%  
  str_detect(" effect") %>% 
  ifelse(
    str_c(
      t6$lab,
      t6$p.value %>%
        gtools::stars.pval() %>% 
        str_remove("\\.") %>% 
        str_remove("\\s")
    ),
    t6$lab
  )


t6 %>% 
  ggplot() +
  geom_pointrange(
    aes(
      xmin = conf.low,
      xmax = conf.high,
      x = estimate,
      y = dep_lab, 
      color = cond,
      group =  cond,
      shape = dep_lab %>% 
        str_detect(" summary")
      ),
    fill = "grey95",
    size = .3,
    data = t6 %>% 
      filter(
        colfac == "Conditional means\n(5pt scale of agreement)"
      ),
    position = position_dodge(width = .2)
  )  +
  geom_text(
    aes(
      x = estimate,
      y = dep_lab, 
      color = cond,
      label = lab
    ),
    nudge_y = -.35,
    size = 2.7,
    data = t6 %>% 
      filter(
        colfac == "Conditional means\n(5pt scale of agreement)" &
          cond == "control"
      )
    )  +
  geom_vline(
    aes(
      xintercept = xint
    ),
    linewidth = .4,
    data = tribble(
      ~xint, ~colfac,
      0, "Correction effects"
    ) %>% 
      mutate(
        colfac  = colfac %>% 
          factor(
            c("Conditional means\n(5pt scale of agreement)", "Correction effects")
          )
      ),
    linetype = "dashed"
  ) +
  geom_label(
    aes(
      x = estimate,
      y = dep_lab,
      label = lab
    ),
    size = 3,
    fill = "grey95",
    nudge_y = .2,
    data = t6 %>% 
      filter(
        colfac != "Conditional means\n(5pt scale of agreement)" &
          country != "Overall"
      ), 
    label.size = 0
  )  +
  geom_text(
    aes(
      x = estimate,
      y = dep_lab, 
      color = cond,
      label = lab
    ),
    nudge_y = .3,
    size = 2.7,
    data = t6 %>% 
      filter(
        colfac == "Conditional means\n(5pt scale of agreement)" &
          cond != "control"
      ),
  )  +
  geom_text(
    aes(
      x, y, label = label, color = cond
    ),
    fontface = "bold.italic",
    size = 3,
    data = tribble(
      ~x, ~y, ~label, ~cond,
      1.35, 3.8, "Treatment", "treatment",
      2.3, 3.7, "Control", "control"
    ) %>%
      mutate(
        colfac = "Conditional means\n(5pt scale of agreement)" %>%
          factor(
            t6$colfac %>%
              factor %>%
              levels
          ),
        country = "Ukraine" %>%
          factor(
            t6$country %>% 
              levels
          ),
        cond = cond %>% 
          factor
      )
  ) +
  geom_linerange(
    aes(
      xmin = conf.low,
      xmax = conf.high,
      x = estimate,
      y = dep_lab, 
      # color = cond
    ),
    data = t6 %>% 
      filter(
        colfac != "Conditional means\n(5pt scale of agreement)"
      ), 
    linewidth = .3
  )  +
  geom_point(
    aes(
      x = estimate, 
      y = dep_lab,
      shape = dep_lab %>% 
        str_detect(" summary"),
      size = country %>% 
        str_detect("Overall")
    ),
    # size = .,
    data = t6 %>% 
      filter(
        colfac != "Conditional means\n(5pt scale of agreement)"
      ),
    fill = "white"
  ) + 
  geom_label(
    aes(
      x = estimate,
      y = dep_lab,
      label = lab
      # color = cond
    ),
    size = 3,
    fill = "transparent",
    nudge_y = .2,
    data = t6 %>% 
      filter(
        colfac != "Conditional means\n(5pt scale of agreement)" &
          dep_lab  %>% 
          str_detect(" summary") &
          country %>% 
          str_detect("Overall") %>% 
          not
      ), 
    family = "Roboto",
    label.size = 0
  )  +
  geom_label(
    aes(
      x = estimate,
      y = dep_lab,
      label = lab
    ),
    size = 3,
    fill = "transparent",
    data = t6 %>% 
      filter(
        colfac != "Conditional means\n(5pt scale of agreement)" &
          dep_lab  %>% 
          str_detect(" summary") &
          country %>% 
          str_detect("Overall") 
      ), 
    label.size = 0
  )  +
  geom_segment(
    aes(x, y, xend = xend, yend = yend),
    linewidth = .5,
    arrow = arrow(
      length = unit(0.15,"cm"), 
      type = "closed"),
    data = tribble(
      ~x, ~y, ~xend, ~yend,
      -.05, .5,  -.4, .5,
      .05, .5,  .4, .5
    ) %>% 
      mutate(
        colfac = "Correction effects" %>%
          factor(
            t6$colfac %>%
              factor %>%
              levels
          ),
        country = "Ukraine" %>%
          factor(
            t6$country %>% 
              levels
          )
      )
  ) +
  geom_text(
    aes(
      x, y, label =  lab
    ),
    size = 2.5,
    hjust = 0,
    fontface = "italic",
    data = tribble(
      ~x, ~y, ~lab,
      -.6, .8, "More accurate",
      .1, .8, "Less accurate",
    ) %>% 
      mutate(
        colfac = "Correction effects" %>%
          factor(
            t6$colfac %>%
              factor %>%
              levels
          ),
        country = "Ukraine" %>%
          factor(
            t6$country %>% 
              levels
          )
      )
  ) +
  facet_grid(
    country ~ colfac,
    scales = "free",
    space = "free_y"
  ) +
  scale_y_discrete(
    breaks = c(
      "Ukraine raises money to buy nuclear weapons", 
      "Ukraine killed Russians civilians in Donbas", 
      "US builds biological weapons in Ukraine",
      "Estonians discriminate against Russians in Narva", 
      "Estonia sabotaged Nord Stream pipeline",
      "Ukraine is not a real country, belongs to Russia", 
      "Estonians defend killer of Daria Dugin",
      "Finland recognizes Donetsk and Luhanks People's Republics", 
      "Meta analytic summary",
      "Ukrainian textbooks endorse Nazi ideology", 
      "Bucha war crime videos are fake",
      "Ukrainian refugees create German teacher shortage", 
      "Ukraine creates fake videos of Russian war crimes",
      "Russia not aggressor in the Ukraine war", 
      "Ukraine grain sales used to buy weapons",
      "Russian military not involved in Bucha massacre"
    ),
    labels = c(
      "Ukraine raises money\nto buy nuclear weapons", 
      "Ukraine killed Russians\ncivilians in Donbas", 
      "US builds biological\nweapons in Ukraine",
      "Estonians discriminate\nagainst Russians in Narva", 
      "Estonia sabotaged\nNord Stream pipeline",
      "Ukraine is not a real\ncountry, belongs to Russia", 
      "Estonians defend\nkiller of Daria Dugin",
      "Finland recognizes Donetsk\nand Luhanks People's Republics", 
      "Meta analytic summary",
      "Ukrainian textbooks\nendorse Nazi ideology", 
      "Bucha war crime\nvideos are fake",
      "Ukrainian refugees create\nGerman teacher shortage", 
      "Ukraine creates fake videos\nof Russian war crimes",
      "Russia not aggressor\nin the Ukraine war", 
      "Ukraine grain sales\nused to buy weapons",
      "Russian military not\ninvolved in Bucha massacre"
    ), 
    # expand = expansion(add = c(.1, .1)) 
  ) +
  scale_shape_manual(
    values = c(19, 23), 
    guide = "none"
  ) +
  scale_size_manual(
    values = c(1.5, 9), 
    guide = "none"
  ) +
  scale_color_brewer(type = "qual", palette = 2) +
  labs(
    x = "",
    y =  ""
  ) 
