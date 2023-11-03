library(tidyverse) 
library(estimatr)
library(magrittr)
library(estimatr)
library(emmeans)
library(broom)
library(metafor)
library(ggtext)


t3 <- "https://github.com/thomasjwood/corrections_eeurope/raw/main/data/t3.rds" %>% 
  url %>% 
  readRDS %>% 
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
    
    # emms <- t3$emms[[1]]
    
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
        
        # emms <- t3$emms[[1]]
        
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
  # estonia
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
  # kazahstan
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
  
  # russia
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
  
  # ukraine
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


# checking how country omissions effect our results


cl <- t5 %>% 
  filter(
    colfac %>% str_detect("Contrasts")
  ) %>% 
  select(-mods) %>% 
  use_series(country)

t7 <- tibble(
  country = cl
  ) %>% 
  mutate(
    tab = country %>% 
      map(
        \(i){
          
          # i <- cl[[2]]
          
          if(
            i %>% equals("Overall")
          ){
            
            t5 %>%
              filter(
                colfac %>% str_detect("Contrasts")
              ) %>% 
              filter(
                country == "Overall"
              ) %>% 
              extract2("data") %>% 
              extract2(1)
          } else {
            
            t5 %>% 
              filter(
                colfac == "Contrasts" &
                country %>% 
                  equals("Overall") %>% 
                  not &
                country %>% 
                  equals(i) %>% 
                  not
              ) %>% 
              use_series(data) %>% 
              list_rbind
          }
          }
      ),
    mods =  tab %>% 
      map(
        \(i)
        rma(
          yi = i$estimate,
          sei = i$std.error
        )
      ),
    est = mods %>% 
      map_dbl(
        \(i)
        i %>% 
          tidy %>% 
          use_series(estimate) %>% 
          extract2(1)
      ),
    se = mods %>% 
      map_dbl(
        \(i)
        i %>% 
          tidy %>% 
          use_series(std.error) %>% 
          extract2(1)
      ),
    pval = mods %>% 
      map_dbl(
        \(i)
        i %>% 
          tidy %>% 
          use_series(p.value) %>% 
          extract2(1)
      ),
    lab = str_c(
      est %>%
        round(3) %>% 
        as.character %>% 
        str_replace(
          "0.", "."
        ),
      pval %>% 
        gtools::stars.pval() %>% 
        str_replace("\\.", ""),
      " (",
      se %>%
        round(3) %>% 
        as.character %>% 
        str_replace(
          "0.", "."
        ),
      ")"
      )
    )

t7 %<>%  
  left_join(
    t7$mods %>% 
      map(
        \(i)
        i %>% 
          broom::tidy(conf.int = T)
      ) %>% 
      list_rbind %>% 
      mutate(
        country = t7$country
      ) %>% 
      select(country, estimate, starts_with("conf")) %>% 
      rename(
        est_num = estimate
      )
    )

t7$perc_diff <- t7$est %>% 
  subtract(
    t7$est[[1]]
  ) %>% 
  divide_by(
    t7$est[[1]]
    ) %>% 
  multiply_by(100)

t7 %<>%
  arrange(
    perc_diff
  ) %>% 
  filter(
    country %>% str_detect("Overall") %>% not
  ) %>% 
  mutate(
    country = country %>% 
      case_match(
        "ukr" ~ "Ukraine",
        "rus" ~ "Russia",
        "blr" ~ "Belarus",
        "est" ~ "Estonia",
        "kaz" ~ "Kazakhstan"
      ) %>% 
      fct_inorder
  )

t8 <- t7 %>% 
  select(-c(tab, mods)) %>% 
  select(-perc_diff) %>% 
  mutate(colfac = "Correction effect") %>% 
  bind_rows(
    t7 %>% 
      select(country, perc_diff) %>% 
      rename(
        est = perc_diff
      ) %>% 
      mutate(colfac = "Percentage difference from overall effect")
  ) %>% 
  mutate(
    colfac = colfac %>% 
      factor(
        c("Correction effect",
          "Percentage difference from overall effect")
      )
  )

t8 %>% 
  ggplot() +
  geom_linerange(
    aes(
      y = country, 
      xmin = conf.low,
      xmax = conf.high,
    ),
    size = .2,
    data = t8 %>% 
      filter(colfac %>% equals("Correction effect"))
  ) +
  geom_point(
    aes(
      y = country, 
      xmin = conf.low,
      xmax = conf.high,
      x =  est, 
    ),
    size = 1.4,
    data = t8 %>% 
      filter(colfac %>% equals("Correction effect"))
  ) +
  geom_vline(
    aes(xintercept = xint),
    data = tibble(
      xint = -.2124,
      colfac = "Correction effect"
    ) %>% 
      mutate(
        colfac = colfac %>% 
          factor(
            t8$colfac %>% 
              levels
          )
      ),
    size = .5,
    linetype = "dotted"
  ) +
  geom_label(
    aes(x = x,
        y = y, 
        label = lab),
    data = tribble(
      ~x, ~y, ~lab, ~colfac,
      -.21, 1.2, "Overall effect =\n-.212* (.086)", "Correction effect"
      ) %>% 
      mutate(
        colfac = colfac %>% 
          factor(
            t8$colfac %>% 
              levels
          )
      ),
    fill = "grey95",
    label.size = 0,
    position = position_nudge(y = -1),
    size = 2.5,
    family = "Roboto",
    # fontface = "bold.italic"
  ) +
  geom_label(
    aes(
      x = est,
      y = country,
      label = lab
      ),
    family = "Roboto",
    label.size = 0,
    fill = "grey95",
    size = 2.5,
    position = position_nudge(y = .35),
    data = t8 %>% 
      filter(colfac %>% equals("Correction effect"))
    ) +
  geom_vline(
    aes(
      xintercept = xint
    ),
    data = t8 %>% 
      filter(colfac %>% equals("Correction effect") %>% not) %>% 
      mutate(
        xint = 0
      ),
    size = .5,
    linetype = "dashed"
    ) +
  geom_col(
    aes(
      y = country,
      x = est
    ),
    fill = "white",
    color = "grey10",
    size=  .25,
    width = .5,
    data = t8 %>% 
      filter(colfac %>% equals("Correction effect") %>% not)
    ) +
  geom_text(
    aes(
      y = country,
      x = est,
      label = est %>% round
    ),
    size=  3,
    data = t8 %>% 
      filter(
        colfac %>% 
          equals("Correction effect") %>% 
          not &
        est > 0
        ),
    position = position_nudge(x = 3),
    family = "Roboto"
    ) +
  geom_text(
    aes(
      y = country,
      x = est,
      label = est %>% round
    ),
    size=  3,
    data = t8 %>% 
      filter(
        colfac %>% 
          equals("Correction effect") %>% 
          not &
        est < 0
        ),
    position = position_nudge(x = -3),
    family = "Roboto"
  ) +
  geom_segment(
    aes(x, y, xend = xend, yend = yend),
    linewidth = .5,
    arrow = arrow(
      length = unit(0.15,"cm"), 
      type = "closed"),
    data = tribble(
      ~x, ~y, ~xend, ~yend,
      -2, .5,  -25, .5,
      2, .5,  25, .5
    ) %>% 
      mutate(
        colfac = "Percentage difference from overall effect" %>%
          factor(
            t8$colfac %>%
              factor %>%
              levels
          )
        )
  ) +
  geom_richtext(
    aes(x, y, label = label),
    data = tribble(
      ~x, ~y, ~label,
      -20, .5,  "Omitting country makes overall<br>correction magnitude ***smaller***" ,
      20, .5,  "Omitting country makes overall<br>correction magnitude ***bigger***" 
    ) %>% 
      mutate(
        colfac = "Percentage difference from overall effect" %>%
          factor(
            t8$colfac %>%
              factor %>%
              levels
          )
      ),
    color = "transparent",
    text.color = "grey10",
    fill = "grey95", 
    label.size = 0,
    lineheight = .8,
    position = position_nudge(y = -.45),
    size = 2.5,
    family = "Roboto"
  ) +
  scale_y_discrete(
    expand = expansion(add = c(1.3, .75))
  ) +
  facet_wrap(
    ~colfac,
    nrow = 1,
    scales = "free_x"
  ) +
  labs(
    x = "",
    y = "Omitted country"
    ) 