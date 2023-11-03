library(tidyverse) 
library(estimatr)
library(magrittr)
library(estimatr)
library(emmeans)
library(broom)
library(metafor)

t1 <- "https://github.com/thomasjwood/corrections_eeurope/raw/main/data/t1.rds" %>% 
  url %>% 
  readRDS

# generate corrections seen indicator

tcs <- t1 %>% 
  select(
    response_id, 
    matches("fc_\\d_treatment"),
    matches("fc_\\d_control")
  ) %>% 
  gather(
    cond, ans, -1
  ) %>% 
  arrange(response_id) %>% 
  filter(ans != "0") 

# keep respondents who have 3 answers

tcs %<>% 
  filter(
    response_id %>% 
      is_in(
        tcs %>% 
          group_by(response_id) %>% 
          tally %>% 
          filter(
            n == 3
          ) %>% 
          use_series(response_id)
      )
  )

# nesting table for modelling

t2 <- tcs %>% 
  group_by(
    response_id
  ) %>% 
  arrange(response_id, cond) %>% 
  summarize(
    cs = cond %>% 
      str_detect("_treatment") %>% 
      sum %>% 
      factor(
        levels = as.character(0:3)
      )
  ) %>% 
  right_join(
    t1 %>% 
      mutate(
        country = case_when(
          blr_survey %>% equals(1) ~ "Belarus",
          est_survey %>% equals(1) ~ "Estonia",
          kaz_survey %>% equals(1) ~ "Kazakhstan",
          rus_survey %>% equals(1) ~ "Russia",
          ukr_survey %>% equals(1) ~ "Urkaine")
      ) %>% 
      select(
        response_id, rus_support, ukr_support,
        country
      )
  ) %>% 
  pivot_longer(
    cols = ends_with("_support"),
    names_to = "support_itm",
    values_to = "ans"
  )

#  include pooled nations

t2 %<>% 
  mutate(country =  "Pooled") %>% 
  bind_rows(t2) %>% 
  group_by(
    country, support_itm
  ) %>% 
  nest

# models

t2$mods <- t2$data %>% 
  map(
    \(i)
    lm_robust(
      ans ~ cs, 
      data = i
    )
  )

t2$emm <- t2$mods %>% 
  map(
    \(i)
    
    emmeans(
      i, 
      pairwise ~ cs
    )
  )


t2$fitted <- t2$emm %>% 
  map(
    \(i)
    
    # i <- t2$emm[[1]]
    
    i$emmeans %>% 
      tidy(conf.int = T)
  )

t3 <- t2 %>% 
  ungroup %>% 
  select(country, support_itm, fitted) %>% 
  pmap(
    \(country, support_itm, fitted)
    
    fitted %>% 
      mutate(
        country, support_itm
      )
  ) %>% 
  list_rbind %>% 
  mutate(
    support_itm = support_itm %>% 
      str_detect("rus") %>% 
      ifelse("...the Russian side...",
             "...the Ukrainian side...")
  )

# ordering countries

t3$country %<>% 
  factor(
    t3 %>% 
      group_by(
        support_itm, country
      ) %>% 
      summarize(mu = estimate %>% mean) %>% 
      mutate(support_itm = support_itm %>% str_extract("[A-Z]")) %>% 
      spread(support_itm, mu) %>% 
      mutate(df = R - U) %>% 
      arrange(df) %>% 
      extract2("country") %>% 
      extract(-3) %>% 
      c("Pooled")
    )

t3 %>% 
  ggplot() +
  geom_linerange(
    aes(
      cs, ymin = conf.low, ymax = conf.high
      ),
    linewidth = .375
    ) +
  geom_point(
    aes(
      cs, estimate
    ),
    shape = 21,
    size = 6,
    fill = "white"
  ) +
  geom_text(
    aes(
      cs, estimate, label = estimate %>% round(2)
    ),
    family = "Roboto",
    size = 1.6,
  ) +
  labs(
    # 
    x = "Corrections seen",
    y = "Agreement (5pt scale)",
    title = "Do you agree with the following statement? I support ... in the war in Ukraine"
  ) +
  facet_grid(
    support_itm ~ country, 
    labeller = label_wrap_gen(width = 10) 
      ) 

library(gt)
library(modelsummary)

t2 %>% 
  filter(
    support_itm %>% 
      str_detect("rus_")
  ) %>% 
  use_series(mods) %>%
  set_attr(
    "names",
    c("Pooled",
      "Ukraine",
      "Russia",
      "Belarus",
      "Kazakhstan",
      "Estonia")
    ) %>% 
  modelsummary::modelsummary(
    stars =  c(
      '*' = .05, 
      '**' = .01,
      '***' = .001
    ), 
    coef_rename = c(
      "Intercept",
      "Corrections seen = 1",
      "Corrections seen = 2",
      "Corrections seen = 3"
      ),  
    fmt = function(i) i %>% round(2),  
    gof_map = c("nobs", "r.squared")
    )

t2 %>% 
  filter(
    support_itm %>% 
      str_detect("ukr_")
  ) %>% 
  use_series(mods) %>%
  set_attr(
    "names",
    c("Pooled",
      "Ukraine",
      "Russia",
      "Belarus",
      "Kazakhstan",
      "Estonia")
  ) %>% 
  modelsummary::modelsummary(
    stars =  c(
      '*' = .05, 
      '**' = .01,
      '***' = .001
    ), 
    coef_rename = c(
      "Intercept",
      "Corrections seen = 1",
      "Corrections seen = 2",
      "Corrections seen = 3"
    ),  
    fmt = function(i) i %>% round(2),  
    gof_map = c("nobs", "r.squared"),
  )


