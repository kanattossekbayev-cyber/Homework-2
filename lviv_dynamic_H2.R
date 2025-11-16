## ================================
## 0. Чистим окружение и пакеты
## ================================
rm(list = ls())

library(deSolve)
library(ggplot2)
library(dplyr)
library(tidyr)

## ================================
## 1. Загружаем параметры Львова
## ================================
## ВАЖНО: файл smpl_joint_lviv50K.Rdata должен лежать в этой же папке
load("smpl_joint_lviv50K.Rdata")

## Проверим, что объект есть
ls()   # должно показывать "smpl.joint.l"

## средние значения параметров
par_means <- apply(smpl.joint.l, 2, mean)

## В исходном файле 21 параметр, нам нужно первые 20 (как в таблице)
par_means <- par_means[1:20]

## даём имена параметрам (в порядке, как в статье / твоей таблице)
names(par_means) <- c(
  "S0", "E0", "On0", "Of0", "A0", "Q0",
  "v0", "w0",
  "lam", "lam0",
  "alp.f", "alp0.f", "irr.alp",
  "rho.n", "rho.f",
  "dlt", "gam",
  "inj.On", "inj.Q"
)

par_means

## ================================
## 2. Распаковываем параметры
## ================================
S0   <- par_means["S0"]
E0   <- par_means["E0"]
On0  <- par_means["On0"]
Of0  <- par_means["Of0"]
A0   <- par_means["A0"]
Q0   <- par_means["Q0"]

lam    <- par_means["lam"]
lam0   <- par_means["lam0"]
alp_f  <- par_means["alp.f"]
alp0_f <- par_means["alp0.f"]
irr_alp<- par_means["irr.alp"]
rho_n  <- par_means["rho.n"]
rho_f  <- par_means["rho.f"]
dlt    <- par_means["dlt"]
gam    <- par_means["gam"]

inj_On <- par_means["inj.On"]
inj_Q  <- par_means["inj.Q"]

## Проверка, что числа адекватные
S0; On0; Q0; lam; alp_f; rho_n; gam; inj_On

## ================================
## 3. Доп. параметры и начальные условия
## ================================
mu_OAT_in  <- 0.25   # скорость перехода с листа ожидания на OAT
mu_OAT_out <- 0.05   # выход из OAT в A
aging_rate <- 0.01   # выход из риск-популяции в E

state0 <- c(
  S  = as.numeric(S0),
  On = as.numeric(On0),
  Of = as.numeric(Of0),
  Q  = as.numeric(Q0),
  B  = 142,   # начальное число пациентов на OAT (из статьи/таблицы)
  A  = as.numeric(A0),
  E  = as.numeric(E0)
)

state0

## ================================
## 4. Функция динамической модели (ОДУ)
## ================================
oat_ode <- function(t, state, pars) {
  with(as.list(c(state, pars)), {
    
    N <- S + On + Of + Q + B + A
    
    ## ---- Инициация OUD (peer effect + спонтанная часть) ----
    peer_component <- lam - lam0
    lambda_eff <- lam0 + peer_component * (On + Of + B) / N
    inc_OUD <- lambda_eff * S
    
    ## ---- Переход на лист ожидания OAT ----
    alpha_eff <- alp0_f + (alp_f - alp0_f) * (1 + irr_alp * (On / N))
    to_waitlist <- alpha_eff * On
    
    ## ---- Переход с листа ожидания в OAT (B) ----
    to_OAT <- mu_OAT_in * Q
    
    ## ---- Dropout с листа ожидания ----
    drop_wait <- dlt * Q
    
    ## ---- Рецидивы и прекращение употребления ----
    relapse_from_A  <- rho_n * A
    relapse_from_Of <- rho_f * Of
    quit_from_On    <- gam * On
    quit_from_Of    <- gam * Of
    
    ## ---- Движение через OAT ----
    quit_from_B    <- mu_OAT_out * B
    relapse_from_B <- rho_f * B
    
    ## ---- Старение / выход из риск-популяции ----
    age_out_S  <- aging_rate * S
    age_out_On <- aging_rate * On
    age_out_Of <- aging_rate * Of
    age_out_A  <- aging_rate * A
    
    ## ---- Уравнения ----
    dS  <- -inc_OUD - age_out_S
    dOn <- inc_OUD + relapse_from_A + relapse_from_Of -
      to_waitlist - quit_from_On - age_out_On
    dOf <- relapse_from_B - quit_from_Of - age_out_Of
    dQ  <- to_waitlist - to_OAT - drop_wait
    dB  <- to_OAT - quit_from_B - relapse_from_B
    dA  <- quit_from_On + quit_from_Of + quit_from_B -
      relapse_from_A - age_out_A
    dE  <- age_out_S + age_out_On + age_out_Of + age_out_A
    
    list(c(dS, dOn, dOf, dQ, dB, dA, dE))
  })
}

## ================================
## 5. Сценарии: baseline vs expansion
## ================================
pars_baseline <- c(
  lam = as.numeric(lam),
  lam0 = as.numeric(lam0),
  alp_f = as.numeric(alp_f),
  alp0_f = as.numeric(alp0_f),
  irr_alp = as.numeric(irr_alp),
  rho_n = as.numeric(rho_n),
  rho_f = as.numeric(rho_f),
  dlt = as.numeric(dlt),
  gam = as.numeric(gam),
  mu_OAT_in = mu_OAT_in,
  mu_OAT_out = mu_OAT_out,
  aging_rate = aging_rate
)

pars_expansion <- pars_baseline
pars_expansion["mu_OAT_in"] <- 0.5               # быстрее набираем на OAT
pars_expansion["alp_f"]     <- as.numeric(alp_f) * 1.5

times <- seq(0, 10, by = 0.1)   # 10 лет, шаг 0.1 года

sim_baseline <- as.data.frame(
  ode(y = state0, times = times, func = oat_ode, parms = pars_baseline)
)
sim_baseline$scenario <- "Baseline"

sim_expansion <- as.data.frame(
  ode(y = state0, times = times, func = oat_ode, parms = pars_expansion)
)
sim_expansion$scenario <- "OAT expansion"

sim_all <- bind_rows(sim_baseline, sim_expansion) %>%
  mutate(
    N_risk   = S + On + Of + Q + B + A,
    coverage = B / (On + Of + B + Q + 1e-9),
    inj_total = On * as.numeric(inj_On) + Q * as.numeric(inj_Q)
  )

## ================================
## 6. Графики
## ================================

## 6.1 Population processes
pop_long <- sim_all %>%
  select(time, scenario, S, On, Q, B, A) %>%
  pivot_longer(cols = c(S, On, Q, B, A),
               names_to = "compartment",
               values_to = "count")

## 6.1 Population processes
p1 <- ggplot(pop_long, aes(x = time, y = count, colour = compartment)) +
  geom_line() +
  facet_wrap(~ scenario, ncol = 1, scales = "free_y") +
  labs(x = "Years since 2016",
       y = "Number of individuals",
       title = "Population processes in Lviv (baseline vs OAT expansion)") +
  theme_bw()
print(p1)

## 6.2 Coverage
p2 <- ggplot(sim_all, aes(x = time, y = coverage, colour = scenario)) +
  geom_line(linewidth = 1) +
  labs(x = "Years since 2016",
       y = "Treatment coverage",
       title = "Long-term effects of OAT expansion on coverage in Lviv") +
  theme_bw()
print(p2)

## 6.3 Injections (proxy for HIV/HCV risk)
p3 <- ggplot(sim_all, aes(x = time, y = inj_total, colour = scenario)) +
  geom_line(linewidth = 1) +
  labs(x = "Years since 2016",
       y = "Total injections per year (proxy)",
       title = "Predictive simulations: injections under baseline vs expansion") +
  theme_bw()
print(p3)
