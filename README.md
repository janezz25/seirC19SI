Model SEIR C19SI
----------------

Prikazan je potek izračuna modela SEIR C19SI, ki je razvit na podlagi
naslednje strukture:

<img src="SEIR_C19SI.svg" width="500" style="display: block; margin: auto;" />

Model SEIR C19SI je deterministični oddelčni model, ki je bil razvit za
primer spremljanja epidemije COVID19 v Sloveniji. Parametri modela so
kalibirirani na podatke sledilnika do datuma 21.12.2020.

Priprava podatkov in parametrov za zagon modela
-----------------------------------------------

### Nastavitev parametra *β*(*t*)

``` r
Bw = c(1.15,0.5,0.47,0.2,0.1,0.85,0.32,0.33,0.54,0.5, 0.64, 0.34, 0.35, 0.34)
casB = c(15,7,8,35,35,30,30,20,32,10,20,22,12,10)
```

### Podatki iz sledilnika

``` r
file_url = "https://raw.githubusercontent.com/slo-covid-19/data/master/csv/stats.csv"

tdat = read.csv(file_url, header = TRUE, stringsAsFactors = FALSE)
tdat = tdat[5:nrow(tdat),]
tdat$date = as.Date(tdat$date)


sidat = data.frame(date = tdat$date, 
                   tests.positive = tdat$tests.positive,
                   state.in_hospital = tdat$state.in_hospital, 
                   state.deceased.todate = tdat$state.deceased.todate,
                   state.icu = tdat$state.icu)

sidat = sidat[!is.na(sidat$date), ]
```

### Parametri modela

``` r
# default parametri
param = fixed_model_parameters_V7_01()

# popravki
param$N = 2100000
param$zacetno_stevilo = 20
param$D_incubation = 5.2
param$D_infectious = 2.9

# cas 1. in 2. vala
param$timeV2 = 100
param$timeV3 = 200

param$p_fatal_a = 2.07 / 100
param$p_fatal_s = 10 / 100

param$p_fatal_a_V2 = 0.2 / 100
param$p_fatal_s_V2 = 30 / 100

param$p_fatal_a_V3 = 2.3 / 100
param$p_fatal_s_V3 = 30 / 100

param$p_ICU = 37 / 100
param$p_ICU_V2 = 20 / 100
param$p_ICU_V3 = 20 / 100

param$p_hosp = 5.2 / 100
```

Izračun modela
--------------

``` r
# časovni interval izračuna
duration_time = 350

# Iz beta in casa izracunamo potek beta(t)
Bfun = Bt_rect_time(duration_time, Bw, casB)

# glavna funkcija za izračun modela
pdat = izracun_modela_V7_01(Bfun, sidat, duration_time, param)

# izločimo okužene iz rezultata
pdat = pdat[pdat$skupine != "okuženi",]
```

Izris rezultatov
----------------

### Izris rezultatov okoli trenutnega datuma

``` r
# datum izvedbe modela
last_date = as.Date(Sys.time(), "%Y-%m-%d")


# definicija lepšega zapisa števil na y osi
point <- format_format(big.mark = " ", decimal.mark = ",", scientific = FALSE)

pg = ggplot(pdat, aes(x=datum, y = stevilo)) +
  
  geom_line(aes(color = skupine), size=2.2) +
  geom_point(aes(x=datum, y=dejansko, color = skupine), size = 3.0) +
  
  scale_y_continuous(labels = point) +
  
  scale_x_date(breaks = seq(min(pdat$datum),max(pdat$datum), 7),
               labels = format(seq(min(pdat$datum),max(pdat$datum), 7), "%m-%d"),
               limits = c(last_date-45, last_date+45)
  ) +
  
  coord_cartesian(ylim = c(0, 3500)) +
  
  scale_fill_manual(values=c("#ff9900", "#cc2900", "#000000", "#8585ad")) + 
  scale_color_manual(values=c("#ff9900", "#cc2900", "#000000", "#8585ad")) +
  
  ylab("število") +
  
  theme_bw() +
  theme(legend.position = "bottom") + 
  theme(axis.text.x = element_text(angle = 45)) +
  theme(text = element_text(size=20)) 

print(pg)
```

![](model_run_01_files/figure-markdown_github/izris%20rezultatov1-1.png)

### Izris rezultatov od začetka

``` r
pg = ggplot(pdat, aes(x=datum, y = stevilo)) +
  
  geom_line(aes(color = skupine), size=2.2) +
  geom_point(aes(x=datum, y=dejansko, color = skupine), size = 3.0) +
  
  scale_y_continuous(labels = point) +
  
  scale_x_date(breaks = seq(min(pdat$datum),max(pdat$datum), 14),
               labels = format(seq(min(pdat$datum),max(pdat$datum), 14), "%m-%d")
  ) +
  
  coord_cartesian(ylim = c(0, 3500)) +
  
  scale_fill_manual(values=c("#ff9900", "#cc2900", "#000000", "#8585ad")) + 
  scale_color_manual(values=c("#ff9900", "#cc2900", "#000000", "#8585ad")) +
  
  ylab("število") +
  
  theme_bw() +
  theme(legend.position = "bottom") + 
  theme(axis.text.x = element_text(angle = 45)) +
  theme(text = element_text(size=20)) 

print(pg)
```

![](model_run_01_files/figure-markdown_github/izris%20rezultatov2-1.png)

### Izris poteka R

### Izris poteka R okoli trenutnega datuma

``` r
# naredimo datume za R
Rdatum = pdat$datum[1]+seq(0,duration_time-1,1)

# sestavimo Rdata frame za izris
Rdata = data.frame(datum = Rdatum, R = Bfun*param$D_infectious)

p2 = ggplot(Rdata, aes(x=datum)) +
  
  geom_line(aes(y=R), size=2.2) +
  
  geom_hline(yintercept=1, linetype="dashed", color = "red") +
  
  scale_x_date(breaks = seq(min(Rdata$datum),max(Rdata$datum), 7),
               labels = format(seq(min(Rdata$datum),max(Rdata$datum), 7), "%m-%d"),
               limits = c(last_date-45, last_date+45) ) +
               
  scale_y_continuous(limits = c(0,3)) +
  
  ylab("R modela") + xlab("datum") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45)) +
  theme(text = element_text(size=20))



print(p2)
```

![](model_run_01_files/figure-markdown_github/izris%20r1-1.png)

### Izris celotnega poteka R

``` r
p2 = ggplot(Rdata, aes(x=datum)) +
  
  geom_line(aes(y=R), size=2.2) +
  
  geom_hline(yintercept=1, linetype="dashed", color = "red") +
  
  scale_x_date(breaks = seq(min(Rdata$datum),max(Rdata$datum), 14),
               labels = format(seq(min(Rdata$datum),max(Rdata$datum), 14), "%m-%d")
               ) +
               
  scale_y_continuous(limits = c(0,4)) +
  
  ylab("R modela") + xlab("datum") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45)) +
  theme(text = element_text(size=20))



print(p2)
```

![](model_run_01_files/figure-markdown_github/izris%20r2-1.png)
