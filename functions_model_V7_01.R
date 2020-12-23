library(deSolve)


##### model ######

seir_model20_07_V1 = function(current_timepoint, state_values, parameters)
{

  # naredimo lokalne spremenljivke
  S = state_values [1]          # susceptibles
  E = state_values [2]          # exposed
  I = state_values [3]          # infectious

  Mild   = state_values [4]        # Recovering mild
  Hosp   = state_values [5]        # recovering hospital
  ICU    = state_values [6]        # recovering ICU
  Dead_s = state_values [7]        # dying symp
  Dead_a = state_values [8]        # dying asymp

  Rmild = state_values [9]        # recovered mild
  Rhosp = state_values [10]       # recovered hospital
  Ricu = state_values [11]        # recovered ICU
  Rdead_s = state_values [12]     # dead
  Rdead_a = state_values [13]     # dead

  #recovered all together Rmild + Rhosp + Ricu

  # parametri
  Bt = parameters[["Bfun"]]  # Bt je funkcija časa

  beta  = Bt[current_timepoint]
  a     = 1/parameters[["D_incubation"]]
  g     = 1/parameters[["D_infectious"]]


  gM     = 1/parameters[["D_recovery_mild"]]

  p_hosp = parameters[["p_hosp"]]
  gHOSP  = 1/parameters[["D_recovery_hosp"]]

  #p_icu  = parameters[["p_ICU"]]
  gICU   = 1/parameters[["D_recovery_ICU"]]



  # prvi in drugi val smrti
  if(current_timepoint <= parameters[["timeV2"]]) {
    p_fatal_a_t = parameters[["p_fatal_a"]]
    p_fatal_s_t = parameters[["p_fatal_s"]]
  }
  else if(current_timepoint <= parameters[["timeV3"]] & current_timepoint > parameters[["timeV2"]]) {
    p_fatal_a_t = parameters[["p_fatal_a_V2"]]
    p_fatal_s_t = parameters[["p_fatal_s_V2"]]
  }
  else {
    p_fatal_a_t = parameters[["p_fatal_a_V3"]]
    p_fatal_s_t = parameters[["p_fatal_s_V3"]]
  }

  # popravek še p_mild
  p_mild = 1 - p_hosp - p_fatal_a_t


  # prvi in drugi val smrti ICU
  if(current_timepoint <= parameters[["timeV2"]]) {
    p_icu_t = parameters[["p_ICU"]]
  }
  else if(current_timepoint <= parameters[["timeV3"]] & current_timepoint > parameters[["timeV2"]]) {
    p_icu_t = parameters[["p_ICU_V2"]]
  }
  else {
    p_icu_t = parameters[["p_ICU_V3"]]
  }


  gF        = 1/parameters[["D_death"]]

  with (
    as.list (parameters),
    {
      # diferencialne enačbe
      dS = (-beta * I ) * S
      dE = (beta * I) * S - (a * E)
      dI = (a * E) - (g * I)

      dMild       =   p_mild * g * I   - gM * Mild
      dHosp       =   p_hosp * g * I   - gHOSP * Hosp
      dICU        =   p_hosp * p_icu_t * g * I   - gICU * ICU
      dDead_s     =   p_fatal_s_t * p_hosp * p_icu_t * g * I   - gICU * Dead_s
      dDead_a     =   p_fatal_a_t * g * I    - gF * Dead_a

      dRmild      =   gM * Mild
      dRhosp      =   gHOSP * Hosp
      dRicu       =   gICU  * ICU
      dRdead_s    =   gICU * Dead_s
      dRdead_a    =   gF * Dead_a



      # rezultati
      results = c (dS, dE, dI, dMild, dHosp, dICU, dDead_s, dDead_a, dRmild, dRhosp, dRicu, dRdead_s, dRdead_a)
      list (results)
    }
  )
}


####### fiksni parametri modela ######

fixed_model_parameters_V7_01 = function() {

  D_incubation      = 5.2
  D_infectious      = 2.9
  D_recovery_mild   = 12
  D_recovery_hosp   = 12
  D_recovery_ICU    = 14
  D_death           = 14


  duration      = 7*12*1e10

  p_hosp        =  0.6 / 100        # delež hospitaliziranih

  p_ICU         = 37 / 100          # delež ICU od hospitaliziranih
  p_ICU_V2      = 20 / 100          # delež ICU od hospitaliziranih v drugem valu
  p_ICU_V3      = 25 / 100          # delež ICU od hospitaliziranih v drugem delu  drugega vala

  p_fatal_s     = 10 / 100          # delež umrlih od ICU
  p_fatal_s_V2  = 3 / 100           # delež umrlih od ICU v 2. valu
  p_fatal_s_V3  = 3 / 100           # delež umrlih od ICU v 2. valu
  p_fatal_a     = 0.22 / 100        # delež smrti asimptomatsko
  p_fatal_a_V2  = 0.01 /100         # delež smrti asimptomatsko drugi val
  p_fatal_a_V3  = 0.01 /100         # delež smrti asimptomatsko drugi val
  timeV2 = 90                       # zacetek drugega vala
  timeV3 = 180                      # drugi del drugega vala

  p_mild        = 1 - p_hosp - p_fatal_a # ni potreben



  parameter_list = list()
  parameter_list$D_incubation = D_incubation
  parameter_list$D_infectious = D_infectious
  parameter_list$D_recovery_mild = D_recovery_mild
  parameter_list$D_recovery_hosp = D_recovery_hosp
  parameter_list$D_recovery_ICU = D_recovery_ICU
  parameter_list$D_death = D_death
  parameter_list$p_hosp = p_hosp
  parameter_list$p_fatal_s = p_fatal_s
  parameter_list$p_fatal_s_V2 = p_fatal_s_V2
  parameter_list$p_fatal_s_V3 = p_fatal_s_V3
  parameter_list$p_fatal_a = p_fatal_a
  parameter_list$p_fatal_a_V2 = p_fatal_a_V2
  parameter_list$p_fatal_a_V3 = p_fatal_a_V3
  parameter_list$timeV2 = timeV2
  parameter_list$timeV3 = timeV3
  parameter_list$p_mild = p_mild
  parameter_list$p_ICU = p_ICU
  parameter_list$p_ICU_V2 = p_ICU_V2
  parameter_list$p_ICU_V3 = p_ICU_V3
  parameter_list$duration=duration

  return(parameter_list)

}


### funkcije Bt #####


Bt_rect_time = function(end_time, Bw, win_len)
{
  if (length(Bw) != length(win_len))
    win_len = rep(10, length(Bw))

  Bfun = rep(0, end_time)
  ct = 1
  for(i in 1:length(Bw)) {

    Bfun[ct:(ct + win_len[i]-1)] = rep(Bw[i], win_len[i])
    ct = ct + win_len[i]
  }

  if (ct < end_time)
    Bfun[ct:end_time] = rep(Bw[length(Bw)], end_time-ct+1)
  else
    Bfun = Bfun[1:end_time]

  return(Bfun)

}



### funkcija za izračun modela s parametri ###


model_seir_V7_01 = function(Bfun, duration_time, param){


  parameter_list = param

  parameter_list$Bfun = Bfun


  # Začetne vrednosti
  N = param$N            # populacija
  X = param$zacetno_stevilo                 # začetno okuženi (20, da ne začnemo na začetku)
  W = N-X                # število dovzetnih
  Z = 0                  # število izpostavljenih



  initial_values = c (S = W/N,
                      E = Z/N,
                      I = X/W,
                      Mild = 0,
                      Hosp = 0,
                      ICU = 0,
                      Dead_s = 0,
                      Dead_a = 0,
                      Rmild = 0,
                      Rhosp = 0,
                      Ricu = 0,
                      Rdead_s = 0,
                      Rdead_a = 0)



  timepoints = seq (1, duration_time, by=1)


  # Simulate SEIR model
  out_seir = lsoda(initial_values, timepoints, seir_model20_07_V1, parameter_list)
  out_seir[,2:14] = round(out_seir[,2:14]*N, digits = 2)


  output = data.frame(out_seir)

  # dnevno stevilo smrti
  output$Rdead_daily = c(0,diff(output$Rdead_s+output$Rdead_a))

  return(output)

}




##### izračun modela z dejanskimi podatki #######


izracun_modela_V7_01 = function(Bfun, sidat, duration_time, param) {

  output = model_seir_V7_01(Bfun, duration_time, param)

  org_st = output$Hosp
  premk_st = output$Hosp
  for (j in 1:length(org_st)) {
    if (j>4)
      premk_st[j-4] = org_st[j]
  }
  output$Hosp = premk_st


  pdat = data.frame(skupine=c(rep("okuženi", nrow(output)),
                              rep("hospitalizirani", nrow(output)),
                              rep("ICU", nrow(output)),
                              rep("umrli: kumulativno", nrow(output)),
                              rep("umrli: dnevno", nrow(output))
                              ),
                    "stevilo"=c(output$I,
                                (output$Hosp+output$ICU),
                                output$ICU,
                                output$Rdead_s + output$Rdead_a,
                                output$Rdead_daily
                                ),
                    dnevi = c(output$time,
                              output$time,
                              output$time,
                              output$time,
                              output$time
                              )
                    )


  prva_smrt = head(output$time[round((output$Rdead_s + output$Rdead_a),0) == 1], n=1)

  datum_smrti = "2020-03-14"
  pdat$datum = as.Date(datum_smrti) + pdat$dnevi - prva_smrt


  # Zložimo podatke, dejanske in simulacijo
  n = sum(pdat$skupine=="hospitalizirani")

  pdat$dejansko = rep(NA, 5*n)

  n_dej = nrow(sidat)
  for (i in seq(1,n_dej,1)) {

    datum = sidat$date[i]

    mask = (pdat$skupine=="okuženi" & pdat$datum == datum)
    pdat$dejansko[mask] = sidat$tests.positive[sidat$date == datum]

    mask = (pdat$skupine=="hospitalizirani" & pdat$datum == datum)
    pdat$dejansko[mask] = sidat$state.in_hospital[sidat$date == datum]

    mask = (pdat$skupine=="umrli: kumulativno" & pdat$datum == datum)
    pdat$dejansko[mask] = sidat$state.deceased.todate[sidat$date == datum]

    mask = (pdat$skupine=="umrli: dnevno" & pdat$datum == datum)
    pdat$dejansko[mask] = sidat$state.deceased.todate[sidat$date == datum]

    mask = (pdat$skupine=="ICU" & pdat$datum == datum)
    pdat$dejansko[mask] = sidat$state.icu[sidat$date == datum]

  }


  # popravimo še dnevno umrle iz podatkov
  mask = (pdat$skupine=="umrli: dnevno")
  pdat$dejansko[mask] = c(0, diff(pdat$dejansko[mask]))

  return(pdat)

}



