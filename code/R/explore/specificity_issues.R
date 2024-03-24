
#run simulation code first

####### spec for alternative model
rapid_cov <- 0
obj <- new_draws %>%
        mutate(combined=rapid_cov*Vaccinee+(1-rapid_cov)*`Outside Window`) 

ind_spec <- (1-rapid_cov)*mean(alternative_spec_obj$ind_spec) + mean(alternative_vax_obj$ind_spec) * rapid_cov
# comb_spec <- obj$combined %>% mean()
comb_spec <- ind_spec

####### sens for alternative model
sens <- alternative_sens_obj$ind_sens %>% mean()


####### get simulated seropositivity
true_seroprev <- 0.1
seropos <- true_seroprev * sens + (1-true_seroprev) * (1-ind_spec)

##### adjust using sens and comb_spec
final <- (seropos + (comb_spec-1))/(sens+comb_spec-1)

print(final)
