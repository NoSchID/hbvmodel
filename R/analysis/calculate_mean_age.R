# Calculate parameters from demographic model output 02/04/19
# Need to run the model first and define outpath

# Get the mean age within different broader age groups
# Need to decide if wanting to combine into 1 year age groups first

outpath <- out

# Mean age of Gambian adults in 1980 (everyone aged 15 years or more)
sum(subset(outpath$pop, time == 1980)[32:201]*ages[31:200])/
  sum(subset(outpath$pop, time == 1980)[32:201])

# Mean age of Gambian adults in 2015 (everyone aged 15 years or more)
sum(subset(outpath$pop, time == 2015)[32:201]*ages[31:200])/
  sum(subset(outpath$pop, time == 2015)[32:201])


# Other age groups

# Mean age in those 60+
sum(subset(outpath$pop, time == 1980)[122:201]*ages[121:200])/
  sum(subset(outpath$pop, time == 1980)[122:201])
sum(subset(outpath$pop, time == 2015)[122:201]*ages[121:200])/
  sum(subset(outpath$pop, time == 2015)[122:201])

# Mean age in those 30+
sum(subset(outpath$pop, time == 1980)[62:201]*ages[61:200])/
  sum(subset(outpath$pop, time == 1980)[62:201])
sum(subset(outpath$pop, time == 2015)[62:201]*ages[61:200])/
  sum(subset(outpath$pop, time == 2015)[62:201])

# Mean age in those <10
sum(subset(outpath$pop, time == 1980)[2:21]*ages[1:20])/
  sum(subset(outpath$pop, time == 1980)[2:21])
sum(subset(outpath$pop, time == 2015)[2:21]*ages[1:20])/
  sum(subset(outpath$pop, time == 2015)[2:21])

# Mean age in those <11
sum(subset(outpath$pop, time == 1980)[2:23]*ages[1:22])/
  sum(subset(outpath$pop, time == 1980)[2:23])
sum(subset(outpath$pop, time == 2015)[2:23]*ages[1:22])/
  sum(subset(outpath$pop, time == 2015)[2:23])

# Mean age in those <15
sum(subset(outpath$pop, time == 1980)[2:31]*ages[1:30])/
  sum(subset(outpath$pop, time == 1980)[2:31])
sum(subset(outpath$pop, time == 2015)[2:31]*ages[1:30])/
  sum(subset(outpath$pop, time == 2015)[2:31])

# Mean age in those <3
sum(subset(outpath$pop, time == 1980)[2:7]*ages[1:6])/
  sum(subset(outpath$pop, time == 1980)[2:7])
sum(subset(outpath$pop, time == 2015)[2:7]*ages[1:6])/
  sum(subset(outpath$pop, time == 2015)[2:7])

# Mean age in those 17+
sum(subset(outpath$pop, time == 1980)[36:201]*ages[35:200])/
  sum(subset(outpath$pop, time == 1980)[36:201])
sum(subset(outpath$pop, time == 2015)[36:201]*ages[35:200])/
  sum(subset(outpath$pop, time == 2015)[36:201])

# Mean age in those 20+
sum(subset(outpath$pop, time == 1980)[42:201]*ages[41:200])/
  sum(subset(outpath$pop, time == 1980)[42:201])
sum(subset(outpath$pop, time == 2015)[42:201]*ages[41:200])/
  sum(subset(outpath$pop, time == 2015)[42:201])


# Mean age in those 70+
sum(subset(outpath$pop, time == 1980)[142:201]*ages[141:200])/
  sum(subset(outpath$pop, time == 1980)[142:201])
sum(subset(outpath$pop, time == 2015)[142:201]*ages[141:200])/
  sum(subset(outpath$pop, time == 2015)[142:201])

# Mean age in those 16-46
sum(subset(outpath$pop, time == 1980)[34:95]*ages[33:94])/
  sum(subset(outpath$pop, time == 1980)[34:95])
sum(subset(outpath$pop, time == 2015)[34:95]*ages[33:94])/
  sum(subset(outpath$pop, time == 2015)[34:95])

# Mean age in those 16-64
sum(subset(outpath$pop, time == 1980)[34:131]*ages[33:130])/
  sum(subset(outpath$pop, time == 1980)[34:131])
sum(subset(outpath$pop, time == 2015)[34:131]*ages[33:130])/
  sum(subset(outpath$pop, time == 2015)[34:131])

# Mean age in those 20-64
sum(subset(outpath$pop, time == 1980)[42:131]*ages[41:130])/
  sum(subset(outpath$pop, time == 1980)[42:131])
sum(subset(outpath$pop, time == 2015)[42:131]*ages[41:130])/
  sum(subset(outpath$pop, time == 2015)[42:131])

# Mean age in those 18-48
sum(subset(outpath$pop, time == 1980)[38:99]*ages[37:98])/
  sum(subset(outpath$pop, time == 1980)[38:99])
sum(subset(outpath$pop, time == 2015)[38:99]*ages[37:98])/
  sum(subset(outpath$pop, time == 2015)[38:99])

# Mean age in those 19-49
sum(subset(outpath$pop, time == 1980)[40:101]*ages[39:100])/
  sum(subset(outpath$pop, time == 1980)[40:101])
sum(subset(outpath$pop, time == 2015)[40:101]*ages[39:100])/
  sum(subset(outpath$pop, time == 2015)[40:101])


# Mean age in the total population in 1974 (Keneba-Manduar survey)
sum(subset(outpath$pop, time == 1974)*ages)/
  sum(subset(outpath$pop, time == 1974))



