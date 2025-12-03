#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Background
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# This is the main paper about using prioritizr for systematic conservation planning: 
# https://conbio.onlinelibrary.wiley.com/doi/full/10.1111/cobi.14376 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set the Working Directory
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# # Ignore this unless you have problems accessing files later (e.g. species.tif)
#
# # If you have any problems it could be related to your path.
# # Check what folder you are working from by using
# getwd()
# # if the folder path you are working from is different from
# # the path where you downloaded `MSc_Protected_Areas then you can change it using the setwd() code below.
# # Make sure you use `/` or `\\` in your path name, otherwise you will get an error.
# setwd("C:/Downloads/MSc_Protected_Areas-main/MSc_Protected_Areas-main/")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load Packages
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(prioritizr)
library(terra)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load the files
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Load species distributions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# We are going to be using Species Distribution Modelling (SDM) data from students last year. 
# We are going to use the binary presence/absence outputs where the probability of a species presence 
# has been translated into a yes/no.
# 
# Each layer describes the spatial distribution of a feature. 
# Here, our feature data correspond to different plant species. 
# Specifically, if the cell (or planning unit) has a 1, a species is predicted to be present in it based on the SDM. 
# If the cell has a 0, then it is predicted to be absent.

# Open file
species = terra::rast("data/species.tif")


### Plot species individual maps
plot(species[[1:3]], nr = 1, axes = FALSE)


### Plot species together
# We have 41 species that we can use to develop a species richness map.
plot(sum(species))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Load Cost Data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# A very important element of protected area planning is the cost. 
# Wherever we put a protected area, there will be a cost involved. 
# 
# Cost can be:
# 
# - no cost: cost is equal across all planning units, which means the prioritisation
# will only rely on the biodiversity data. However, a protected area plan that doesn't 
# account for cost (and therefore people) is not likely to be very helpful.
#   
# - directly estimated: we can evaluate the economic losses of 
# converting land to  protected land. This can be through the cost of purchasing land, 
# the staff needed to patrol, the revenue from tourism, and the loss 
# of opportunities (loss of mining land, grazing land, agricultural land, and even homes).
# 
# - indirectly estimated: we can estimate what land is most costly. 
# For example, the human footprint index is often used to see where people 
# use land. We then make the assumption that areas that people use are 
# most costly and created the greatest conflict, and should be avoided as a place 
# to put a protected area.
# 
# Here, it's an exercise, so we use the Human footprint layer as a proxy for an indirect cost estimate, 
# The layer was downloaded from this paper:
#   https://www.nature.com/articles/s41597-022-01284-8#Sec12


# load the file
cost0 = rast("data/cost.tif")

# Plot the file
par(mfrow = c(1,2)) # put 2 plots side by side
plot(cost0, main = "Cost as Human footprint")
plot(sum(species), main="Species richness")
par(mfrow = c(1,1)) 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Dataset Resolution
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Notice that cost is very fine scale while the species data are coarse. 
# We need to lower the resolution of the cost data so that it is the same as the species data.

# plot the data
par(mfrow=c(1,2)) # split into 2 plots
plot(cost0, main = "Cost High Res")
cost = terra::aggregate(cost0, fact = 12, fun = "median",  na.rm=TRUE)
plot(cost, main = "Cost Low Res")
par(mfrow=c(1,1)) # make sure the next plot is just 1 plot

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Spatial extent sanity checks
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Projection - Do they look the same?
crs(cost, describe=T)
crs(species, describe=T)

## Extent - Do they look the same?
ext(cost)
ext(species)

## Resolution - Do they look the same?
res(cost)
res(species)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Solvers
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# test the installation
library(highs)

# if there's an error or it doesn't load, install HiGHS solver, 
# but this should have been done before the practical
# install.packages("highs", repos = "https://cran.rstudio.com/")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Spatial Prioritisation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Formulate the problem

# I want to protect as many species as possible within my budget by setting cells aside for protection, 
# aiming for 30% of the range of each species
# 
# The things that we need to track are:
# •	Budget: the total cost of the prioritization will represent a 30% of the total land value in the study area. Given this budget, we want the prioritization to increase feature representation, as much as possible,
# •	Features: each feature (i.e. species) would, ideally, have 30% of its distribution covered by the prioritization. In this scenario, we can either purchase all of the land inside a given planning unit, or none of the land inside a given planning unit.
# •	Planning Units: These are sites that we are considering to protect or not. We have the option to protect or not each site. In this case, each site represents one pixel of the raster
# 
# We then need to translate this to code:
# •	We start by creating a new problem() that will use a
# •	minimum shortfall objective to spend the budget (via add_min_shortfall_objective()), while getting as close as possible to a
# •	relative targets of 30% of the range of each species(via add_relative_targets()),
# •	binary decisions to purchase+protect or not (via add_binary_decisions()), and
# •	near-optimal solutions (i.e., 10% from optimality) using the
# •	best solver installed on our computer (via add_default_solver()).



# calculate budget
budget <- terra::global(cost, "sum", na.rm = TRUE)[[1]] * 0.3

# create problem
p1 <-
  problem(cost, features = species) %>%
  add_min_shortfall_objective(budget) %>%
  add_relative_targets(0.3) %>%
  add_binary_decisions() %>%
  add_default_solver(gap = 0.1, verbose = FALSE)

# print problem
print(p1)

## Solve the problem
s1 <- solve(p1)


## Get the time it took to solve in seconds
print(attr(s1, "runtime"))


## Get the status of the solver i.e. was it optimal or not?
print(attr(s1, "status"))


## Plot
plot(s1, main = "Solution", axes = FALSE)


## How many planning units were protected in this solution?
eval_n_summary(p1, s1)

## How much does it cost?
eval_cost_summary(p1, s1)


## Was each species protected 30%?
p1_target_coverage <- eval_target_coverage_summary(p1, s1)
print(p1_target_coverage)


## What percentage of species had at least 30% of their range protected?
print(mean(p1_target_coverage$met) * 100)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Exercise 2
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# # What would be the solution if we require a minimum of 70% of the species range to be protected?

# Hint: change the values in add_relative_targets

# calculate budget
budget <- terra::global(cost, "sum", na.rm = TRUE)[[1]] * 0.3

# create problem
pt <-
  problem(cost, features = species) %>%
  add_min_shortfall_objective(budget) %>%
  add_relative_targets(0.7) %>%
  add_binary_decisions() %>%
  add_default_solver(gap = 0.1, verbose = FALSE)

# print problem
print(pt)

# solve
st <- solve(pt)

#calculate percentage meeting target
pt_target_coverage <- eval_target_coverage_summary(pt, st)
print(pt_target_coverage)
print(mean(pt_target_coverage$met) * 100)

# Note that a lot of species still have 70% of their range protected… 
# How might this be possible when we are only protecting what we can purchase with 30% of our budget 
# (which is close to 30% of the area of Madagascar)? 
# This is because many species will have a range that is smaller than ~30% of Madagascar, 
# so protecting ~30% will protect >30% of many species ranges, especially if they are narrow range species.

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Adding existing PAs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Madagascar already had a Protected Area network. 
# Do the areas we are protecting overlap with the existing protected areas? 
# Perhaps we should develop protected areas that complement and enhance the existing ones.

# load protected areas
PAs = terra::rast("data/locked_in.tif")

# add the land in
PAs = ifel((!is.na(cost0) & is.na(PAs)),0,PAs )

# Re-scale
par(mfrow = c(1,2))
plot(PAs, main="High Res PAs")
PAs = terra::aggregate(PAs, fact = 12, fun = "modal",  na.rm=TRUE) # rescale
plot(PAs, main="Low Res PAs")
par(mfrow = c(1,1))

# Notice that some PAs get swallowed up by the resolution change.

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Modify the problem to add PAs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# The function `add_locked_in_constraints()` 
# allows us to add in areas that are "locked in" 
# i.e. already conserved or protected. 

# create new problem where protected areas are locked in
p2 <-  p1 %>%
  add_locked_in_constraints(PAs) # this is the only new line

# solve the problem
s2 <- solve(p2)

# plot the solution
par(mfrow=c(1,2)) # split plot into 2 panels
plot(s1, main = "without PAs", axes = FALSE)
plot(s2, main = "With PAs", axes = FALSE)
par(mfrow=c(1,1)) # Go back to one plot

# calculate the number of selected planning units in the solution
eval_n_summary(p2, s2)

# calculate the total cost of the solution
eval_cost_summary(p2, s2)

# Note that the scenario 2 is more expensive - why do you think that is?
# 
# It is likely the first scenario has no constraints and protects the cheapest sites first 
# that meet the biodiversity target. However, some of these sites are already protected in scenario 2, 
# so the algorithm has to consider other sites to meet the conservation target. 
# The current PA system may also be sub optimal (for instance, protecting some ecosystems/species more than others), 
# so the algorithm has to compensate by protecting more expensive sites 
# to meet the same biodiversity targets for under-represented ecosystems/species.


# calculate target coverage for the solution
p2_target_coverage <- eval_target_coverage_summary(p2, s2)
print(p2_target_coverage)

# check percentage of the features that have their target met given the solution
print(mean(p2_target_coverage$met) * 100)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Add boundary penalty
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# We can further modify the problem by adding penalties that punish overly fragmented solutions 
# (via add_boundary_penalties()). Here we will use a penalty factor (i.e., boundary length modifier) of 0.003,
# and an edge factor of 50% so that planning units that occur on the outer edge of the study area 
# are not overly penalized.


# create a new problem with boundary penalties added to it
p3 <-
  p2 %>%
  add_boundary_penalties(penalty = 0.003, edge_factor = 0.5)

# solve the problem
s3 <- solve(p3)

# plot the solution
par(mfrow=c(1,3))# split into 3 plots side by side
plot(s1, main = "Baseline", axes = FALSE)
plot(s2, main = "With PAs", axes = FALSE)
plot(s3, main = "With Boundary layer", axes = FALSE)
par(mfrow=c(1,1))# go back to only 1 plot

# calculate target coverage for the solution
p3_target_coverage <- eval_target_coverage_summary(p3, s3)
print(p3_target_coverage)

# check the percentage of the features that have their target met given the solution 
print(mean(p3_target_coverage$met) * 100)


# What is the first thing you notice? 
# It’s much longer… This is because it has to compute a lot more relationships 
# between cells/planning units.
# 
# What is the second thing you notice? 
# The protected areas are huge. 
# It’s important to consider whether very large protected areas make sense in human landscapes, 
# as they could impact livelihoods and displace people. However, large protected areas can also be easier 
# to manage and patrol, and can protect wide-ranging species better. 
# It’s important to think about these things when designing protected area networks.




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Calculate importance scores
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Irreplaceability tells you how important each cell is to the solution: 
# if it is red, you absolutely need to protect that cell. 
# If it is beige, it is interchangeable with other cells. 
# If you have a lot of common species (like we do), then most of them will be beige, 
# but if you have rare narrow range endemics, then those cells will become red 
# (because they only occur in a few locations where they must be protected otherwise they will be lost).


rc <-
  p3 %>%
  eval_ferrier_importance(s3)

# plot the total importance scores
## note that gray cells are not selected by the prioritisation
plot(
  rc[["total"]], main = "Importance scores", axes = FALSE,
  breaks = c(0, 1e-10, 0.005, 0.01, 0.025),
  col = c("#e5e5e5", "#fff7ec", "#fc8d59", "#7f0000")
)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Exercise 3: lock out some areas
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# As an exercise, you will “lock out” some areas from protection. 
# These are places like cities that should be ignored from the prioritisation.
# 
# Go find the locked_out.tif raster in your data directory 
# and use terra::rast to import it. Then use the add_locked_out_constraints() 
# to lock out your raster from the p2 problem we formulated. 
# If you want you can also plot the locked_out data layer



locked_out = terra::rast("data/locked_out.tif")

# make sure no PAs are locked out
locked_out[sum(PAs,locked_out) == 2] = 0 

# create new problem with boundary penalties added to it
p4 <-
  p2 %>%
  add_locked_out_constraints(locked_out) 
  
# solve the problem
s4 <- solve(p4)

# plot the solution
par(mfrow=c(2,2))# split into 4 plots side by side
plot(s1, main = "Baseline", axes = FALSE)
plot(s2, main = "With PAs", axes = FALSE)
plot(s4, main = "With PAs + Avoid cities", axes = FALSE)
plot(s3, main = "With PAs + Boundary layer", axes = FALSE)
par(mfrow=c(1,1))# go back to only 1 plot

# calculate target coverage for the solution
p4_target_coverage <- eval_target_coverage_summary(p4, s4)
print(p4_target_coverage)

# check the percentage of the features that have their target met given the solution
print(mean(p4_target_coverage$met) * 100)

# Notice that in reality it doesn’t change much because we are using the same data 
# to derive the cost layer and the locked out layer. However locking areas out of the prioritisation 
# means that some “expensive” but biodiverse sites cannot be protected, and need to be replaced 
# with many cheaper planning units, meaning that more land has to be protected to compensate.


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Exercise 4: MARXAN
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Marxan addresses conservation problems the opposite way round to prioritizr. 
# It tries to find the cheapest solution (that will impact the least people) while meeting the biodiversity target 
# (here to protect 30% of the range of each species). 
# This can be done by replacing add_min_shortfall_objective(budget) 
# with add_min_set_objective() in the p1 solution. 
# We no longer rely on a budget.

py <-
  problem(cost, features = species) %>%
  # add_min_shortfall_objective(budget) %>% # replace this
  add_min_set_objective() %>% # with this
  add_relative_targets(0.3) %>%
  add_binary_decisions() %>%
  add_locked_in_constraints(PAs) %>%
  add_default_solver(gap = 0.1, verbose = FALSE)

sy <- solve(py)


# plot the solution
par(mfrow=c(1,2))# split into 2 plots side by side
plot(s2, main = "maximise biodiversity", axes = FALSE)
plot(sy, main = "minimise cost", axes = FALSE)
par(mfrow=c(1,1))# go back to only 1 plot


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Exercise 5: Ignore cost
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# So far we have used cost as the human footprint layer, 
# but what if we want to just ignore cost and look only at the biodiversity? 
# We are going to create a new cost layer called cost_null where cost is set to 1 on land, 
# and NA everywhere else. 
# We will then replace cost with cost_null in problem(cost, features = species).
# How does it compare to the layer with cost?
  
# Hint: replace all cost with cost_null in budget and pr


# make the cost 1 on land and NA at sea
cost_null = ifel(!is.na(cost),1,NA)

# calculate budget
budget <- terra::global(cost_null, "sum", na.rm = TRUE)[[1]] * 0.3# replace cost with cost_null here

# create problem
pr <-
  problem(cost_null, features = species) %>% # replacee cost with cost_null in problem(cost, features = species) %>% 
  add_min_shortfall_objective(budget) %>%
  add_relative_targets(0.3) %>%
  add_binary_decisions() %>%
  add_default_solver(gap = 0.1, verbose = FALSE)

sr <- solve(pr)

# changed the plot colours because why not
# plot the solution
par(mfrow=c(1,2))# split into 2 plots side by side
plot(s1, main = "with cost", axes = FALSE)
plot(sr, main = "without cost", axes = FALSE)
par(mfrow=c(1,1))# go back to only 1 plot