# pid_mom / pid_dad missing
expect_that(load_individuals(pid = c(0L, 1L), is_male = c(0L, 0L)), throws_error('argument "pid_mom" is missing, with no default'))
expect_that(load_individuals(pid = c(0L, 1L), is_male = c(0L, 0L), pid_mom = c(0L, 0L)), throws_error('argument "pid_dad" is missing, with no default'))

# Simple run
db_males <- data.frame(pid = c(1L, 2L, 3L, 4L), pid_dad = c(0L, 1L, 1L, 0L))
pop <- with(db_males, load_individuals(pid = pid, 
                                       is_male = rep(TRUE, nrow(db_males)), 
                                       pid_mom = rep(0L, nrow(db_males)), 
                                       pid_dad = pid_dad, 
                                       progress = FALSE))

# pop_size
expect_identical(pop_size(pop), 4L)

# print.popr_population
expect_identical(capture_output(print(pop)), "Population with 4 individuals")

# `[[`
expect_that(pop[[48L]], throws_error('Individual not found'))
expect_s3_class(pop[[1L]], "popr_individual")
expect_s3_class(pop[[1L]], "externalptr")

