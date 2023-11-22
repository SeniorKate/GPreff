case_data_diagnostic <- function(summarised_case_data,
                                 target_dates,
                                 smoothing = c("kernel","wmean","gam")) {

    smoothing <- match.arg(smoothing)
    # get test type proportion matrices
    PCR_prop_matrix <- pivot_test_type_prop_to_wide_matrix(
        summarised_case_data, 'prop_PCR', target_dates)
    RAT_prop_matrix <- pivot_test_type_prop_to_wide_matrix(
        summarised_case_data, 'prop_RAT', target_dates)

    # make validity matrices to mark days with no data (for either test type) as
    # invalid. Discontinuation of RAT reporting is hard coded
    PCR_valid_mat <- PCR_prop_matrix > 0
    RAT_valid_mat <- make_RAT_validity_matrix(RAT_matrix)
    # the above line hard codes end dates of RAT reporting, but also combine
    # that with observed days of non-reporting here
    RAT_valid_mat <- RAT_valid_mat & RAT_prop_matrix > 0

    # make sure invalid RAT or PCR days have NA proportion. This is because when
    # there is no case from either test type, the other test type proportion goes to
    # 1, creating spikes in test type proportion. We need to overcome these spikes
    # by smoothing empirical proportions around them
    RAT_prop_matrix[!RAT_valid_mat | !PCR_valid_mat] <- NA
    # do the same for PCR
    PCR_prop_matrix[!RAT_valid_mat | !PCR_valid_mat] <- NA

    # if any column of prop matrix is entirely NA (because RAT had been entirely
    # switched off for this period), this wouldn't work, we need to manually
    # override it to a constant proportion value (this value doesn't matter since
    # it's constant)
    PCR_prop_matrix <- apply(PCR_prop_matrix,
                             2,
                             FUN = function(x){
                                 if (sum(x, na.rm = TRUE) == 0) {
                                     rep(1,length(x))
                                 } else {
                                     x
                                 } # the sum condition checks if the column is all NA
                             })

    RAT_prop_matrix <- apply(RAT_prop_matrix,
                             2,
                             FUN = function(x){
                                 if (sum(x, na.rm = TRUE) == 0) {
                                     rep(1,length(x))
                                 } else {
                                     x
                                 } # the sum condition checks if the column is all NA
                             })

    # apply smoothing
    if (smoothing == "kernel") {
        # smooth through the proportion. this applies a generic Gaussian kernel smoothing
        # for each column of the proportion matrices
        RAT_prop_matrix <- apply(RAT_prop_matrix,2,
                                 FUN = function(x){
                                     t_idx <- time(x)
                                     ksmooth(t_idx,x,"normal",5)$y
                                 })

        # do the smoothing for PCR
        PCR_prop_matrix <- apply(PCR_prop_matrix,2,
                                 FUN = function(x){
                                     t_idx <- time(x)
                                     ksmooth(t_idx,x,"normal",5)$y
                                 })

    } else if (smoothing == "wmean") {

        # get unique week index
        week_indx <- lubridate::week(target_dates) |>
            as.factor() |>
            as.integer()

        # for each column of the proportion matrices
        RAT_prop_matrix <- apply(RAT_prop_matrix,2,
                                 FUN = function(x){
                                     tapply(x, week_indx,
                                            FUN = function(y) rep(mean(y, na.rm = TRUE), length(y))
                                     ) |> unlist()
                                 })

        # do the smoothing for PCR
        PCR_prop_matrix <- apply(PCR_prop_matrix,2,
                                 FUN = function(x){
                                     tapply(x, week_indx,
                                            FUN = function(y) rep(mean(y, na.rm = TRUE), length(y))
                                     ) |> unlist()
                                 })
        # now apply smoothing
        RAT_prop_matrix <- apply(RAT_prop_matrix,2,
                                 FUN = function(x){
                                     t_idx <- time(x)
                                     ksmooth(t_idx,x,"normal",5)$y
                                 })

        # do the smoothing for PCR
        PCR_prop_matrix <- apply(PCR_prop_matrix,2,
                                 FUN = function(x){
                                     t_idx <- time(x)
                                     ksmooth(t_idx,x,"normal",5)$y
                                 })

    } else if (smoothing == "gam") {
        RAT_prop_matrix <- apply(RAT_prop_matrix,2,
                                 FUN = function(x){
                                     t_idx <- time(x)
                                     mod <- mgcv::gam(x ~ s(t_idx),
                                                      family = "binomial")
                                     predict(mod,newdata = list(t_idx = t_idx), type = "response")
                                 })

        # do the smoothing for PCR
        PCR_prop_matrix <- apply(PCR_prop_matrix,2,
                                 FUN = function(x){
                                     t_idx <- time(x)
                                     mod <- mgcv::gam(x ~ s(t_idx),
                                                      family = "binomial")
                                     predict(mod,newdata = list(t_idx = t_idx), type = "response")
                                 })
    }

    # clamp proportion and hold as constant for end of RAT reporting. This
    # ensures we don't extrapolate something crazy for PCR proportion for this
    # period, thus creating spikes in infection inference.
    # use zoo package for carrying forward last obs to replace NAs
    PCR_prop_matrix <- zoo::na.locf(PCR_prop_matrix)
    RAT_prop_matrix <- zoo::na.locf(RAT_prop_matrix)

    # get rid of any remaing NAs in the middle of proportion
    PCR_prop_matrix <- zoo::na.approx(PCR_prop_matrix)
    RAT_prop_matrix <- zoo::na.approx(RAT_prop_matrix)

    # get rid of any remaing NAs at the start of proportion
    PCR_prop_matrix <- zoo::na.locf(PCR_prop_matrix,fromLast = TRUE)
    RAT_prop_matrix <- zoo::na.locf(RAT_prop_matrix,fromLast = TRUE)

    # NOTE WELL that all these smoothing and na replacement operations
    # produces a lot of INFERRED test type proportion corrections for the
    # case data. All of these are to help with stability of the inferred
    # infection trajectory, but they will NOT allow you to precisely
    # reconstruct the observed case time series using the inferred and
    # smoothed proportions, so don't freak out if case matrix / proportion
    # matrix does not give you the same answer as the observed total number
    # of cases combined between the two test types.

    module(PCR_valid_mat,
           RAT_valid_mat,
           PCR_prop_matrix,
           RAT_prop_matrix)

    }
