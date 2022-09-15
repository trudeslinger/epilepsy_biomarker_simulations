#!/usr/bin/env Rscript
#
# Simulate 100 studies with different group sizes and effect sizes
# Save histograms and ROC-curves for the first iteration (group size == 400)
#
# W.M. Otte (w.m.otte@umcutrecht.nl) 
# G. Slinger (g.slinger-2@umcutrecht.nl)

################################################################################

# load required packages
library( "effectsize" )
library( "ggpubr" )
library( "cowplot" )

################################################################################
# BEGIN FUNCTIONS
################################################################################

number_ticks <- function( n ) { function( limits ) pretty( limits, n ) }

###
# simulate biomarker
###
sim_biomarker <- function( N, mean_control, mean_intervention, sd )
{
  # sample sizes
  n_C <- N
  n_I <- N
  
  # get data
  y_C <- abs( rnorm( n_C, mean_control, sd ) )
  y_I <- abs( rnorm( n_I, mean_intervention, sd ) )
  
  # groups
  group_C <- rep( "Control", n_C )
  group_I <- rep( "Intervention", n_I )
  
  df <- data.frame( group = as.factor( c( group_C, group_I ) ), y = c( y_C, y_I ) )
  
  return( df )
}

###
# get 2x2 table
###
get_two_by_two_table <- function( predictor, response, threshold = 0.5 )
{
  # continuous predictor -> binary predictor, based on threshold
  pred_bin <- rep.int( 0, times = length( predictor ) )
  pred_bin[ predictor > threshold ] <- 1
  
  TP <- sum( pred_bin == 1 & response == "Intervention" )
  TN <- sum( pred_bin == 0 & response == "Control" )
  
  FP <- sum( pred_bin == 1 & response == "Control" )
  FN <- sum( pred_bin == 0 & response == "Intervention" )
  
  # 2x2 Table
  mtx <- as.table( matrix( c( TP, FP, FN, TN ), nrow = 2, byrow = TRUE,
                           dimnames=list( predicted = c( "Intervention", "Control" ),
                                          reference = c( "Intervention", "Control" ) ) ) ) 
  
  return( mtx )
}

###
# get performance with 95% confidence intervals
###
get_performance <- function( predictor, y, threshold = 0.5 )
{
  # get roc curve, levels = 'Control', 'Epilepsy', direction: controls < cases
  roc.rf <- pROC::roc( response = y, predictor = predictor, levels = base::levels( as.factor( y ) ), direction = "<" )
  
  # get mean + 95% CIs for AUC
  aucRaw <- pROC::ci.auc( roc.rf )
  
  # add AUC to data.frame
  aucMean <- aucRaw[ 2 ]
  aucLower <- aucRaw[ 1 ]
  aucUpper <- aucRaw[ 3 ]
  
  # binarize
  tdata <- rep( 0, length( y ) )
  tdata[ predictor >= threshold ] <- 1
  tdata <- as.factor( tdata )
  treference <- as.factor( as.integer( as.factor( y ) ) - 1 )
  
  # if only single outcome in data, no 2x2 table could be constructed
  if( length( unique( levels( tdata ) ) ) < 2 )
  {
    return( NA )
  }
  if( length( unique( levels( treference ) ) ) < 2 )
  {
    return( NA )
  }
  
  # get 2x2 table
  mtx <- get_two_by_two_table( predictor = predictor, response = y )
  
  TP <- mtx[ "Intervention", "Intervention" ]
  TN <- mtx[ "Control", "Control" ]
  
  FP <- mtx[ "Intervention", "Control" ]
  FN <- mtx[ "Control", "Intervention" ]
  
  # total
  N = TP + TN + FP + FN
  
  ## sensitivity with 95% confidence interval 
  sensMean <- caret::sensitivity( mtx )
  sens_errors <- sqrt( caret::sensitivity( mtx ) * ( 1 - caret::sensitivity( mtx ) ) / sum( mtx[ , 1 ] ) )
  sensLower <- caret::sensitivity( mtx ) - 1.96 * sens_errors
  sensUpper <- caret::sensitivity( mtx ) + 1.96 * sens_errors
  
  ## specificity with 95% confidence interval 
  specMean <- caret::specificity( mtx )
  spec_errors <- sqrt( caret::specificity( mtx ) * ( 1 - caret::specificity( mtx ) ) / sum( mtx[ , 2 ] ) )
  specLower <- caret::specificity( mtx ) - 1.96 * spec_errors
  specUpper <- caret::specificity( mtx ) + 1.96 * spec_errors
  
  ## positive predictive value with 95% confidence interval 
  ppvMean <- caret::posPredValue( mtx )
  ppv_errors <- sqrt( caret::posPredValue( mtx ) * ( 1 - caret::posPredValue( mtx ) ) / sum( mtx[ 1 , ] ) )
  ppvLower <- caret::posPredValue( mtx ) - 1.96 * ppv_errors
  ppvUpper <- caret::posPredValue( mtx ) + 1.96 * ppv_errors
  
  ## negative predictive value with 95% confidence interval 
  npvMean <- caret::negPredValue( mtx )
  npv_errors <- sqrt( caret::negPredValue( mtx ) * ( 1 - caret::negPredValue( mtx ) ) / sum( mtx[ 2 , ] ) )
  npvLower <- caret::negPredValue( mtx ) - 1.96 * npv_errors
  npvUpper <- caret::negPredValue( mtx ) + 1.96 * npv_errors
  
  ## accuracy with 95% confidence interval
  additional <- caret::confusionMatrix( mtx )
  acc <- additional$overall[ c( "Accuracy", "AccuracyLower", "AccuracyUpper" ) ]
  accMean <- acc[ 1 ]
  accLower <- acc[ 2 ]
  accUpper <- acc[ 3 ]
  
  if( sensLower < 0 )
    sensLower <- 0
  
  if( specLower < 0 )
    specLower <- 0
  
  if( ppvLower < 0 )
    ppvLower <- 0
  
  if( npvLower < 0 )
    npvLower <- 0
  
  # collect output
  out <- cbind( threshold, N, TP, TN, FP, FN,
                aucMean = round( aucMean, 2 ), aucLower = round( aucLower, 2 ), aucUpper = round( aucUpper, 2 ), 
                accMean = round( accMean, 2 ), accLower = round( accLower, 2 ), accUpper = round( accUpper, 2 ),  
                sensMean = round( sensMean, 2 ), sensLower = round( sensLower, 2 ), sensUpper = round( sensUpper, 2 ),  
                specMean = round( specMean, 2 ), specLower = round( specLower, 2 ), specUpper = round( specUpper, 2 ), 
                ppvMean = round( ppvMean, 2 ), ppvLower = round( ppvLower, 2 ), ppvUpper = round( ppvUpper, 2 ), 
                npvMean = round( npvMean, 2 ), npvLower = round( npvLower, 2 ), npvUpper = round( npvUpper, 2 ) ) 
  
  rownames( out ) <- NULL
  out <- data.frame( out )
  
  return( out )
}

################################################################################
# END FUNCTIONS
################################################################################

# outdir 
outdir <- "out.sim"
dir.create( outdir, showWarnings = FALSE )

# seed
set.seed( 123 )

# number of iterations 
niter <- 100

# baseline mean + sd
mean_control <- 100
sd <- 15

# Cohen d's to test [small, medium, large, very large]
cohens_ds <- c( 0.25, 0.50, 0.75, 1.0, 1.25, 1.50, 1.75, 2.00, 2.25, 2.5 )

# group sample sizes [single group]
N_trains <- c( 25, 50, 100, 200, 400, 800 )

# containers
all <- NULL
roc_auc <- NULL
roc_ci <- NULL

# plot vector containers and counters
plots <- vector( mode = 'list', length = length( cohens_ds ) )
plots_roc <- vector( mode = 'list', length = length( cohens_ds ) )
counter <- 1
counter_roc <- 1

# get performance + create plots
# loop over sample size
for( N_train in N_trains )
{
  # loop over effect size
  for( cohens_d in cohens_ds )
  {
    # loop over iterations
    for( iter in 1:niter )
    {
      # mean intervention group
      mean_intervention <- mean_control + ( sd * cohens_d ) 
      
      # simulate data
      df_train <- sim_biomarker( N = N_train, mean_control, mean_intervention, sd )
      df_test <-  sim_biomarker( N = 500, mean_control, mean_intervention, sd )
 
      # get histogram example for first iteration
      if( iter == 1 )
      {
        hbins <- 30
        if( N_train < 100 )
          hbins <- floor( N_train / 4 )
        
        # plot
        p <- ggpubr::gghistogram( df_train, x = "y", add = "mean", bins = hbins,
                                  rug = TRUE, color = "group", fill = "group", 
                                  ylab = "Subjects", xlab = "Outcome",
                                  palette = c( "#ce8080", "#5698a3" ),
                                  binwidth = 3.5 ) +
          scale_x_continuous( breaks = number_ticks( 8 ) ) +
          scale_y_continuous( breaks = number_ticks( 8 ) )
        
        # save subset for multi plot
        if( N_train == 400 )
        {
          plots[[ counter ]] <- list( N_train = N_train, cohens_d = cohens_d, p = p )
          counter <- counter + 1
        }
      }
      
      # fit prediction model
      model <- glm( group ~ y, data = df_train, family = "binomial" )
      sum <- summary( model )
      z_value <- sum$coefficients[ "y", "z value" ]
      p_value <- sum$coefficients[ "y", "Pr(>|z|)" ]
      
      # predict
      df_test$pred <- predict( model, newdat = df_test, type = "response" )
      
      # performance
      perf <- get_performance( pred = df_test$pred, y = df_test$group )
      perf$z_value <- z_value
      perf$p_value <- p_value
      perf$iter <- iter
      perf$cohens_d <- cohens_d
      perf$N_train <- N_train

      # save ROC-curve for first iteration
      if( iter == 1 & N_train == 400 )
      {
        # get roc curve, levels = 'Control', 'Intervention', direction: controls < cases
        roc <- pROC::roc( response = df_test$group, predictor = df_test$pred, 
                          levels = base::levels( as.factor( df_test$group ) ), direction = "<" )
        
        # prepare plotting dataframe (with sens/spec as percentage)
        data <- data.frame( spec = 100 * roc$specificities, sens = 100 * roc$sensitivities )
        
        # roc curve with inverted specificity axis
        p_roc <- ggplot( data = data, aes( x = spec, y = sens ) ) +
          geom_line( size = 0.8, colour = "#5698a3" ) +
          geom_abline( intercept = 100, slope = 1, size = 1, linetype = 3, colour = "grey50" ) +
          xlab( "Specificity" ) +				
          ylab( "Sensitivity" ) +
          scale_y_continuous( breaks = number_ticks( 8 ) ) +
          scale_x_reverse( breaks = number_ticks( 8 ) ) +
          theme_classic() +
          theme( aspect.ratio = 1 )
        
        # add to container
        plots_roc[[ counter_roc ]] <- list( N_train = N_train, cohens_d = cohens_d, p = p_roc )
        roc_auc <- rbind( roc_auc, roc$auc[ 1 ] )
        roc_ci <- rbind( roc_ci, pROC::ci.auc( roc, conf.level = 0.95 ) )
        counter_roc <- counter_roc + 1
      }
      
      # add to container
      all <- rbind( all, perf )
    }
  }
}

# combine histograms + roc curves 
# N fixed at 400 subjects per group
p_comb_hist_roc_n400 <- plot_grid( 
  
  # plots at Cohen's d 0.25
  plots[[ 1 ]]$p + 
    theme( legend.position = "none" ) +
    theme( axis.title = element_text( size = 12, colour = "black" ), axis.text = element_text( size = 11, colour = "black" ) ) +
    coord_cartesian( xlim = c( 30, 190 ), ylim = c( 0, 70 ) ) +
    annotate( "text", label = paste0( "Cohen's d: ", plots[[ 1 ]]$cohens_d ), x = 56, y = 56, size = 4, colour = "gray30" ),
  
  plots_roc[[ 1 ]]$p + 
    theme( legend.position = "none" ) +
    theme( axis.title = element_text( size = 12, colour = "black" ), axis.text = element_text( size = 11, colour = "black" ) ) +
    annotate( "text", label = paste0( "AUC: ", format( round( roc_ci[ 1, 2 ], digits = 2 ), nsmall = 2 ),
                                      " (CI: ", format( round( roc_ci[ 1, 1 ], digits = 2), nsmall = 2 ),
                                      "-", format( round( roc_ci[ 1, 3 ], digits = 2), nsmall = 2 ),
                                      ")"), x = 30, y = 5, size = 4, colour = "gray30" ),
  
  # plots at Cohen's d 0.75
  plots[[ 3 ]]$p + 
    theme( legend.position = "none" ) +
    theme( axis.title = element_text( size = 12, colour = "black" ), axis.text = element_text( size = 11, colour = "black" ) ) +
    coord_cartesian( xlim = c( 30, 190 ), ylim = c( 0, 70 ) ) +
    annotate( "text", label = paste0( "Cohen's d: ", plots[[ 3 ]]$cohens_d ), x = 56, y = 56, size = 4, colour = "gray30" ),
  
  plots_roc[[ 3 ]]$p + 
    theme( legend.position = "none" ) +
    theme( axis.title = element_text( size = 12, colour = "black" ), axis.text = element_text( size = 11, colour = "black" ) ) +
    annotate( "text", label = paste0( "AUC: ", format( round( roc_ci[ 3, 2 ], digits = 2 ), nsmall = 2 ),
                                      " (CI: ", format( round( roc_ci[ 3, 1 ], digits = 2), nsmall = 2 ),
                                      "-", format( round( roc_ci[ 3, 3 ], digits = 2), nsmall = 2 ),
                                      ")"), x = 30, y = 5, size = 4, colour = "gray30" ),
  
  # plots at Cohen's d 1.25
  plots[[ 5 ]]$p + 
    theme( legend.position = "none" ) +
    theme( axis.title = element_text( size = 12, colour = "black" ), axis.text = element_text( size = 11, colour = "black" ) ) +
    coord_cartesian( xlim = c( 30, 190 ), ylim = c( 0, 70 ) ) +
    annotate( "text", label = paste0( "Cohen's d: ", plots[[ 5 ]]$cohens_d ), x = 56, y = 56, size = 4, colour = "gray30" ),
  
  plots_roc[[ 5 ]]$p + 
    theme( legend.position = "none" ) +
    theme( axis.title = element_text( size = 12, colour = "black" ), axis.text = element_text( size = 11, colour = "black" ) ) +
    annotate( "text", label = paste0( "AUC: ", format( round( roc_ci[ 5, 2 ], digits = 2 ), nsmall = 2 ),
                                      " (CI: ", format( round( roc_ci[ 5, 1 ], digits = 2), nsmall = 2 ),
                                      "-", format( round( roc_ci[ 5, 3 ], digits = 2), nsmall = 2 ),
                                      ")"), x = 30, y = 5, size = 4, colour = "gray30" ),
  
  # plots at Cohen's d 2.50
  plots[[ 10 ]]$p +
    theme( legend.position = "none" ) +
    theme( axis.title = element_text( size = 12, colour = "black" ), axis.text = element_text( size = 11, colour = "black" ) ) +
    coord_cartesian( xlim = c( 30, 190 ), ylim = c( 0, 70 ) ) +
    annotate( "text", label = paste0( "Cohen's d: ", plots[[ 10 ]]$cohens_d ), x =56, y = 56, size = 4, colour = "gray30" ),
  
  plots_roc[[ 10 ]]$p + 
    theme( legend.position = "none" ) +
    theme( axis.title = element_text( size = 12, colour = "black" ), axis.text = element_text( size = 11, colour = "black" ) ) +
    annotate( "text", label = paste0( "AUC: ", format( round( roc_ci[ 10, 2 ], digits = 2 ), nsmall = 2 ),
                                      " (CI: ", format( round( roc_ci[ 10, 1 ], digits = 2), nsmall = 2 ),
                                      "-", format( round( roc_ci[ 10, 3 ], digits = 2), nsmall = 2 ),
                                      ")"), x = 30, y = 5, size = 4, colour = "gray30" ),
  
  # plot layout
  labels = c( "A", "", "B", "", "C", "", "D", "" ), 
  label_size = 12,
  rel_widths = 1, 
  rel_heights = 1,
  nrow = 4, ncol = 2 ) +
  
  theme( plot.background = element_rect( fill = "white", colour = NA ) )

# save plot
save_plot( p_comb_hist_roc_n400, file = paste0( outdir, "/FACET_hist_roc_n400.png" ), ncol = 2, nrow = 4, base_width = 4.5, base_height = 3.7 )

# write simulated data to file
write.csv( all, file = paste0( outdir, "/simulated_perf.csv" ) )
