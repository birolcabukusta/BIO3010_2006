# -------------------------------------------------
# Function to compute volcano statistics per metabolite
# -------------------------------------------------

compute_volcano_stats <- function(long_data,
                                  metabolite_list,
                                  pseudocount_fc = 1) {
   
   # Prepare result table
   result <- data.frame(
      metabolite    = metabolite_list,
      mean_control  = NA_real_,
      mean_ad       = NA_real_,
      log2FC        = NA_real_,
      p_value       = NA_real_,
      stringsAsFactors = FALSE
   )
   
   # Loop over metabolites
   for (i in seq_along(metabolite_list)) {
      
      m <- metabolite_list[i]
      
      values <- long_data$value[long_data$metabolite == m]
      groups <- long_data$group[long_data$metabolite == m]
      
      control_vals <- values[groups == "Control"]
      ad_vals      <- values[groups == "AD"]
      
      result$mean_control[i] <- mean(control_vals, na.rm = TRUE)
      result$mean_ad[i]      <- mean(ad_vals, na.rm = TRUE)
      
      result$log2FC[i] <- log2(
         (result$mean_ad[i] + pseudocount_fc) /
            (result$mean_control[i] + pseudocount_fc)
      )
      
      result$p_value[i] <- suppressWarnings(
         t.test(ad_vals, control_vals)$p.value
      )
   }
   
   return(result)
}
