PredictLabels = function(HIBAG.preds)
{
  
  # create data frame to store results
  name = c("sample.id", "prediction")
  preds = data.frame(matrix(ncol = length(name), nrow = nrow(HIBAG.preds)))
  colnames(preds) = name
  
  # alleles for DQ7 & DQ8. We are only interested in DR3, DR4, DQ7 and DQ8
  DQ7 = c("03:01", "03:04")
  DQ8 = c("03:02", "03:05")
  
  
  for (i in 1:nrow(HIBAG.preds))
  {
    preds["sample.id"][i, ] = HIBAG.preds["sample.id"][i, ]
    
    # in case of allele1 = 01:01 for DRB1, switch positions of allele1 and allele2
    if (substr(HIBAG.preds["HLA.DRB1_allele1"][i, ], 1, 2) == "01")
    {
      temp = HIBAG.preds["HLA.DRB1_allele1"][i, ]
      HIBAG.preds["HLA.DRB1_allele1"][i, ] = HIBAG.preds["HLA.DRB1_allele2"][i, ]
      HIBAG.preds["HLA.DRB1_allele2"][i, ] = temp
    }
    
    DQB1.allele1 = "x"
    DQB1.allele2 = "x"
    DRB1.allele1 = "DRx"
    DRB1.allele2 = "x"
    
    pred.allele1 = "x"
    pred.allele2 = "x"
    
    
    ######################################################################################################################
    # match for allele1
    ######################################################################################################################
    if (substr(HIBAG.preds["HLA.DRB1_allele1"][i, ], 1, 2) == "04") DRB1.allele1 = "DR4"
    else if (substr(HIBAG.preds["HLA.DRB1_allele1"][i, ], 1, 2) == "03") DRB1.allele1 = "DR3"
    
    pred.allele1 = DRB1.allele1
    
    if (HIBAG.preds["HLA.DQB1_allele1"][i, ] %in% DQ7) DQB1.allele1 = "DQ7"
    else if (HIBAG.preds["HLA.DQB1_allele1"][i, ] %in% DQ8) DQB1.allele1 = "DQ8"
    
    
    if (DRB1.allele1 == "DR4" & DQB1.allele1 == "DQ7") pred.allele1 = paste(DRB1.allele1, DQB1.allele1, sep = "-")
    else if (DRB1.allele1 == "DR4" & DQB1.allele1 == "DQ8") pred.allele1 = paste(DRB1.allele1, DQB1.allele1, sep = "-")
    
    ######################################################################################################################
    # match for allele2
    ######################################################################################################################
    if (substr(HIBAG.preds["HLA.DRB1_allele2"][i, ], 1, 2) == "04") DRB1.allele2 = "DR4"
    else if (substr(HIBAG.preds["HLA.DRB1_allele2"][i, ], 1, 2) == "03") DRB1.allele2 = "DR3"
    
    pred.allele2 = DRB1.allele2
    
    if (HIBAG.preds["HLA.DQB1_allele2"][i, ] %in% DQ7) DQB1.allele2 = "DQ7"
    else if (HIBAG.preds["HLA.DQB1_allele2"][i, ] %in% DQ8) DQB1.allele2 = "DQ8"
    
    
    if (DRB1.allele2 == "DR4" & DQB1.allele2 == "DQ7") pred.allele2 = paste(DRB1.allele2, DQB1.allele2, sep = "-")
    else if (DRB1.allele2 == "DR4" & DQB1.allele2 == "DQ8") pred.allele2 = paste(DRB1.allele2, DQB1.allele2, sep = "-")
    
    preds["prediction"][i, ] = paste(pred.allele1, pred.allele2, sep = "/")
    
  }
  
  return(preds)
  
}
