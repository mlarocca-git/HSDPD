fit.sdpd.model <- function(series, model, check=FALSE, two.stage=FALSE, NAcovs="pairwise.complete.obs", parallelize=FALSE) 
{
  
  if(parallelize){
    ## costruzione lista di oggetti sdpd-serie su gruppi parallelizzati
    df.data <- parallelize.sdpd.data(series, model=model)
  }
  else df.data <- list(series)

  ## estimate the H-SDPD model parameters on each group of pixels
  res    <- df.data |>
    map(fit.sdpd.procedure, model=model, check=check, two.stage=two.stage, NAcovs=NAcovs) |>
    map(fit.sdpd.assemble) |>
    bind_rows()
  
  ## output
  res
}