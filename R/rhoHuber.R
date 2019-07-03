rhoHuber = function(x,c=2.1){
  # x is a univariate sample
  # c is the tuning constant
  # output is rho(x)
  #
  rho = (x/c)^2
  rho[rho > 1] = 1
  1.54^2*rho
}
