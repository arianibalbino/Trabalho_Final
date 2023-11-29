# clean workspace
rm(list=ls())

library(openxlsx)

# Carregue a biblioteca
# Determinando o número de blocos
#Mínima diferença de importância prática (padronizada): (d∗ = δ∗/σ) = 0.5
# D(diferença máxima que desejamos detectar) / σ = 0.5
#Significância desejada: α = 0.05
# Potência mínima desejada (para o caso d = d∗): π = 1 − β = 0.8-->β = 0.2

#Numero de níveis
a = 2

#phi^2 = (b/2a).(D / σ)^2
#phi^2 = (b/4).(0.5)^2
#phi = 0,25 raiz(b)

#Numerador dos graus de liberdade
#v1 = a -1 = 1
v1 = 1

#V2 = Número de blocos
#V2 = b-1

#Usando 
#b = 64   phi = 2

#64 blocos   usando N suficientemente grande por exemplo = 30.

#Fazer 10 de cada primeiro

# Equipe C
## Config 1
recpars1 <- list(name = "recombination_blxAlphaBeta", alpha = 0, beta = 0)
mutpars1 <- list(name = "mutation_rand", f = 4)

## Config 2
recpars2 <- list(name = "recombination_linear")
mutpars2 <- list(name = "mutation_rand", f = 1.5)

# Obter vetor com as dimensões a serem utilizadas.
# Criar vetor com 64 valores entre 2 e 150
dimensoes <- round(seq(2, 150, length.out = 64))

# Garantir que os valores sejam inteiros positivos
dimensoes <- pmax(dimensoes, 1)

Nblocos = 64   #64
Ntestes = 10   #10

############################################
# Para a configuração 1:
############################################

# Iniciar vetor para Configuração 1
dados_conf <- numeric(0)

# Inicializar uma matriz vazia
matriz_resultados <- matrix(numeric(0), nrow = 0, ncol = Ntestes)

for (j in 1:Nblocos) {
  dim <- dimensoes[j]
  selpars <- list(name = "selection_standard")
  stopcrit <- list(names = "stop_maxeval", maxevals = 5000 * dim, maxiter = 100 * dim)
  probpars <- list(name = "fn", xmin = rep(-5, dim), xmax = rep(10, dim))
  popsize <- 5 * dim
  
  ######################################################################
  # Gera função de Rosenbrock de uma dada dimensão  # Comparar para 2 e 150
  
  suppressPackageStartupMessages(library(smoof))
  fn <- function(X){
    if(!is.matrix(X)) X <- matrix(X, nrow = 1) # <- if a single vector is passed as X
    Y <- apply(X, MARGIN = 1, FUN = smoof::makeRosenbrockFunction(dimensions = dim))
    return(Y)
  }
  ####################################################################
  
  # Inicializar um vetor vazio para execução
  vetor <- numeric(0)
  
  for (i in 1:Ntestes) {
    suppressPackageStartupMessages(library(ExpDE))
    # Run algorithm on problem:
    out <- ExpDE(mutpars = mutpars1,
                 recpars = recpars1,
                 popsize = popsize,
                 selpars = selpars,
                 stopcrit = stopcrit,
                 probpars = probpars,
                 showpars = list(show.iters = "dots", showevery = 20))
    # Store observation:
    
    vetor <- c(vetor, out$Fbest)
  }
  
  matriz_resultados <- rbind(matriz_resultados, vetor)
  media <- mean(vetor)
  dados_conf <- c(dados_conf, media)
}

############################################
# Para a configuração 2:
############################################

# Iniciar vetor para Configuração 2
dados_conf2 <- numeric(0)

# Inicializar uma matriz vazia
matriz_resultados2 <- matrix(numeric(0), nrow = 0, ncol = Ntestes)

for (j in 1:Nblocos) {
  dim <- dimensoes[j]
  selpars <- list(name = "selection_standard")
  stopcrit <- list(names = "stop_maxeval", maxevals = 5000 * dim, maxiter = 100 * dim)
  probpars <- list(name = "fn", xmin = rep(-5, dim), xmax = rep(10, dim))
  popsize <- 5 * dim
  
  ######################################################################
  # Gera função de Rosenbrock de uma dada dimensão  # Comparar para 2 e 150
  
  suppressPackageStartupMessages(library(smoof))
  fn <- function(X){
    if(!is.matrix(X)) X <- matrix(X, nrow = 1) # <- if a single vector is passed as X
    Y <- apply(X, MARGIN = 1, FUN = smoof::makeRosenbrockFunction(dimensions = dim))
    return(Y)
  }
  ####################################################################
  
  # Inicializar um vetor vazio para execução
  vetor2 <- numeric(0)
  
  for (i in 1:Ntestes) {
    suppressPackageStartupMessages(library(ExpDE))
    # Run algorithm on problem:
    out <- ExpDE(mutpars = mutpars2,
                 recpars = recpars2,
                 popsize = popsize,
                 selpars = selpars,
                 stopcrit = stopcrit,
                 probpars = probpars,
                 showpars = list(show.iters = "dots", showevery = 20))
    # Store observation:
    
    vetor2 <- c(vetor2, out$Fbest)
  }
  
  matriz_resultados2 <- rbind(matriz_resultados2, vetor2)
  media2 <- mean(vetor2)
  dados_conf2 <- c(dados_conf2, media2)
}







