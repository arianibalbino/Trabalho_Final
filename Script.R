# clean workspace
rm(list=ls())

# Lendo os dados do arquivo CSV
dadosEspeleo <- read.csv('/home/sofia/Documentos/UFMG/Scripts Radio/Dados.csv/Recebida_Espeleo(1).csv')
dadosBase <- read.csv('/home/sofia/Documentos/UFMG/Scripts Radio/Dados.csv/Recebida_Base(2).csv')

############################################################
#Modelos Log-Normal
############################################################

# Função logNormal
logNormal <- function(Pt_dBm, Pr_d0, d0, d, n) {
  PL_d0 <- Pr_d0 - Pt_dBm
  
  # Perda de trajetória que inclui os ganhos das antenas
  PL <- PL_d0 + 10 * n * log10(d / d0)
  
  # Potência recebida em dBm em d metros
  Pr_dBm <- Pt_dBm - PL
  
  return(Pr_dBm)
}

# Função logNormalShadowing
logNormalShadowing <- function(Pt_dBm, Pr_d0, d0, d, sigma, n) {
  PL_d0 <- Pr_d0 - Pt_dBm
  
  # Variável aleatória distribuição normal
  X <- sigma * rnorm(length(d))
  
  # Perda de trajetória 
  PL <- PL_d0 + 10 * n * log10(d / d0) - X
  
  # Potência recebida em dBm em d metros
  Pr_dBm <- Pt_dBm - PL
  
  return(Pr_dBm)
}

n <- 3.2432432432432434
n2 <- 3.823823823823824
sigma <- 4.509778472159174
sigma2 <- 4.889798254864195
Pr_d0 <- dadosEspeleo$Potencia_dBm[1]
d0 <- dadosEspeleo$Distancias[1]
Pr_d02 <- dadosBase$Potencia_dBm[1]
d02 <- dadosBase$Distancias[1]
Pt_dBm = 28
  
# Para d1
d1 <- seq(1, 125, length.out = length(dadosEspeleo$Potencia_dBm))
  
# Para d2
d2 <- seq(1, 125, length.out = length(dadosBase$Potencia_dBm))
    
# Chamadas de função em R
Pr_LogSha <- logNormalShadowing(Pt_dBm, Pr_d0, d0, d1, sigma, n)
Pr_Log <- logNormal(Pt_dBm, Pr_d0, d0, d1, n)
Pr_LogSha2 <- logNormalShadowing(Pt_dBm, Pr_d02, d02, d2, sigma2, n2)
Pr_Log2 <- logNormal(Pt_dBm, Pr_d02, d02, d2, n2)

#############################################################################

# Criar o gráfico de dispersão para dadosEspeleo
plot(dadosEspeleo$Distancias, dadosEspeleo$Potencia_dBm, type = "o", col = "green", pch = 14, main = "Distância vs Potência Recebida", xlab = "Distância", ylab = "Potência (dBm)")

# Adicionar pontos para dadosBase em vermelho
points(dadosBase$Distancias, dadosBase$Potencia_dBm, type = "o", col = "orange", pch = 14)

# Adicionar curva ajustada para dadosEspeleo
lines(d1, Pr_Log, col = "blue", lty = 2)

# Adicionar curva ajustada para dadosBase
lines(d2, Pr_Log2, col = "red", lty = 2)

# Adicionar legenda
legend("topright", legend = c("Espeleo", "Base", "Pr_Log", "Pr_Log2"), col = c("green", "orange", "blue", "red"), pch = c(16, 16, NA, NA), lty = c(NA, NA, 1, 1))
  
##############################################
#Modelo de dois raios
##############################################

dois_raios_modelo <- function(d, R, u, ht, hr, G_los, G_gr, band) {
  d_ref <- sqrt((ht + hr)^2 + d^2)  # Distância ao longo do caminho refletido
  d_los <- sqrt((ht - hr)^2 + d^2)  # Distância ao longo do caminho de linha de visão (LOS)
  
  tworayloss <- 0
  freespaceloss <- 0
  freq <- as.integer(band)
  lam <- 3 * 10^2 / freq  # Comprimento de onda da frequência especificada
  
  # Calcule os valores de perda de caminho para os modelos de Dois Raios e Espaço Livre
  for (i in seq_along(lam)) {
    phi <- 2 * pi * (d_ref - d_los) / lam[i]
    loscoef <- sqrt(G_los) / d_los
    reflcoef <- R * sqrt(G_gr) * exp(-1i * phi) / d_ref
    rs <- lam[i] * (loscoef + reflcoef) / (4 * pi)
    
    tworayloss <- tworayloss + 10 * log10(abs(rs)^2)  # Modelo de Dois Raios
    freespaceloss <- freespaceloss + 20 * log10(lam[i] / (4 * pi * d_los))  # Modelo de Espaço Livre (FSPL)
  }
  
  # Normalize as perdas de caminho em relação ao valor inicial
  freespace_u <- 10 * length(freq) * log10(u^2) + freespaceloss
  norm <- max(tworayloss[1], freespace_u[1])
  
  tworayloss <- tworayloss - norm
  return(tworayloss)
}
R <- -0.9999999999999999
u <- 25.740047715355534
R2 <- -0.9999999999999999
u2 <- 99.740047715355534
ht <- 1.1999999999999997
hr <- 0.5124139298720781
G_los <- 84
G_gr <- 6
band <- '912'

Pr_TwoRay <- dois_raios_modelo(d1, R, u, ht, hr, G_los, G_gr, c(band))
Pr_TwoRay2 <- dois_raios_modelo(d2, R2, u2, ht, hr, G_los, G_gr, c(band))

#############################################################################

# Criar o gráfico de dispersão para dadosEspeleo
plot(dadosEspeleo$Distancias, dadosEspeleo$Potencia_dBm, type = "o", col = "green", pch = 14, main = "Distância vs Potência Recebida", xlab = "Distância", ylab = "Potência (dBm)")

# Adicionar pontos para dadosBase em vermelho
points(dadosBase$Distancias, dadosBase$Potencia_dBm, type = "o", col = "orange", pch = 14)

# Adicionar curva ajustada para dadosEspeleo
lines(d1, Pr_TwoRay, col = "blue", lty = 2)

# Adicionar curva ajustada para dadosBase
lines(d2, Pr_TwoRay2, col = "red", lty = 2)

# Adicionar legenda
legend("topright", legend = c("Espeleo", "Base", "Pr_TwoRay", "Pr_TwoRay"), col = c("green", "orange", "blue", "red"), pch = c(16, 16, NA, NA), lty = c(NA, NA, 1, 1))

##############################################################################
#Análises Estatísticas
#####################################################

# Organize os resultados dos modelos em um dataframe
resultados1 <- data.frame(
  Modelo = c(rep("Pr_LogSha", length(Pr_LogSha)),
             rep("Pr_Log", length(Pr_Log)),
             rep("Pr_TwoRay", length(Pr_TwoRay))),
  Potencia_dBm = c(Pr_LogSha, Pr_Log, Pr_TwoRay)
)

# Execute o teste ANOVA
modelo_anova1 <- aov(Potencia_dBm ~ Modelo, data = resultados1)

# Examine os resultados do teste ANOVA
summary(modelo_anova1)

# Organize os resultados dos modelos em um dataframe
resultados2 <- data.frame(
  Modelo = c(rep("Pr_LogSha2", length(Pr_LogSha2)),
             rep("Pr_Log2", length(Pr_Log2)),
             rep("Pr_TwoRay2", length(Pr_TwoRay2))),
  Potencia_dBm = c(Pr_LogSha2, Pr_Log2,Pr_TwoRay2)
)

# Execute o teste ANOVA
modelo_anova2 <- aov(Potencia_dBm ~ Modelo, data = resultados2)

# Examine os resultados do teste ANOVA
summary(modelo_anova2)



