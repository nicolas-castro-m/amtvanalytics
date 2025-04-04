% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PCA.R
\name{pca_realized}
\alias{pca_realized}
\title{Realiza un Análisis de Componentes Principales (PCA)}
\usage{
pca_realized(datos, escala = TRUE, n_comp = NULL)
}
\arguments{
\item{datos}{base de datos a emplear para el proceso}

\item{escala}{Booleano: si TRUE, estandariza las variables}

\item{n_comp}{Numero de componentes principales a extraer (opcional)}
}
\value{
Una lista con tres elementos:
\item{prcomp_obj}{Objeto original de prcomp().}
\item{componentes}{Matriz con las coordenadas de los datos en los nuevos ejes principales.}
\item{var_exp}{Vector con la varianza explicada por cada componente.}
}
\description{
Realiza un Análisis de Componentes Principales (PCA)
}
