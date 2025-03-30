#' Realiza un Análisis de Componentes Principales (PCA)
#'
#' @param datos base de datos a emplear para el proceso
#' @param escala Booleano: si TRUE, estandariza las variables
#' @param n_comp Numero de componentes principales a extraer (opcional)
#' @return Una lista con tres elementos:
#'   \item{prcomp_obj}{Objeto original de prcomp().}
#'   \item{componentes}{Matriz con las coordenadas de los datos en los nuevos ejes principales.}
#'   \item{var_exp}{Vector con la varianza explicada por cada componente.}
#' @export

pca_realized <- function(datos, escala = TRUE, n_comp = NULL) {
  if (!all(sapply(datos, is.numeric))) {
    stop("Las variables deben ser numéricas para realizar PCA")
  }

  componentesp <- prcomp(datos, center = TRUE, scale. = escala)

  if (is.null(n_comp)) {
    n_comp <- min(ncol(datos), nrow(datos) - 1)
  } else {
    n_comp <- min(n_comp, ncol(componentesp$x))
  }

  var_exp <- (componentesp$sdev^2) / sum(componentesp$sdev^2)

  return(list(
    prcomp_obj = componentesp,
    componentes = componentesp$x[, 1:n_comp, drop = FALSE],
    var_exp = var_exp[1:n_comp]
  ))
}

#' Calcula la contribucion absoluta de las variables a los componentes principales
#'
#' @param pca_obj Objeto resultante de la función prcomp()
#' @return Una matriz con los valores de contribución de cada variable a cada componente principal
#' @export

contrib_pca <- function(pca_obj) {
  cargas_cuadradas <- pca_obj$rotation^2
  varianza_exp <- pca_obj$sdev^2
  contrib <- sweep(cargas_cuadradas, 2, varianza_exp, "*")
  return(contrib)
}

#' Obtiene las cargas factoriales del PCA
#'
#' @param pca_obj Objeto resultante de la función prcomp()
#' @return Una matriz con las cargas factoriales de cada variable en cada componente principal.
#' @export

cargas_pca <- function(pca_obj) {
  return(pca_obj$rotation)
}

#' Calcula la contribución relativa (coseno cuadrado) de las variables a los componentes principales
#'
#' @param pca_obj Objeto resultante de la funcion prcomp()
#' @return Matriz con los cosenos cuadrados de cada variable en cada componente.
#' @export

cos2_pca <- function(pca_obj) {
  cargas_cuadradas <- pca_obj$rotation^2
  contrib_relativa <- sweep(cargas_cuadradas, 2, colSums(cargas_cuadradas), "/")
  return(contrib_relativa)
}

#' Grafica las cargas factoriales del PCA (Biplot)
#'
#' @importFrom ggplot2 aes geom_point geom_text geom_hline geom_vline labs theme_minimal
#' @param pca_obj Objeto resultante de la funcion prcomp()
#' @param comp_x Numero del componente principal a graficar en el eje X (por defecto, PC1).
#' @param comp_y Numero del componente principal a graficar en el eje Y (por defecto, PC2).
#' @return Un grafico que muestra las cargas factoriales de cada variable en los componentes seleccionados.
#' @export

graficar_cargas_pca <- function(pca_obj, comp_x = 1, comp_y = 2) {
  # Extraer las cargas de los componentes seleccionados
  cargas <- as.data.frame(pca_obj$rotation[, c(comp_x, comp_y)])
  cargas$variable <- rownames(cargas)

  # Crear gráfico con ggplot2
  ggplot2::ggplot(cargas, ggplot2::aes(x = .data[[names(cargas)[1]]],
                                       y = .data[[names(cargas)[2]]],
                                       label = variable)) +
    ggplot2::geom_point(color = "blue", size = 3) +
    ggplot2::geom_text(vjust = -0.5, size = 4) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
    ggplot2::labs(
      title = "Biplot de Cargas Factoriales",
      x = paste("Componente Principal", comp_x),
      y = paste("Componente Principal", comp_y)
    ) +
    ggplot2::theme_minimal()
}

#' Análisis completo de Componentes Principales (PCA)
#'
#' @importFrom ggplot2 aes geom_point geom_text geom_hline geom_vline labs theme_minimal
#' @param datos Base de datos con variables numericas a analizar.
#' @param escala Booleano: si TRUE, estandariza las variables antes del PCA.
#' @param n_comp Número de componentes principales a extraer (opcional).
#' @param graficar Booleano: si TRUE, genera gráficos de biplot y scree plot.
#' @return Una lista con detalles completos del análisis de PCA, incluyendo gráficos opcionales.
#' @export

analisis_pca_completo <- function(datos, escala = TRUE, n_comp = NULL, graficar = TRUE) {

  pca_resultado <- pca_realized(datos, escala, n_comp)
  pca_obj <- pca_resultado$prcomp_obj

  cargas <- cargas_pca(pca_obj)
  contrib_abs <- contrib_pca(pca_obj)
  contrib_relativa <- cos2_pca(pca_obj)

  if (graficar) {
    print(graficar_cargas_pca(pca_obj))
  }

  if (graficar) {
    var_exp_df <- data.frame(Componente = seq_along(pca_resultado$var_exp),
                             Varianza = pca_resultado$var_exp)

    scree_plot <- ggplot2::ggplot(var_exp_df, aes(x = Componente, y = Varianza)) +
      ggplot2::geom_bar(stat = "identity", fill = "skyblue") +
      ggplot2::geom_line(aes(group = 1), color = "red", size = 1) +
      ggplot2::geom_point(color = "red", size = 3) +
      ggplot2::labs(title = "Scree Plot - Varianza Explicada",
                    x = "Componente Principal",
                    y = "Proporción de Varianza Explicada") +
      ggplot2::theme_minimal()

    print(scree_plot)
  }

  return(list(
    resumen = summary(pca_obj),
    componentes = pca_resultado$componentes,
    var_exp = pca_resultado$var_exp,
    cargas = cargas,
    contribucion_absoluta = contrib_abs,
    contribucion_relativa = contrib_relativa
  ))
}
