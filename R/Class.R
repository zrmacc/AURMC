# Purpose: Class for AURMC output.
# Updated: 2022-08-20

#' AURMC Object
#'
#' @slot Arm0 Results for the reference arm.
#' @slot Arm1 Results for the treatment arm.
#' @slot Contrast Contrasts. 
#' @name AURMC-class
#' @rdname AURMC-class
#' @exportClass AURMC
setClass(
  Class = "AURMC",
  representation = representation(
    Arm0 = "data.frame",
    Arm1 = "data.frame",
    Contrast = "data.frame"
  )
)


#' Print Method for AURMC Object
#'
#' Print method for objects of class \code{AURMC}.
#'
#' @param x An object of class \code{AURMC}.
#' @param ... Unused.
#' @export
print.AURMC <- function(x, ...) {
  
  disp <- function(y) {
    out <- y
    if (is.numeric(y)) {
      dec_part <- (y %% 1)
      if (max(dec_part, na.rm = TRUE) > 0) {
        out <- signif(y, digits = 3)
      }
    }
    return(out)
  }
  
  # Arm 0.
  cat("Arm 0:\n")
  arm0 <- x@Arm0
  arm0[, ] <- lapply(arm0, disp)
  show(arm0)
  cat("\n\n")
  
  # Arm 1.
  cat("Arm 1:\n")
  arm1 <- x@Arm1
  arm1[, ] <- lapply(arm1, disp)
  show(arm1)
  cat("\n\n")
  
  # Contrasts.
  cat("Contrast:\n")
  contrast <- x@Contrast
  contrast[, ] <- lapply(contrast, disp)
  show(contrast)
  cat("\n\n")
}


#' Show Method for AURMC Object
#'
#' @param object An object of class \code{AURMC}.
#' @rdname AURMC-method
#' @importFrom methods show
setMethod(
  f = "show",
  signature = c(object = "AURMC"),
  definition = function(object) {
    print.AURMC(x = object)
  }
)