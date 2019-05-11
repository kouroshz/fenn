# Scaled Linear Discriminant Analysis (slda)
slda is an R package implementing a scaled linear discriminant model. The discriminant directions and the their corresponding optimal weights are identified by optimizing a Von Neumann entropy peanalized distance metric learning probelm. The model can be applied to single class (unsupervied) or multi class data, bridging the gap between PCA, LDA and distance metric learning. For details


## Installation
To install the package directly from github, you need the devtools libray.
```{R}
library(devtools)
install_github("kouroshz/slda", local = FALSE)
```
## Usage
fit <- slda(y~., data)

p <- predict(fit, newdata)