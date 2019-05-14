---
bibliography: references.bib
---

# Free Energy Nearest Neighbour (FENN)
FENN is an R package implementing a Von Neumann Entropy based distance metric learning model. The method identifies informative and non-informative features in the data. It scales and projects the data onto a subspace of maximal separation between classes. The corresponding optimal weights on each dimension are identified by optimizing the Von Neumann entropy. The model can be applied to single class (unsupervised) or multi class data, bridging the gap between PCA, LDA and distance metric learning. 


## Installation
To install the package directly from GitHub, you need the devtools library.
```{R}
library(devtools)
install_github("kouroshz/fenn", local = FALSE)
```
## Usage
```
fit <- fenn(y~., data)
p <- predict(fit, newdata)
```

## Example
```{R}
library(fenn)
Iris <- data.frame(rbind(iris3[,,1], iris3[,,2], iris3[,,3]),Sp = rep(c("s","c","v"), rep(50,3)))
train <- sample(1:150, 75)
table(Iris$Sp[train])
z <- fenn(Sp ~ ., Iris, prior = c(1,1,1)/3, subset = train)
predict(z, Iris[-train, ])$class
print(z)
plot(z)
```
# Citation
---
references:
- id: fenner2012a
  title: One-click science marketing
  author:
  - family: Fenner
    given: Martin
  container-title: Nature Materials
  volume: 11
  URL: 'http://dx.doi.org/10.1038/nmat3283'
  DOI: 10.1038/nmat3283
  issue: 4
  publisher: Nature Publishing Group
  page: 261-263
  type: article-journal
  issued:
    year: 2012
    month: 3
---
